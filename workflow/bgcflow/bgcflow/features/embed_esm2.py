import argparse
import json
import logging
import os
from pathlib import Path

import numpy as np

logging.basicConfig(
    format="%(levelname)-8s %(asctime)s   %(message)s",
    datefmt="%d/%m %H:%M:%S",
    level=logging.DEBUG,
)


def load_esm_model(model_name, device):
    """Load ESM2 model and alphabet. Isolated here so tests can mock it."""
    import esm  # imported lazily — not available in test env without torch

    model, alphabet = esm.pretrained.load_model_and_alphabet(model_name)
    model.eval()
    model = model.to(device)
    return model, alphabet


def _resolve_device(device_str):
    if device_str == "auto":
        try:
            import torch

            return "cuda" if torch.cuda.is_available() else "cpu"
        except ImportError:
            return "cpu"
    return device_str


def _load_checkpoint(checkpoint_dir, model_name):
    """
    Load a checkpoint dir's already-embedded {id: vector} map, discarding it
    (and any manifest) if it was built with a different model — embeddings
    from two different models must never mix in one resumed run.
    """
    checkpoint_dir = Path(checkpoint_dir)
    manifest_path = checkpoint_dir / "manifest.json"
    shard_paths = sorted(checkpoint_dir.glob("shard_*.npz"))

    stale = False
    if manifest_path.exists():
        manifest = json.loads(manifest_path.read_text())
        if manifest.get("model_name") != model_name:
            stale = True
    elif shard_paths:
        stale = True  # shards with no manifest — can't trust their provenance

    if stale:
        logging.warning(
            f"Checkpoint in {checkpoint_dir} was built with a different model "
            "— discarding and starting fresh"
        )
        for shard_path in shard_paths:
            shard_path.unlink()
        shard_paths = []

    checkpoint_dir.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text(json.dumps({"model_name": model_name}))

    done_map = {}
    for shard_path in shard_paths:
        data = np.load(shard_path, allow_pickle=True)
        for label, vec in zip(data["ids"], data["vectors"]):
            done_map[str(label)] = vec

    return done_map, len(shard_paths)


def _write_checkpoint_shard(checkpoint_dir, shard_index, ids, vectors):
    """Atomically write one checkpoint shard (temp file + rename) so a killed
    process never leaves a half-written, corrupt shard behind."""
    checkpoint_dir = Path(checkpoint_dir)
    final_path = checkpoint_dir / f"shard_{shard_index:06d}.npz"
    # np.savez_compressed appends ".npz" if the name doesn't already end with
    # it, so the temp name must end in ".npz" too or the later rename 404s.
    tmp_path = checkpoint_dir / f".tmp_shard_{shard_index:06d}.npz"
    np.savez_compressed(
        tmp_path,
        ids=np.array(ids, dtype=object),
        vectors=np.stack(vectors).astype(np.float32),
    )
    os.replace(tmp_path, final_path)


def embed_sequences(
    sequences,
    model_name="esm2_t30_150M_UR50D",
    batch_size=8,
    device="auto",
    checkpoint_dir=None,
    checkpoint_every=50,
):
    """
    Embed a list of (id, sequence) tuples with ESM2.

    If checkpoint_dir is given, already-embedded ids (from a prior, possibly
    interrupted run) are loaded and skipped, and progress is written to disk
    every `checkpoint_every` batches — so a killed/interrupted run can be
    resumed by calling again with the same checkpoint_dir and model_name.

    Returns (ids: list[str], vectors: np.ndarray[float32, L2-normalized]),
    always in the same order as the input `sequences`, regardless of resume.
    Empty input returns ([], np.zeros((0, dim), float32)).
    """
    if not sequences:
        logging.warning("No sequences to embed — returning empty arrays")
        return [], np.zeros(
            (0, 480), dtype=np.float32
        )  # 480 = esm2_t30_150M_UR50D; update if model changes

    done_map = {}
    next_shard_index = 0
    if checkpoint_dir is not None:
        done_map, next_shard_index = _load_checkpoint(checkpoint_dir, model_name)

    remaining = [(label, seq) for label, seq in sequences if label not in done_map]

    if remaining:
        if done_map:
            logging.info(
                f"Resuming: {len(done_map)} already embedded, {len(remaining)} remaining"
            )

        device = _resolve_device(device)
        model, alphabet = load_esm_model(model_name, device)
        batch_converter = alphabet.get_batch_converter()

        import torch  # imported after model load so mock ordering works

        pending_ids = []
        pending_vectors = []
        batches_since_checkpoint = 0

        for i in range(0, len(remaining), batch_size):
            batch = remaining[i : i + batch_size]
            batch_labels, _, batch_tokens = batch_converter(batch)
            batch_tokens = batch_tokens.to(device)

            with torch.no_grad():
                results = model(batch_tokens, repr_layers=[model.num_layers])

            token_repr = results["representations"][model.num_layers]

            for j, (label, _) in enumerate(batch):
                # mean-pool over non-padding, non-special tokens (positions 1 to -1)
                vec = token_repr[j, 1:-1].mean(0).cpu().numpy().astype(np.float32)
                norm = np.linalg.norm(vec)
                if norm > 0:
                    vec = vec / norm
                done_map[label] = vec
                pending_ids.append(label)
                pending_vectors.append(vec)

            batches_since_checkpoint += 1
            is_last_batch = i + batch_size >= len(remaining)
            if checkpoint_dir is not None and (
                batches_since_checkpoint >= checkpoint_every or is_last_batch
            ):
                if pending_ids:
                    _write_checkpoint_shard(
                        checkpoint_dir, next_shard_index, pending_ids, pending_vectors
                    )
                    next_shard_index += 1
                pending_ids = []
                pending_vectors = []
                batches_since_checkpoint = 0

    ids_out = [label for label, _ in sequences]
    vectors_out = np.stack([done_map[label] for label in ids_out], axis=0)
    return ids_out, vectors_out


def embed_fasta_to_npz(
    fasta_path,
    output_path,
    model_name="esm2_t30_150M_UR50D",
    batch_size=8,
    device="auto",
):
    """Read a FASTA file, embed all sequences, save as .npz."""
    from Bio import SeqIO

    sequences = []
    with open(fasta_path) as fh:
        for rec in SeqIO.parse(fh, "fasta"):
            sequences.append((rec.id, str(rec.seq)))

    if not sequences:
        logging.warning(f"Empty FASTA {fasta_path} — writing empty .npz")
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(
            output_path,
            ids=np.array([], dtype=object),
            vectors=np.zeros(
                (0, 480), dtype=np.float32
            ),  # 480 = esm2_t30_150M_UR50D; update if model changes
        )
        return

    ids, vectors = embed_sequences(
        sequences, model_name=model_name, batch_size=batch_size, device=device
    )
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(output_path, ids=np.array(ids, dtype=object), vectors=vectors)
    logging.info(f"Saved {len(ids)} embeddings to {output_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Embed FASTA sequences with frozen ESM2"
    )
    parser.add_argument(
        "--fasta", required=True, metavar="FAA", help="Input FASTA file"
    )
    parser.add_argument(
        "--output", required=True, metavar="NPZ", help="Output .npz file"
    )
    parser.add_argument("--model", default="esm2_t30_150M_UR50D", metavar="MODEL_ID")
    parser.add_argument("--batch-size", type=int, default=8, metavar="N")
    parser.add_argument("--device", default="auto", choices=["auto", "cpu", "cuda"])
    args = parser.parse_args()
    embed_fasta_to_npz(
        args.fasta,
        args.output,
        model_name=args.model,
        batch_size=args.batch_size,
        device=args.device,
    )
