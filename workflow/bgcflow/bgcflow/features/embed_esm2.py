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


def _length_bucketed_batches(sequences, batch_size, max_tokens_per_batch):
    """
    Group (id, seq) pairs into batches respecting both batch_size and a
    token budget (batch_size_so_far * longest_seq_in_batch).

    ESM pads every sequence in a batch to the batch's longest member, and
    self-attention memory scales with the square of that padded length —
    protein datasets like MIBiG mix ~100-600aa domains with multi-thousand-
    residue megasynthases (PKS/NRPS), so batching by count alone can put a
    5000-residue outlier in the same batch as several short sequences and
    force all of them through a catastrophically large forward pass. A
    single sequence that alone exceeds the budget still gets its own batch
    — it's processed, never silently dropped.
    """
    ordered = sorted(sequences, key=lambda pair: len(pair[1]))
    batches = []
    current = []
    current_max_len = 0
    for item in ordered:
        item_len = len(item[1])
        candidate_max_len = max(current_max_len, item_len)
        candidate_size = len(current) + 1
        over_budget = candidate_size * candidate_max_len > max_tokens_per_batch
        if current and (candidate_size > batch_size or over_budget):
            batches.append(current)
            current = [item]
            current_max_len = item_len
        else:
            current.append(item)
            current_max_len = candidate_max_len
    if current:
        batches.append(current)
    return batches


def embed_sequences(
    sequences,
    model_name="esm2_t30_150M_UR50D",
    batch_size=8,
    device="auto",
    checkpoint_dir=None,
    checkpoint_every=50,
    max_tokens_per_batch=None,
):
    """
    Embed a list of (id, sequence) tuples with ESM2.

    Sequences are grouped into batches by length (see
    _length_bucketed_batches) so a handful of very long outliers don't force
    huge padding onto batches of otherwise-short sequences. max_tokens_per_batch
    defaults to batch_size * 1024 when not given.

    If checkpoint_dir is given, already-embedded ids (from a prior, possibly
    interrupted run) are loaded and skipped, and progress is written to disk
    every `checkpoint_every` batches — so a killed/interrupted run can be
    resumed by calling again with the same checkpoint_dir and model_name.

    Returns (ids: list[str], vectors: np.ndarray[float32, L2-normalized]),
    always in the same order as the input `sequences`, regardless of resume.
    Empty input returns ([], np.zeros((0, dim), float32)).
    """
    if max_tokens_per_batch is None:
        max_tokens_per_batch = batch_size * 1024
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

        batches = _length_bucketed_batches(remaining, batch_size, max_tokens_per_batch)
        for batch_index, batch in enumerate(batches):
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
            is_last_batch = batch_index == len(batches) - 1
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
    max_tokens_per_batch=None,
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
        sequences,
        model_name=model_name,
        batch_size=batch_size,
        device=device,
        max_tokens_per_batch=max_tokens_per_batch,
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
    parser.add_argument(
        "--max-tokens-per-batch",
        type=int,
        default=None,
        metavar="N",
        help="Token budget per batch (default: batch_size * 1024) — caps padded "
        "batch size so long outlier sequences don't blow up memory",
    )
    args = parser.parse_args()
    embed_fasta_to_npz(
        args.fasta,
        args.output,
        model_name=args.model,
        batch_size=args.batch_size,
        device=args.device,
        max_tokens_per_batch=args.max_tokens_per_batch,
    )
