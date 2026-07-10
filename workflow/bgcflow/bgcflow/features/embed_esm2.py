import argparse
import logging
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


def embed_sequences(sequences, model_name="esm2_t30_150M_UR50D", batch_size=8, device="auto"):
    """
    Embed a list of (id, sequence) tuples with ESM2.

    Returns (ids: list[str], vectors: np.ndarray[float32, L2-normalized]).
    Empty input returns ([], np.zeros((0, dim), float32)).
    """
    if not sequences:
        logging.warning("No sequences to embed — returning empty arrays")
        return [], np.zeros((0, 480), dtype=np.float32)  # 480 = esm2_t30_150M_UR50D; update if model changes

    device = _resolve_device(device)
    model, alphabet = load_esm_model(model_name, device)
    batch_converter = alphabet.get_batch_converter()

    import torch  # imported after model load so mock ordering works

    ids_out = []
    vectors_out = []

    for i in range(0, len(sequences), batch_size):
        batch = sequences[i : i + batch_size]
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
            ids_out.append(label)
            vectors_out.append(vec)

    return ids_out, np.stack(vectors_out, axis=0)


def embed_fasta_to_npz(fasta_path, output_path, model_name="esm2_t30_150M_UR50D", batch_size=8, device="auto"):
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
            vectors=np.zeros((0, 480), dtype=np.float32),  # 480 = esm2_t30_150M_UR50D; update if model changes
        )
        return

    ids, vectors = embed_sequences(sequences, model_name=model_name, batch_size=batch_size, device=device)
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(output_path, ids=np.array(ids, dtype=object), vectors=vectors)
    logging.info(f"Saved {len(ids)} embeddings to {output_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Embed FASTA sequences with frozen ESM2")
    parser.add_argument("--fasta", required=True, metavar="FAA", help="Input FASTA file")
    parser.add_argument("--output", required=True, metavar="NPZ", help="Output .npz file")
    parser.add_argument("--model", default="esm2_t30_150M_UR50D", metavar="MODEL_ID")
    parser.add_argument("--batch-size", type=int, default=8, metavar="N")
    parser.add_argument("--device", default="auto", choices=["auto", "cpu", "cuda"])
    args = parser.parse_args()
    embed_fasta_to_npz(args.fasta, args.output, model_name=args.model, batch_size=args.batch_size, device=args.device)
