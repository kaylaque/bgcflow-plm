import argparse
import logging
from pathlib import Path

import numpy as np
import pandas as pd

logging.basicConfig(
    format="%(levelname)-8s %(asctime)s   %(message)s",
    datefmt="%d/%m %H:%M:%S",
    level=logging.DEBUG,
)

_CONFIDENCE_BANDS = [
    (0.0, 0.2, "high"),
    (0.2, 0.4, "medium"),
    (0.4, float("inf"), "low"),
]


def _confidence_band(distance):
    for lo, hi, label in _CONFIDENCE_BANDS:
        if lo <= distance < hi:
            return label
    return "low"


def _parse_header(header):
    """Parse >{strain_id}|{region_id}|{protein_id} header into components."""
    parts = header.split("|", 2)
    if len(parts) == 3:
        return parts[0], parts[1], parts[2]
    return header, "unknown", header


def query_index(embeddings_paths, index_path, meta_path, top_k=10, distance_threshold=0.5, out_per_protein=None, out_per_region=None, expected_model=None):
    """
    Query FAISS index with query embeddings and emit novelty TSVs.

    Per-protein columns: strain_id, region_id, protein_id, top_k_mibig_ids,
                         cosine_distances, nearest_distance, confidence_band, novelty_flag
    Per-region columns: strain_id, region_id, n_core_proteins, frac_with_close_hit,
                        region_novelty_score, flagged_novel
    """
    import faiss

    index = faiss.read_index(str(index_path))
    meta_df = pd.read_parquet(meta_path)
    logging.info(f"Index: {index.ntotal} vectors; metadata: {len(meta_df)} rows")

    if expected_model is not None:
        version_file = Path(index_path).parent / "model_version.txt"
        if version_file.exists():
            indexed_model = None
            for line in version_file.read_text().splitlines():
                if line.startswith("model="):
                    indexed_model = line.split("=", 1)[1].strip()
            if indexed_model and indexed_model != expected_model:
                logging.warning(
                    f"Model mismatch: index was built with '{indexed_model}' "
                    f"but query embeddings use '{expected_model}'. "
                    f"Distances will be meaningless if embedding dims differ."
                )

    per_protein_rows = []

    for npz_path in embeddings_paths:
        data = np.load(npz_path, allow_pickle=True)
        ids = list(data["ids"])
        vectors = data["vectors"]

        if len(ids) == 0:
            logging.warning(f"Empty embeddings file: {npz_path} — skipping")
            continue

        if vectors.shape[1] != index.d:
            raise ValueError(
                f"Embedding dim mismatch: vectors have dim {vectors.shape[1]} "
                f"but FAISS index expects dim {index.d}. "
                f"Ensure the same ESM2 model is used for both index build and query embedding."
            )

        logging.info(f"Querying {len(ids)} proteins from {npz_path}")
        distances_ip, indices = index.search(vectors.astype(np.float32), top_k)
        cosine_distances = 1.0 - distances_ip  # inner product on L2-normalized → cosine sim

        for i, protein_header in enumerate(ids):
            strain_id, region_id, protein_id = _parse_header(protein_header)
            top_k_rows = indices[i]
            top_k_dists = cosine_distances[i]

            top_k_mibig_ids = "|".join(
                meta_df.iloc[r]["mibig_id"] if r >= 0 else "none" for r in top_k_rows
            )
            top_k_dist_str = "|".join(f"{d:.4f}" for d in top_k_dists)
            nearest = float(top_k_dists[0])
            band = _confidence_band(nearest)
            novel = nearest > distance_threshold

            per_protein_rows.append(
                {
                    "strain_id": strain_id,
                    "region_id": region_id,
                    "protein_id": protein_id,
                    "top_k_mibig_ids": top_k_mibig_ids,
                    "cosine_distances": top_k_dist_str,
                    "nearest_distance": nearest,
                    "confidence_band": band,
                    "novelty_flag": novel,
                }
            )

    per_protein_df = pd.DataFrame(
        per_protein_rows,
        columns=[
            "strain_id",
            "region_id",
            "protein_id",
            "top_k_mibig_ids",
            "cosine_distances",
            "nearest_distance",
            "confidence_band",
            "novelty_flag",
        ],
    )

    if out_per_protein:
        Path(out_per_protein).parent.mkdir(parents=True, exist_ok=True)
        per_protein_df.to_csv(out_per_protein, sep="\t", index=False)
        logging.info(f"Per-protein TSV: {len(per_protein_df)} rows → {out_per_protein}")

    # Aggregate per region
    region_rows = []
    for (strain_id, region_id), grp in per_protein_df.groupby(["strain_id", "region_id"]):
        n = len(grp)
        n_close = int((~grp["novelty_flag"]).sum())
        frac = n_close / n if n > 0 else float("nan")
        score = 1.0 - frac if n > 0 else float("nan")
        flagged = score > 0.5 if n > 0 else False
        region_rows.append(
            {
                "strain_id": strain_id,
                "region_id": region_id,
                "n_core_proteins": n,
                "frac_with_close_hit": frac,
                "region_novelty_score": score,
                "flagged_novel": flagged,
            }
        )

    region_df = pd.DataFrame(
        region_rows,
        columns=[
            "strain_id",
            "region_id",
            "n_core_proteins",
            "frac_with_close_hit",
            "region_novelty_score",
            "flagged_novel",
        ],
    )

    if out_per_region:
        Path(out_per_region).parent.mkdir(parents=True, exist_ok=True)
        region_df.to_csv(out_per_region, sep="\t", index=False)
        logging.info(f"Per-region TSV: {len(region_df)} rows → {out_per_region}")

    return per_protein_df, region_df


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Query MIBiG FAISS index with BGC core protein embeddings"
    )
    parser.add_argument(
        "--embeddings",
        nargs="+",
        required=True,
        metavar="NPZ",
        help="One or more .npz embedding files (each from embed_esm2.py)",
    )
    parser.add_argument("--index", required=True, metavar="FAISS", help="FAISS index file")
    parser.add_argument("--meta", required=True, metavar="PARQUET", help="metadata.parquet")
    parser.add_argument("--top-k", type=int, default=10, metavar="N")
    parser.add_argument("--distance-threshold", type=float, default=0.5, metavar="F")
    parser.add_argument("--out-per-protein", required=True, metavar="TSV")
    parser.add_argument("--out-per-region", required=True, metavar="TSV")
    parser.add_argument(
        "--expected-model",
        default=None,
        metavar="MODEL_ID",
        help="Model used to embed query sequences; checked against model_version.txt for mismatch warning",
    )
    args = parser.parse_args()
    query_index(
        args.embeddings,
        args.index,
        args.meta,
        top_k=args.top_k,
        distance_threshold=args.distance_threshold,
        out_per_protein=args.out_per_protein,
        out_per_region=args.out_per_region,
        expected_model=args.expected_model,
    )
