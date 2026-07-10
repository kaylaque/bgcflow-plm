import argparse
import hashlib
import logging
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import SeqIO

logging.basicConfig(
    format="%(levelname)-8s %(asctime)s   %(message)s",
    datefmt="%d/%m %H:%M:%S",
    level=logging.DEBUG,
)


def _resolve_protein_id(feature):
    for key in ("protein_id", "locus_tag", "gene"):
        if key in feature.qualifiers:
            return feature.qualifiers[key][0]
    return None


def extract_core_proteins_from_dir(mibig_dir):
    """
    Walk mibig_dir for *.gbk files, extract biosynthetic-core CDS.

    Returns list of dicts: {mibig_id, protein_id, product_class, sequence, sequence_hash}
    """
    mibig_dir = Path(mibig_dir)
    proteins = []
    seen_ids = set()

    for gbk_path in sorted(mibig_dir.rglob("*.gbk")):
        mibig_id = gbk_path.stem
        logging.info(f"Scanning {gbk_path}")
        with open(gbk_path) as fh:
            for record in SeqIO.parse(fh, "genbank"):
                product_class = "unknown"
                for feat in record.features:
                    if feat.type == "cluster" or feat.type == "region":
                        product_class = feat.qualifiers.get("product", ["unknown"])[0]
                        break

                for feature in record.features:
                    if feature.type != "CDS":
                        continue
                    gene_kind = feature.qualifiers.get("gene_kind", [""])[0]
                    if gene_kind != "biosynthetic":
                        continue

                    protein_id = _resolve_protein_id(feature)
                    if protein_id is None:
                        logging.warning(f"CDS in {gbk_path} has no identifier — skipping")
                        continue

                    translation = feature.qualifiers.get("translation", [None])[0]
                    if translation is None:
                        logging.warning(f"CDS {protein_id} missing translation — skipping")
                        continue

                    unique_key = f"{mibig_id}::{protein_id}"
                    if unique_key in seen_ids:
                        continue
                    seen_ids.add(unique_key)

                    seq_hash = hashlib.md5(translation.encode()).hexdigest()
                    proteins.append(
                        {
                            "mibig_id": mibig_id,
                            "protein_id": protein_id,
                            "product_class": product_class,
                            "sequence": translation,
                            "sequence_hash": seq_hash,
                        }
                    )

    return proteins


def build_index(mibig_dir, output_dir, model_name="esm2_t30_150M_UR50D", batch_size=8, device="auto"):
    """
    Build FAISS index from MIBiG GBK files and save to output_dir.

    Writes:
      - index.faiss
      - metadata.parquet  (row order == FAISS row order)
      - model_version.txt
    """
    import sys as _sys
    from pathlib import Path as _Path
    _sys.path.insert(0, str(_Path(__file__).resolve().parent))

    import faiss
    from embed_esm2 import embed_sequences  # noqa: E402

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    proteins = extract_core_proteins_from_dir(mibig_dir)
    if not proteins:
        raise ValueError(f"No biosynthetic-core proteins found in {mibig_dir}")

    logging.info(f"Found {len(proteins)} core proteins across MIBiG GBKs")

    sequences = [(p["mibig_id"] + "::" + p["protein_id"], p["sequence"]) for p in proteins]
    ids, vectors = embed_sequences(sequences, model_name=model_name, batch_size=batch_size, device=device)

    dim = vectors.shape[1]
    index = faiss.IndexFlatIP(dim)
    index.add(vectors)
    faiss.write_index(index, str(output_dir / "index.faiss"))
    logging.info(f"FAISS index built: {index.ntotal} vectors, dim={dim}")

    meta_df = pd.DataFrame(
        {
            "mibig_id": [p["mibig_id"] for p in proteins],
            "protein_id": [p["protein_id"] for p in proteins],
            "product_class": [p["product_class"] for p in proteins],
            "sequence_hash": [p["sequence_hash"] for p in proteins],
        }
    )
    meta_df.to_parquet(output_dir / "metadata.parquet", index=False)
    logging.info(f"Metadata written: {len(meta_df)} rows")

    version_txt = output_dir / "model_version.txt"
    try:
        import esm as _esm

        esm_ver = getattr(_esm, "__version__", "unknown")
    except ImportError:
        esm_ver = "unknown"
    version_txt.write_text(f"model={model_name}\nfair-esm={esm_ver}\n")
    logging.info(f"Model version: {model_name} (fair-esm {esm_ver})")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Build ESM2 FAISS index from MIBiG GBK files"
    )
    parser.add_argument(
        "--mibig-dir",
        required=True,
        metavar="DIR",
        help="Directory containing MIBiG GBK files (recursively scanned)",
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        metavar="DIR",
        help="Output directory for index.faiss, metadata.parquet, model_version.txt",
    )
    parser.add_argument("--model", default="esm2_t30_150M_UR50D", metavar="MODEL_ID")
    parser.add_argument("--batch-size", type=int, default=8, metavar="N")
    parser.add_argument("--device", default="auto", choices=["auto", "cpu", "cuda"])
    args = parser.parse_args()
    build_index(
        args.mibig_dir,
        args.output_dir,
        model_name=args.model,
        batch_size=args.batch_size,
        device=args.device,
    )
