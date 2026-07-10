import argparse
import logging
from pathlib import Path

from Bio import SeqIO

logging.basicConfig(
    format="%(levelname)-8s %(asctime)s   %(message)s",
    datefmt="%d/%m %H:%M:%S",
    level=logging.DEBUG,
)


def extract_core_proteins(input_dir, strains, output_file):
    """
    Extract biosynthetic-core CDS translations from antiSMASH region GBKs.

    Scans all *.region*.gbk files in input_dir, keeps only CDS features where
    gene_kind == "biosynthetic", and writes a FASTA with headers
    >{strains}|{region_id}|{protein_id}.  Emits an empty FASTA (exit 0) when
    no core proteins are found so the DAG does not fail on unusual strains.
    """
    input_dir = Path(input_dir)
    output_file = Path(output_file)
    output_file.parent.mkdir(parents=True, exist_ok=True)

    records_written = 0

    with open(output_file, "w") as out_fh:
        for gbk_path in sorted(input_dir.glob("*.region*.gbk")):
            region_id = gbk_path.stem
            logging.info(f"Processing {gbk_path}")
            with open(gbk_path) as fh:
                for record in SeqIO.parse(fh, "genbank"):
                    for feature in record.features:
                        if feature.type != "CDS":
                            continue
                        gene_kind = feature.qualifiers.get("gene_kind", [""])[0]
                        if gene_kind != "biosynthetic":
                            continue

                        protein_id = _resolve_protein_id(feature)
                        if protein_id is None:
                            logging.warning(
                                f"CDS in {gbk_path} has no protein_id/locus_tag/gene — skipping"
                            )
                            continue

                        translation = feature.qualifiers.get("translation", [None])[0]
                        if translation is None:
                            logging.warning(
                                f"CDS {protein_id} in {gbk_path} missing translation qualifier — skipping"
                            )
                            continue

                        header = f">{strains}|{region_id}|{protein_id}"
                        out_fh.write(f"{header}\n{translation}\n")
                        records_written += 1

    if records_written == 0:
        logging.warning(
            f"No biosynthetic-core proteins found in {input_dir} for strain {strains}"
        )
    else:
        logging.info(f"Wrote {records_written} core protein records to {output_file}")


def _resolve_protein_id(feature):
    for key in ("protein_id", "locus_tag", "gene"):
        if key in feature.qualifiers:
            return feature.qualifiers[key][0]
    return None


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extract biosynthetic-core proteins from antiSMASH region GBKs"
    )
    parser.add_argument(
        "--input",
        required=True,
        metavar="DIR",
        help="antiSMASH strain output directory (contains *.region*.gbk files)",
    )
    parser.add_argument(
        "--strains",
        required=True,
        metavar="STRAIN_ID",
        help="Strain identifier embedded in FASTA headers",
    )
    parser.add_argument(
        "--output",
        required=True,
        metavar="FAA",
        help="Output FASTA file path",
    )
    args = parser.parse_args()
    extract_core_proteins(args.input, args.strains, args.output)
