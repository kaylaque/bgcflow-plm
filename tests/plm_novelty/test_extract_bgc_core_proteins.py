"""Tests for extract_bgc_core_proteins.py — run without conda/antismash."""
from pathlib import Path

import pytest
from Bio import SeqIO

from workflow.bgcflow.bgcflow.features.extract_bgc_core_proteins import extract_core_proteins


def test_filters_biosynthetic_only(mini_region_dir, tmp_path):
    """Only the gene_kind=biosynthetic CDS should appear in output."""
    out = tmp_path / "out.faa"
    extract_core_proteins(str(mini_region_dir), "strain1", str(out))
    records = list(SeqIO.parse(str(out), "fasta"))
    assert len(records) == 1
    assert "AAA00001.1" in records[0].id


def test_id_fallback_locus_tag(tmp_path):
    """CDS with only locus_tag (no protein_id) should still be emitted."""
    gbk_content = """\
LOCUS       NODE_LT.region001        100 bp    DNA     linear   BCT 01-JAN-2023
DEFINITION  test locus_tag fallback.
ACCESSION   NODE_LT.region001
VERSION     NODE_LT.region001
FEATURES             Location/Qualifiers
     CDS             1..100
                     /locus_tag="LTAG_0001"
                     /gene_kind="biosynthetic"
                     /translation="MFALLBACK"
ORIGIN
        1 atgtttatca aagagctgat gcgcttaatg tttatcaaag agctgatgcg cttaatgttt
//
"""
    gbk_file = tmp_path / "NODE_LT.region001.gbk"
    gbk_file.write_text(gbk_content)
    out = tmp_path / "out.faa"
    extract_core_proteins(str(tmp_path), "test_strain", str(out))
    records = list(SeqIO.parse(str(out), "fasta"))
    assert len(records) == 1
    assert "LTAG_0001" in records[0].id


def test_id_fallback_gene(tmp_path):
    """CDS with only gene qualifier should be emitted."""
    gbk_content = """\
LOCUS       NODE_GN.region001        100 bp    DNA     linear   BCT 01-JAN-2023
DEFINITION  test gene fallback.
ACCESSION   NODE_GN.region001
VERSION     NODE_GN.region001
FEATURES             Location/Qualifiers
     CDS             1..100
                     /gene="terpS"
                     /gene_kind="biosynthetic"
                     /translation="MGENEFALLBACK"
ORIGIN
        1 atgtttatca aagagctgat gcgcttaatg tttatcaaag agctgatgcg cttaatgttt
//
"""
    gbk_file = tmp_path / "NODE_GN.region001.gbk"
    gbk_file.write_text(gbk_content)
    out = tmp_path / "out.faa"
    extract_core_proteins(str(tmp_path), "test_strain", str(out))
    records = list(SeqIO.parse(str(out), "fasta"))
    assert len(records) == 1
    assert "terpS" in records[0].id


def test_empty_region_writes_empty_faa_exit_0(empty_region_dir, tmp_path):
    """Empty-core region must produce a zero-record FASTA and not raise."""
    out = tmp_path / "empty.faa"
    extract_core_proteins(str(empty_region_dir), "empty_strain", str(out))
    assert out.exists()
    records = list(SeqIO.parse(str(out), "fasta"))
    assert len(records) == 0


def test_malformed_cds_skipped_not_crashed(malformed_region_dir, tmp_path):
    """CDS missing translation qualifier must be skipped; exit 0."""
    out = tmp_path / "malformed.faa"
    extract_core_proteins(str(malformed_region_dir), "mal_strain", str(out))
    assert out.exists()
    records = list(SeqIO.parse(str(out), "fasta"))
    # malformed_region.gbk has 2 biosynthetic CDS: first lacks translation, second has it
    assert len(records) == 1
    assert "CCC00002.1" in records[0].id


def test_fasta_header_format(mini_region_dir, tmp_path):
    """Header must match >{strains}|{region_id}|{protein_id}."""
    out = tmp_path / "out.faa"
    extract_core_proteins(str(mini_region_dir), "myStrain", str(out))
    records = list(SeqIO.parse(str(out), "fasta"))
    assert len(records) == 1
    header = records[0].id
    parts = header.split("|")
    assert len(parts) == 3, f"Expected 3 pipe-separated parts, got: {header}"
    assert parts[0] == "myStrain"
    assert "region" in parts[1].lower() or parts[1] != ""
    assert parts[2] == "AAA00001.1"
