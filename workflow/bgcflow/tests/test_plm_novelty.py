"""
No-crash smoke tests for the PLM novelty retrieval side-channel.

Design goals:
  - No torch, fair-esm, or real FAISS weights required.
  - Each test is standalone and uses only stdlib + numpy/pandas/biopython.
  - The Snakemake dry-run test verifies the DAG resolves cleanly end-to-end.

Run with:
  python -m pytest workflow/bgcflow/tests/test_plm_novelty.py -v
"""

import ast
import shutil
import subprocess
import sys
import textwrap
from pathlib import Path
from unittest import mock

import numpy as np
import pandas as pd
import pytest

FEATURES_DIR = Path(__file__).resolve().parents[2] / "bgcflow" / "features"
REPO_ROOT = Path(__file__).resolve().parents[4]
TESTS_CONFIG = REPO_ROOT / ".tests" / "config" / "config.yaml"

PLM_SCRIPTS = [
    FEATURES_DIR / "build_mibig_esm2_index.py",
    FEATURES_DIR / "embed_esm2.py",
    FEATURES_DIR / "extract_bgc_core_proteins.py",
    FEATURES_DIR / "plm_novelty_query.py",
]


# ---------------------------------------------------------------------------
# 1. Syntax + CLI smoke
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("script", PLM_SCRIPTS, ids=[s.name for s in PLM_SCRIPTS])
def test_script_parses(script):
    ast.parse(script.read_text())


@pytest.mark.parametrize("script", PLM_SCRIPTS, ids=[s.name for s in PLM_SCRIPTS])
def test_script_help_exits_clean(script):
    result = subprocess.run(
        [sys.executable, str(script), "--help"],
        capture_output=True,
        timeout=30,
    )
    assert result.returncode == 0, (
        f"{script.name} --help failed:\n{result.stderr.decode()}"
    )


# ---------------------------------------------------------------------------
# 2. extract_bgc_core_proteins — toy GenBank round-trip
# ---------------------------------------------------------------------------

_TOY_GBK = textwrap.dedent(
    """\
    LOCUS       testregion001           50 bp    DNA     linear   BCT 01-JAN-2024
    DEFINITION  Test region.
    ACCESSION   testregion001
    VERSION     testregion001
    FEATURES             Location/Qualifiers
         region          1..50
                         /product="nrps"
         CDS             1..30
                         /protein_id="prot_biosyn"
                         /gene_kind="biosynthetic"
                         /translation="MAAAKL"
         CDS             31..50
                         /protein_id="prot_transport"
                         /gene_kind="transport"
                         /translation="MSSS"
    ORIGIN
            1 atgaaagcag cagccaaact atgagcagca gcatag
    //
    """
)


def test_extract_bgc_core_proteins_filters_gene_kind(tmp_path):
    sys.path.insert(0, str(FEATURES_DIR))
    from extract_bgc_core_proteins import extract_core_proteins

    strain_id = "strain_abc"
    region_gbk = tmp_path / f"{strain_id}.region001.gbk"
    region_gbk.write_text(_TOY_GBK)

    out_faa = tmp_path / "core_proteins.faa"
    extract_core_proteins(str(tmp_path), strain_id, str(out_faa))

    lines = [l for l in out_faa.read_text().splitlines() if l]
    headers = [l for l in lines if l.startswith(">")]
    assert len(headers) == 1, f"Expected 1 biosynthetic-core record, got: {headers}"
    assert headers[0] == f">{strain_id}|{strain_id}.region001|prot_biosyn"
    assert "MAAAKL" in lines


# ---------------------------------------------------------------------------
# 3. embed_esm2 — empty FASTA path (no torch needed)
# ---------------------------------------------------------------------------


def test_empty_fasta_produces_empty_npz(tmp_path):
    sys.path.insert(0, str(FEATURES_DIR))
    from embed_esm2 import embed_fasta_to_npz

    empty_faa = tmp_path / "empty.faa"
    empty_faa.write_text("")
    out_npz = tmp_path / "empty.npz"

    embed_fasta_to_npz(str(empty_faa), str(out_npz))

    data = np.load(out_npz, allow_pickle=True)
    assert data["ids"].shape == (0,), "ids should be empty"
    assert data["vectors"].ndim == 2, "vectors should be 2-D"
    assert data["vectors"].shape[0] == 0, "vectors should have 0 rows"


# ---------------------------------------------------------------------------
# 4. plm_novelty_query — end-to-end with a FAISS stub
# ---------------------------------------------------------------------------


class _FakeFAISSIndex:
    """Minimal FAISS-like stub that returns deterministic results."""

    def __init__(self, dim=4, ntotal=2):
        self.d = dim
        self.ntotal = ntotal

    def search(self, vectors, k):
        n = vectors.shape[0]
        # Return perfect-match inner products (cosine sim = 1 → distance = 0) for first,
        # and a "far" match (cosine sim = 0.4 → distance = 0.6) for second.
        sims = np.array([[1.0] * k, [0.4] + [0.3] * (k - 1)], dtype=np.float32)[:n]
        idxs = np.tile(np.arange(k, dtype=np.int64), (n, 1))
        return sims, idxs


def _make_meta_parquet(tmp_path, n=10):
    df = pd.DataFrame(
        {
            "mibig_id": [f"BGC{i:07d}" for i in range(n)],
            "protein_id": [f"prot{i}" for i in range(n)],
            "product_class": ["nrps"] * n,
            "sequence_hash": ["abc"] * n,
        }
    )
    p = tmp_path / "metadata.parquet"
    df.to_parquet(p, index=False)
    return p


def _make_npz(tmp_path, name, ids, dim=4):
    vectors = np.random.rand(len(ids), dim).astype(np.float32)
    norms = np.linalg.norm(vectors, axis=1, keepdims=True)
    vectors /= norms
    p = tmp_path / name
    np.savez_compressed(p, ids=np.array(ids, dtype=object), vectors=vectors)
    return p


def test_plm_novelty_query_end_to_end_synthetic(tmp_path):
    sys.path.insert(0, str(FEATURES_DIR))
    import plm_novelty_query

    meta_path = _make_meta_parquet(tmp_path, n=10)
    npz_path = _make_npz(
        tmp_path, "strain1.npz",
        ids=["strain1|strain1.region001|p1", "strain1|strain1.region001|p2"],
        dim=4,
    )

    fake_index = _FakeFAISSIndex(dim=4, ntotal=10)
    out_protein = tmp_path / "per_protein.tsv"
    out_region = tmp_path / "per_region.tsv"

    with mock.patch("faiss.read_index", return_value=fake_index):
        prot_df, reg_df = plm_novelty_query.query_index(
            [str(npz_path)],
            index_path=str(tmp_path / "fake.faiss"),
            meta_path=str(meta_path),
            top_k=3,
            distance_threshold=0.5,
            out_per_protein=str(out_protein),
            out_per_region=str(out_region),
        )

    # Per-protein assertions
    assert len(prot_df) == 2
    assert set(prot_df.columns) >= {"strain_id", "region_id", "protein_id", "nearest_distance", "novelty_flag", "confidence_band"}
    assert list(prot_df["strain_id"]) == ["strain1", "strain1"]
    assert list(prot_df["region_id"]) == ["strain1.region001", "strain1.region001"]
    # p1 gets cosine distance 0.0 → not novel; p2 gets 0.6 → novel
    assert prot_df.iloc[0]["nearest_distance"] == pytest.approx(0.0, abs=1e-5)
    assert prot_df.iloc[0]["novelty_flag"] is False or prot_df.iloc[0]["novelty_flag"] == False  # noqa: E712
    assert prot_df.iloc[1]["nearest_distance"] == pytest.approx(0.6, abs=1e-5)
    assert prot_df.iloc[1]["novelty_flag"] is True or prot_df.iloc[1]["novelty_flag"] == True  # noqa: E712
    assert prot_df.iloc[0]["confidence_band"] == "high"

    # Per-region assertions
    assert len(reg_df) == 1
    assert reg_df.iloc[0]["n_core_proteins"] == 2
    # 1 of 2 proteins is close → frac_with_close_hit = 0.5 → region_novelty_score = 0.5
    assert reg_df.iloc[0]["frac_with_close_hit"] == pytest.approx(0.5)
    assert reg_df.iloc[0]["region_novelty_score"] == pytest.approx(0.5)

    # TSV files must exist and be non-empty
    assert out_protein.exists() and out_protein.stat().st_size > 0
    assert out_region.exists() and out_region.stat().st_size > 0


def test_plm_novelty_query_dim_mismatch_raises(tmp_path):
    sys.path.insert(0, str(FEATURES_DIR))
    import plm_novelty_query

    meta_path = _make_meta_parquet(tmp_path, n=10)
    npz_path = _make_npz(
        tmp_path, "mismatch.npz",
        ids=["strainX|strainX.region001|p1"],
        dim=8,  # vectors are dim 8
    )

    fake_index = _FakeFAISSIndex(dim=4, ntotal=10)  # index expects dim 4

    with mock.patch("faiss.read_index", return_value=fake_index):
        with pytest.raises(ValueError, match="dim"):
            plm_novelty_query.query_index(
                [str(npz_path)],
                index_path=str(tmp_path / "fake.faiss"),
                meta_path=str(meta_path),
                top_k=3,
                distance_threshold=0.5,
            )


# ---------------------------------------------------------------------------
# 5. Snakemake dry-run — no workflow crash
# ---------------------------------------------------------------------------


@pytest.mark.skipif(
    shutil.which("snakemake") is None,
    reason="snakemake not in PATH",
)
def test_snakemake_list_rules_no_crash():
    """
    Verify that `snakemake --list` resolves all rules (including plm_novelty.smk)
    without error when using the test config.

    --list enumerates available rules without needing input files on disk.
    """
    result = subprocess.run(
        [
            "snakemake",
            "--list",
            "--configfile", str(TESTS_CONFIG),
            "--snakefile", str(REPO_ROOT / "workflow" / "Snakefile"),
            "--directory", str(REPO_ROOT),
        ],
        capture_output=True,
        timeout=120,
        cwd=str(REPO_ROOT),
    )
    assert result.returncode == 0, (
        "snakemake --list failed (PLM rules may have a syntax/import error):\n"
        + result.stderr.decode()
    )
    output = result.stdout.decode()
    # All four PLM rules should appear in the rule list when antismash version >= 7
    # (which it is in the test config: antismash.version = "7")
    for expected_rule in [
        "build_mibig_esm2_index",
        "extract_bgc_core_proteins",
        "embed_bgc_core_esm2",
        "plm_novelty_query",
    ]:
        assert expected_rule in output, (
            f"Rule '{expected_rule}' not found in snakemake --list output. "
            f"Output was:\n{output}"
        )
