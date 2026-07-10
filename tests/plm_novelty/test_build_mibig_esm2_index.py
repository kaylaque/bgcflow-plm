"""Tests for build_mibig_esm2_index.py — ESM2 model is mocked."""
import sys
from pathlib import Path
from unittest.mock import MagicMock, patch

import numpy as np
import pandas as pd
import pytest

DIM = 8


def _make_stub_embeddings():
    """Patch embed_sequences to return deterministic vectors."""

    def fake_embed_sequences(sequences, model_name="", batch_size=8, device="auto"):
        n = len(sequences)
        ids = [s[0] for s in sequences]
        vecs = np.random.default_rng(42).random((n, DIM)).astype(np.float32)
        norms = np.linalg.norm(vecs, axis=1, keepdims=True)
        vecs = vecs / norms
        return ids, vecs

    return fake_embed_sequences


@pytest.fixture(autouse=True)
def mock_embed():
    with patch(
        "workflow.bgcflow.bgcflow.features.build_mibig_esm2_index.embed_sequences",
        side_effect=_make_stub_embeddings(),
    ):
        yield


def test_index_and_metadata_row_count_match(mini_mibig_dir, tmp_path):
    """FAISS index ntotal must equal metadata.parquet row count (row-order invariant)."""
    import faiss

    from workflow.bgcflow.bgcflow.features.build_mibig_esm2_index import build_index

    out_dir = tmp_path / "index_out"
    build_index(str(mini_mibig_dir), str(out_dir), model_name="stub")

    index = faiss.read_index(str(out_dir / "index.faiss"))
    meta = pd.read_parquet(out_dir / "metadata.parquet")
    assert index.ntotal == len(meta), (
        f"Row-order invariant violated: index has {index.ntotal} vectors "
        f"but metadata has {len(meta)} rows"
    )
    # mini_mibig has 2 GBKs with 2 biosynthetic core CDS each → 4 total
    assert index.ntotal == 4


def test_model_version_txt_contains_model_id(mini_mibig_dir, tmp_path):
    """model_version.txt must include the model name passed in."""
    from workflow.bgcflow.bgcflow.features.build_mibig_esm2_index import build_index

    out_dir = tmp_path / "index_out2"
    build_index(str(mini_mibig_dir), str(out_dir), model_name="esm2_t12_35M_UR50D")
    version_txt = (out_dir / "model_version.txt").read_text()
    assert "esm2_t12_35M_UR50D" in version_txt


def test_metadata_columns(mini_mibig_dir, tmp_path):
    """metadata.parquet must have the required columns."""
    from workflow.bgcflow.bgcflow.features.build_mibig_esm2_index import build_index

    out_dir = tmp_path / "index_out3"
    build_index(str(mini_mibig_dir), str(out_dir))
    meta = pd.read_parquet(out_dir / "metadata.parquet")
    for col in ("mibig_id", "protein_id", "product_class", "sequence_hash"):
        assert col in meta.columns, f"Missing column: {col}"
