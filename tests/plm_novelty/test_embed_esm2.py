"""Tests for embed_esm2.py — ESM2 model is mocked (no torch/fair-esm in CI)."""
import sys
from pathlib import Path
from unittest.mock import MagicMock, patch

import numpy as np
import pytest

DIM = 8  # stub embedding dimension


def _make_stub_model(dim=DIM):
    """Return a mock (model, alphabet) pair emitting deterministic dim-vectors."""
    import types

    def batch_converter(batch):
        labels = [b[0] for b in batch]
        seqs = [b[1] for b in batch]
        # fake tokens: shape (batch, len+2) filled with zeros
        max_len = max(len(s) for s in seqs) + 2
        tokens = np.zeros((len(batch), max_len), dtype=np.int64)
        return labels, seqs, tokens

    alphabet_mock = MagicMock()
    alphabet_mock.get_batch_converter.return_value = batch_converter

    class FakeResult:
        def __init__(self, n, length, dim):
            arr = np.ones((n, length, dim), dtype=np.float32) * 0.5
            import torch
            self.representations = {1: torch.tensor(arr)}

    class FakeModel:
        num_layers = 1

        def eval(self):
            return self

        def to(self, device):
            return self

        def __call__(self, tokens, repr_layers=None):
            n, seq_len = tokens.shape
            return FakeResult(n, seq_len, dim)

    return FakeModel(), alphabet_mock


@pytest.fixture(autouse=True)
def mock_esm_load():
    """Patch load_esm_model in embed_esm2 before each test."""
    model, alphabet = _make_stub_model()
    with patch(
        "workflow.bgcflow.bgcflow.features.embed_esm2.load_esm_model",
        return_value=(model, alphabet),
    ):
        yield


def test_npz_keys_dtype_normalized(tmp_path):
    """Output .npz has keys 'ids' and 'vectors', float32, L2-normalized."""
    import torch

    fasta = tmp_path / "test.faa"
    fasta.write_text(">seq1\nMFIXEDSEQ\n>seq2\nMANOTHERONE\n")
    out = tmp_path / "out.npz"

    from workflow.bgcflow.bgcflow.features.embed_esm2 import embed_fasta_to_npz

    with patch("workflow.bgcflow.bgcflow.features.embed_esm2.torch", torch):
        embed_fasta_to_npz(str(fasta), str(out), model_name="stub", batch_size=8)

    data = np.load(str(out), allow_pickle=True)
    assert set(data.files) >= {"ids", "vectors"}
    assert data["vectors"].dtype == np.float32
    norms = np.linalg.norm(data["vectors"], axis=1)
    np.testing.assert_allclose(norms, 1.0, atol=1e-5)


def test_empty_fasta_writes_empty_npz(tmp_path):
    """Empty FASTA → empty .npz arrays, no exception."""
    fasta = tmp_path / "empty.faa"
    fasta.write_text("")
    out = tmp_path / "empty_out.npz"

    from workflow.bgcflow.bgcflow.features.embed_esm2 import embed_fasta_to_npz

    embed_fasta_to_npz(str(fasta), str(out))
    data = np.load(str(out), allow_pickle=True)
    assert len(data["ids"]) == 0
    assert data["vectors"].shape[0] == 0


def test_batch_larger_than_record_count(tmp_path):
    """batch_size > len(sequences) must not crash."""
    import torch

    fasta = tmp_path / "one.faa"
    fasta.write_text(">single\nMONLYONE\n")
    out = tmp_path / "one.npz"

    from workflow.bgcflow.bgcflow.features.embed_esm2 import embed_fasta_to_npz

    with patch("workflow.bgcflow.bgcflow.features.embed_esm2.torch", torch):
        embed_fasta_to_npz(str(fasta), str(out), batch_size=100)

    data = np.load(str(out), allow_pickle=True)
    assert len(data["ids"]) == 1
