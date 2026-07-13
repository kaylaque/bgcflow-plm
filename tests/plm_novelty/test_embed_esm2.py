"""Tests for embed_esm2.py — ESM2 model is mocked (no torch/fair-esm in CI)."""

import json
import sys
from pathlib import Path
from unittest.mock import MagicMock, patch

import numpy as np
import pytest

DIM = 8  # stub embedding dimension


def _make_stub_model(dim=DIM, call_log=None):
    """Return a mock (model, alphabet) pair emitting deterministic dim-vectors.

    If call_log is given, each __call__ appends the batch's labels to it —
    lets tests assert exactly which sequences were (re-)embedded.
    """
    import types

    def batch_converter(batch):
        import torch

        labels = [b[0] for b in batch]
        seqs = [b[1] for b in batch]
        if call_log is not None:
            call_log.append(list(labels))
        # fake tokens: shape (batch, len+2) filled with zeros
        max_len = max(len(s) for s in seqs) + 2
        tokens = torch.zeros((len(batch), max_len), dtype=torch.int64)
        return labels, seqs, tokens

    alphabet_mock = MagicMock()
    alphabet_mock.get_batch_converter.return_value = batch_converter

    class FakeModel:
        num_layers = 1

        def eval(self):
            return self

        def to(self, device):
            return self

        def __call__(self, tokens, repr_layers=None):
            import torch

            n, seq_len = tokens.shape
            arr = torch.ones((n, seq_len, dim), dtype=torch.float32) * 0.5
            return {"representations": {1: arr}}

    return FakeModel(), alphabet_mock


def _write_checkpoint_shard(checkpoint_dir, shard_name, ids, vectors, model_name):
    """Pre-populate a checkpoint dir as if a prior run wrote it, for resume tests."""
    checkpoint_dir = Path(checkpoint_dir)
    checkpoint_dir.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        checkpoint_dir / shard_name,
        ids=np.array(ids, dtype=object),
        vectors=np.stack(vectors).astype(np.float32),
    )
    (checkpoint_dir / "manifest.json").write_text(
        json.dumps({"model_name": model_name})
    )


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
    fasta = tmp_path / "test.faa"
    fasta.write_text(">seq1\nMFIXEDSEQ\n>seq2\nMANOTHERONE\n")
    out = tmp_path / "out.npz"

    from workflow.bgcflow.bgcflow.features.embed_esm2 import embed_fasta_to_npz

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
    fasta = tmp_path / "one.faa"
    fasta.write_text(">single\nMONLYONE\n")
    out = tmp_path / "one.npz"

    from workflow.bgcflow.bgcflow.features.embed_esm2 import embed_fasta_to_npz

    embed_fasta_to_npz(str(fasta), str(out), batch_size=100)

    data = np.load(str(out), allow_pickle=True)
    assert len(data["ids"]) == 1


def test_resume_skips_already_checkpointed_sequences(tmp_path):
    """A sequence already present in the checkpoint must not be re-embedded."""
    from workflow.bgcflow.bgcflow.features.embed_esm2 import embed_sequences

    checkpoint_dir = tmp_path / "ckpt"
    prechecked_vector = np.full(DIM, 0.25, dtype=np.float32)
    prechecked_vector /= np.linalg.norm(prechecked_vector)
    _write_checkpoint_shard(
        checkpoint_dir,
        "shard_000000.npz",
        ids=["seq1"],
        vectors=[prechecked_vector],
        model_name="stub",
    )

    call_log = []
    model, alphabet = _make_stub_model(call_log=call_log)
    with patch(
        "workflow.bgcflow.bgcflow.features.embed_esm2.load_esm_model",
        return_value=(model, alphabet),
    ):
        ids, vectors = embed_sequences(
            [("seq1", "MFIXEDSEQ"), ("seq2", "MANOTHERONE")],
            model_name="stub",
            batch_size=8,
            checkpoint_dir=str(checkpoint_dir),
        )

    all_embedded_labels = [label for batch in call_log for label in batch]
    assert "seq1" not in all_embedded_labels, "already-checkpointed id was re-embedded"
    assert "seq2" in all_embedded_labels

    assert ids == ["seq1", "seq2"]
    np.testing.assert_allclose(vectors[0], prechecked_vector, atol=1e-6)


def test_resume_preserves_original_input_order(tmp_path):
    """Output order must match input order, even when the middle id was checkpointed."""
    from workflow.bgcflow.bgcflow.features.embed_esm2 import embed_sequences

    checkpoint_dir = tmp_path / "ckpt"
    prechecked_vector = np.full(DIM, 0.25, dtype=np.float32)
    prechecked_vector /= np.linalg.norm(prechecked_vector)
    _write_checkpoint_shard(
        checkpoint_dir,
        "shard_000000.npz",
        ids=["seq2"],
        vectors=[prechecked_vector],
        model_name="stub",
    )

    model, alphabet = _make_stub_model()
    with patch(
        "workflow.bgcflow.bgcflow.features.embed_esm2.load_esm_model",
        return_value=(model, alphabet),
    ):
        ids, vectors = embed_sequences(
            [("seq1", "AAA"), ("seq2", "BBB"), ("seq3", "CCC")],
            model_name="stub",
            batch_size=8,
            checkpoint_dir=str(checkpoint_dir),
        )

    assert ids == ["seq1", "seq2", "seq3"]
    np.testing.assert_allclose(vectors[1], prechecked_vector, atol=1e-6)


def test_checkpoint_shards_written_during_run(tmp_path):
    """Checkpoint shards must appear on disk before the run finishes (crash-safety)."""
    from workflow.bgcflow.bgcflow.features.embed_esm2 import embed_sequences

    checkpoint_dir = tmp_path / "ckpt"
    model, alphabet = _make_stub_model()
    with patch(
        "workflow.bgcflow.bgcflow.features.embed_esm2.load_esm_model",
        return_value=(model, alphabet),
    ):
        embed_sequences(
            [("s1", "AAA"), ("s2", "BBB"), ("s3", "CCC")],
            model_name="stub",
            batch_size=1,
            checkpoint_dir=str(checkpoint_dir),
            checkpoint_every=1,
        )

    shard_files = sorted(checkpoint_dir.glob("shard_*.npz"))
    assert len(shard_files) >= 3, "expected one shard per batch at checkpoint_every=1"
    all_ids = set()
    for shard in shard_files:
        data = np.load(shard, allow_pickle=True)
        all_ids.update(data["ids"].tolist())
    assert all_ids == {"s1", "s2", "s3"}


def test_checkpoint_model_mismatch_discards_stale_checkpoint(tmp_path):
    """A checkpoint built with a different model must be ignored, not merged in."""
    from workflow.bgcflow.bgcflow.features.embed_esm2 import embed_sequences

    checkpoint_dir = tmp_path / "ckpt"
    stale_vector = np.full(DIM, 0.9, dtype=np.float32)
    stale_vector /= np.linalg.norm(stale_vector)
    _write_checkpoint_shard(
        checkpoint_dir,
        "shard_000000.npz",
        ids=["seq1"],
        vectors=[stale_vector],
        model_name="some-other-model",
    )

    call_log = []
    model, alphabet = _make_stub_model(call_log=call_log)
    with patch(
        "workflow.bgcflow.bgcflow.features.embed_esm2.load_esm_model",
        return_value=(model, alphabet),
    ):
        ids, vectors = embed_sequences(
            [("seq1", "MFIXEDSEQ")],
            model_name="stub",
            batch_size=8,
            checkpoint_dir=str(checkpoint_dir),
        )

    all_embedded_labels = [label for batch in call_log for label in batch]
    assert (
        "seq1" in all_embedded_labels
    ), "stale checkpoint from a different model was reused"
    manifest = json.loads((checkpoint_dir / "manifest.json").read_text())
    assert manifest["model_name"] == "stub"


def test_long_sequence_isolated_from_short_batch():
    """A very long sequence must not drag short sequences into a huge padded batch.

    Self-attention memory scales with the square of the padded batch length,
    so grouping one 2000-residue outlier with 7 short sequences would force
    all 8 through a 2000-token-wide forward pass.
    """
    from workflow.bgcflow.bgcflow.features.embed_esm2 import embed_sequences

    long_seq = ("long1", "M" * 2000)
    short_seqs = [(f"short{i}", "M" * 10) for i in range(7)]
    sequences = short_seqs + [long_seq]

    call_log = []
    model, alphabet = _make_stub_model(call_log=call_log)
    with patch(
        "workflow.bgcflow.bgcflow.features.embed_esm2.load_esm_model",
        return_value=(model, alphabet),
    ):
        ids, vectors = embed_sequences(
            sequences, model_name="stub", batch_size=8, max_tokens_per_batch=500
        )

    long_batch = next(b for b in call_log if "long1" in b)
    assert long_batch == [
        "long1"
    ], f"long sequence was batched with others: {long_batch}"
    assert len(ids) == 8


def test_short_sequences_still_batch_up_to_batch_size():
    """Without long outliers, batching still respects batch_size as before."""
    from workflow.bgcflow.bgcflow.features.embed_esm2 import embed_sequences

    sequences = [(f"s{i}", "M" * 20) for i in range(17)]

    call_log = []
    model, alphabet = _make_stub_model(call_log=call_log)
    with patch(
        "workflow.bgcflow.bgcflow.features.embed_esm2.load_esm_model",
        return_value=(model, alphabet),
    ):
        embed_sequences(
            sequences, model_name="stub", batch_size=8, max_tokens_per_batch=100_000
        )

    batch_sizes = sorted((len(b) for b in call_log), reverse=True)
    assert batch_sizes == [8, 8, 1], f"unexpected batch grouping: {batch_sizes}"


def test_oversized_single_sequence_processed_alone():
    """A sequence longer than the token budget by itself must still be processed, not skipped."""
    from workflow.bgcflow.bgcflow.features.embed_esm2 import embed_sequences

    huge_seq = ("huge1", "M" * 5000)
    model, alphabet = _make_stub_model()
    with patch(
        "workflow.bgcflow.bgcflow.features.embed_esm2.load_esm_model",
        return_value=(model, alphabet),
    ):
        ids, vectors = embed_sequences(
            [huge_seq], model_name="stub", batch_size=8, max_tokens_per_batch=100
        )

    assert ids == ["huge1"]
    assert vectors.shape[0] == 1
