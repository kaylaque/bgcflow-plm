"""Tests for plm_novelty_query.py — uses fixture-built FAISS index (no ESM2)."""
import math
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

DIM = 8
N_MIBIG = 4  # matches mini_mibig fixture (2 GBKs × 2 proteins)


@pytest.fixture
def fixture_index(tmp_path):
    """Build a tiny FAISS IndexFlatIP from known vectors, save index + metadata."""
    import faiss

    rng = np.random.default_rng(0)
    vecs = rng.random((N_MIBIG, DIM)).astype(np.float32)
    vecs /= np.linalg.norm(vecs, axis=1, keepdims=True)

    index = faiss.IndexFlatIP(DIM)
    index.add(vecs)
    index_path = tmp_path / "index.faiss"
    faiss.write_index(index, str(index_path))

    meta = pd.DataFrame(
        {
            "mibig_id": [f"BGC000000{i}" for i in range(N_MIBIG)],
            "protein_id": [f"BGC000000{i}_001" for i in range(N_MIBIG)],
            "product_class": ["terpene"] * N_MIBIG,
            "sequence_hash": [f"hash{i}" for i in range(N_MIBIG)],
        }
    )
    meta_path = tmp_path / "metadata.parquet"
    meta.to_parquet(meta_path, index=False)

    return index_path, meta_path, vecs  # expose raw vecs for query crafting


def _save_npz(path, ids, vectors):
    np.savez_compressed(path, ids=np.array(ids, dtype=object), vectors=vectors.astype(np.float32))


def test_identical_vector_not_novel(fixture_index, tmp_path):
    """Query vector identical to an indexed one → nearest_distance≈0, novelty_flag=False."""
    from workflow.bgcflow.bgcflow.features.plm_novelty_query import query_index

    index_path, meta_path, vecs = fixture_index
    query_vec = vecs[[0]]  # identical to first MIBiG vector

    npz = tmp_path / "q.npz"
    _save_npz(npz, ["strain1|region001|prot001"], query_vec)

    per_protein, _ = query_index(
        [str(npz)],
        str(index_path),
        str(meta_path),
        top_k=1,
        distance_threshold=0.5,
    )
    assert len(per_protein) == 1
    assert per_protein.iloc[0]["nearest_distance"] < 0.01
    assert per_protein.iloc[0]["novelty_flag"] == False


def test_orthogonal_vector_novel(fixture_index, tmp_path):
    """Orthogonal query vector → nearest_distance≈1, novelty_flag=True at threshold 0.5."""
    import faiss

    from workflow.bgcflow.bgcflow.features.plm_novelty_query import query_index

    index_path, meta_path, vecs = fixture_index

    # Build an orthogonal vector by using the null space complement
    rng = np.random.default_rng(99)
    q = rng.random((1, DIM)).astype(np.float32)
    # Zero out projection onto all MIBiG vectors
    for v in vecs:
        q -= np.dot(q, v) * v
    norm = np.linalg.norm(q)
    if norm < 1e-6:
        # Fallback: use a vector very far from all MIBiG vecs
        q = np.zeros((1, DIM), dtype=np.float32)
        q[0, 0] = 1.0
    else:
        q /= norm

    npz = tmp_path / "orth.npz"
    _save_npz(npz, ["strain1|region002|prot002"], q)

    per_protein, _ = query_index(
        [str(npz)],
        str(index_path),
        str(meta_path),
        top_k=1,
        distance_threshold=0.5,
    )
    assert per_protein.iloc[0]["nearest_distance"] > 0.5
    assert per_protein.iloc[0]["novelty_flag"] == True


def test_region_summary_math(fixture_index, tmp_path):
    """4 core proteins in a region, 2 with close hits → frac=0.5, score=0.5."""
    from workflow.bgcflow.bgcflow.features.plm_novelty_query import query_index

    index_path, meta_path, vecs = fixture_index

    # 2 identical (close hit) + 2 orthogonal (novel)
    close = vecs[:2]  # guaranteed close hits
    far = np.zeros((2, DIM), dtype=np.float32)
    far[0, 0] = 1.0
    far[1, 1] = 1.0
    for v in vecs:
        far -= np.outer(np.einsum("ij,j->i", far, v), v)
    for i in range(len(far)):
        n = np.linalg.norm(far[i])
        if n > 1e-6:
            far[i] /= n
    all_vecs = np.concatenate([close, far], axis=0)

    ids = [
        "strainX|regionA|p1",
        "strainX|regionA|p2",
        "strainX|regionA|p3",
        "strainX|regionA|p4",
    ]
    npz = tmp_path / "region.npz"
    _save_npz(npz, ids, all_vecs)

    per_protein, per_region = query_index(
        [str(npz)],
        str(index_path),
        str(meta_path),
        top_k=1,
        distance_threshold=0.5,
    )

    assert len(per_region) == 1
    row = per_region.iloc[0]
    assert row["n_core_proteins"] == 4
    assert 0.0 <= row["frac_with_close_hit"] <= 1.0
    assert abs(row["region_novelty_score"] - (1.0 - row["frac_with_close_hit"])) < 1e-6


def test_empty_npz_does_not_crash(fixture_index, tmp_path):
    """Empty embeddings file → zero rows in output, no exception."""
    from workflow.bgcflow.bgcflow.features.plm_novelty_query import query_index

    index_path, meta_path, _ = fixture_index
    empty_npz = tmp_path / "empty.npz"
    np.savez_compressed(
        str(empty_npz),
        ids=np.array([], dtype=object),
        vectors=np.zeros((0, DIM), dtype=np.float32),
    )
    per_protein, per_region = query_index(
        [str(empty_npz)],
        str(index_path),
        str(meta_path),
        top_k=1,
        distance_threshold=0.5,
    )
    assert len(per_protein) == 0
    assert len(per_region) == 0


def test_multi_strain_keyed_by_strain_id(fixture_index, tmp_path):
    """Two .npz files from different strains → both strain_ids in per_protein."""
    from workflow.bgcflow.bgcflow.features.plm_novelty_query import query_index

    index_path, meta_path, vecs = fixture_index

    npz1 = tmp_path / "s1.npz"
    npz2 = tmp_path / "s2.npz"
    _save_npz(npz1, ["strain_A|reg001|p1"], vecs[[0]])
    _save_npz(npz2, ["strain_B|reg001|p1"], vecs[[1]])

    per_protein, _ = query_index(
        [str(npz1), str(npz2)],
        str(index_path),
        str(meta_path),
        top_k=1,
        distance_threshold=0.5,
    )
    strain_ids = set(per_protein["strain_id"])
    assert "strain_A" in strain_ids
    assert "strain_B" in strain_ids


def test_output_schema_per_protein(fixture_index, tmp_path):
    """Per-protein TSV must have the exact design §7 column set."""
    from workflow.bgcflow.bgcflow.features.plm_novelty_query import query_index

    index_path, meta_path, vecs = fixture_index
    npz = tmp_path / "schema.npz"
    _save_npz(npz, ["strainZ|region001|prot1"], vecs[[0]])

    out_pp = tmp_path / "per_protein.tsv"
    out_pr = tmp_path / "per_region.tsv"
    query_index(
        [str(npz)],
        str(index_path),
        str(meta_path),
        top_k=1,
        distance_threshold=0.5,
        out_per_protein=str(out_pp),
        out_per_region=str(out_pr),
    )

    pp = pd.read_csv(out_pp, sep="\t")
    expected_pp = {
        "strain_id",
        "region_id",
        "protein_id",
        "top_k_mibig_ids",
        "cosine_distances",
        "nearest_distance",
        "confidence_band",
        "novelty_flag",
    }
    assert set(pp.columns) == expected_pp, f"Unexpected columns: {set(pp.columns)}"

    pr = pd.read_csv(out_pr, sep="\t")
    expected_pr = {
        "strain_id",
        "region_id",
        "n_core_proteins",
        "frac_with_close_hit",
        "region_novelty_score",
        "flagged_novel",
    }
    assert set(pr.columns) == expected_pr, f"Unexpected columns: {set(pr.columns)}"
