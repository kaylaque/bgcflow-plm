# PLM Retrieval Side-Channel for BGCFlow

**Adding a protein-language-model augmentation layer at function annotation and comparative genomics**

Status: Design + 1-day build plan · Scope: prototype plumbing, not validated science

---

## 1. Motivation

BGCFlow's recall is capped by its reference databases. Every stage that matters for novelty — antiSMASH HMMs, MIBiG comparison, Pfam, GTDB — recognizes only what resembles the already-known. Out-of-distribution BGC variants ("variasi diluar nalar") are structurally invisible, because homology-matching *is* the ceiling.

Protein language models (PLMs) widen the reachable radius: ESM-style embeddings encode remote homology, so a protein doing the same chemistry at <20% sequence identity can still be retrieved. The goal here is **not** to replace any deterministic call. It is to bolt on a *novelty-tolerant retrieval side-channel* that surfaces what the pipeline would otherwise miss, while preserving the audit trail bioinformatics reproducibility requires.

The framing analogy is RAG done correctly: the database stays the grounded source of truth; the PLM is the *retriever* that searches it by meaning instead of by exact sequence. Retrieval raises transparency (you keep a named reference to point at); it must never lower it by making an opaque call inside the critical path.

---

## 2. Design principle — the validity firewall

> The PLM is an **evidence/ranking layer that sits alongside the deterministic DAG, never inside it.**

- It **consumes** existing antiSMASH output; it does not modify the DAG upstream of itself.
- It **never overwrites** an antiSMASH call. Its output is a *separate novelty axis*.
- Every output is auditable: nearest known neighbor(s), embedding distance, confidence band, provenance.
- Determinism is preserved by pinning the model weights + version as a pipeline dependency (inference is deterministic given fixed weights).

This is the "wrap, don't replace / learned novelty layer" strategy realized as a reproducible workflow rule — which is precisely the part no existing tool has shipped.

---

## 3. Where PLM work maps onto BGCFlow stages

| BGCFlow stage | Existing PLM research | Replace / augment | 1-day feasibility |
|---|---|---|---|
| Data selection (QC, dereplication) | none meaningful (CheckM2 / GTDB-Tk / dRep) | — | skip |
| **Function annotation** | PLMSearch, TM-Vec/DeepBLAST, npj 2026 embedding retrieval, PLM-eXplain | **augment** | **strong — inference-only** ← primary |
| Phylogenetic analysis | deterministic (autoMLST / GTDB-Tk); PLM phylogeny niche & contested | — | skip |
| Genome mining (detection + class) | CoreFinder, BGC-Finder, BGC-Prophet, BiGCARP, DeepBGC, GECCO; BGC-MAC/MAP, BGCat | mostly **replace** | bad — needs training + GPU + benchmark |
| **Comparative genomics** (family clustering) | BiG-SCAPE/BiG-SLICE deterministic; embedding cosine distance groups divergent families | **augment** | feasible as a freebie ← secondary |

**Read of the map:** genome mining is where everyone piles in and the one stage that cannot be done credibly in a day. Function annotation is the highest-value, lowest-risk insertion point, and it is the least crowded when framed as *retrieval* rather than *replacement*. Comparative genomics comes almost for free from the same embeddings.

---

## 4. Insertion point A — Function annotation (primary)

**What it does:** for each antiSMASH-detected region, retrieve the closest known biosynthetic proteins in MIBiG by embedding proximity, and flag proteins/regions whose nearest known neighbor is far away (candidate novelty).

**Pipeline:**

1. Read existing antiSMASH region GBKs (BGCFlow already produces these). No DAG changes.
2. Extract **biosynthetic-core** protein sequences per region; drop regulatory / transport / other genes (the BGCat trick — roughly halves compute and removes incidental noise).
3. Embed with **frozen ESM2** (inference only, no training).
   - GPU: `esm2_t33_650M_UR50D`
   - CPU fallback: `esm2_t30_150M` or `esm2_t12_35M`
4. Query each embedding (cosine / FAISS) against a **precomputed MIBiG reference embedding index** (built once, cached).
5. Emit per-protein top-k neighbors + distance + confidence band + novelty flag.

**False-positive control (the caveat that killed ClusterFinder):** individual lucky hits are unreliable and coincidental co-occurrence is real. Mitigate exactly as the npj 2026 paper did — **aggregate to region/strain level and rank by the *fraction* of core proteins with close hits**, not by single protein hits.

---

## 5. Insertion point B — Comparative genomics (freebie)

Once every core protein is embedded, you already hold everything needed to group **divergent BGC families** that BiG-SCAPE would split apart:

- Mean-pool core-protein embeddings → one region-level vector.
- Cluster regions by embedding cosine distance → a novelty-tolerant complement to BiG-SCAPE GCFs.

**Do not build this on day 1.** Design the rule so this is a one-function add on top of the same embedding cache (near-zero extra compute). It is the second output of the same machinery.

---

## 6. One-day build plan

Scope for the day: a **working, wired-in retrieval rule** that produces the side-channel TSV. Not a validated detector.

**Sequence of work:**

1. **Build MIBiG reference index (do first, cache forever).**
   Download MIBiG (~2,500 BGCs) → extract core proteins → embed once → build FAISS index + neighbor metadata table. ~1–2 hrs, then never repeated. This is the only real bottleneck.
2. **Query-side extraction.** GBK → biosynthetic-core FASTA per region (filter gene kinds).
3. **Embedding step.** Frozen ESM2 batch inference over query proteins → vectors.
4. **Retrieval + scoring.** FAISS top-k vs MIBiG → distances → per-region novelty aggregation (fraction rule).
5. **Output writer.** Emit `bgc_novelty_retrieval.tsv` + a small per-region summary.
6. **Snakemake integration.** One additive rule in `workflow/rules/`, config toggle, model version pinned. Consumes antiSMASH output, writes to a new results subdir.

**Rule sketch (spec for Claude Code to scaffold):**

```
rule plm_novelty_retrieval:
    input:
        gbk   = "data/interim/antismash/{strain}/",       # existing BGCFlow output
        index = "resources/mibig_esm2_faiss.index",       # prebuilt, cached
        meta  = "resources/mibig_neighbor_meta.parquet"
    output:
        tsv     = "data/processed/plm_novelty/{strain}/bgc_novelty_retrieval.tsv",
        summary = "data/processed/plm_novelty/{strain}/region_novelty_summary.tsv"
    params:
        model      = config["plm"]["model"],              # pinned, e.g. esm2_t33_650M
        top_k      = config["plm"]["top_k"],
        core_only  = True,                                # BGCat-style gene filter
        novelty_thr= config["plm"]["distance_threshold"]
    script:
        "../scripts/plm_novelty_retrieval.py"
```

**Where the 20x Max Claude sub pays off:** the day's real cost is plumbing, not ideas — GBK→FASTA extraction, FAISS wiring, batching, the output schema, Snakemake config integration. That is exactly the high-volume boilerplate to hand to Claude Code first thing. Point it at BGCFlow's existing rule structure and let it scaffold steps 2–6 against this spec.

---

## 7. Output schema

`bgc_novelty_retrieval.tsv` (per query core protein):

| column | meaning |
|---|---|
| strain_id | source genome |
| region_id | antiSMASH region |
| protein_id | query core protein |
| top_k_mibig_ids | nearest known neighbors (provenance) |
| cosine_distances | distance to each neighbor |
| nearest_distance | min distance (novelty signal) |
| confidence_band | high / medium / low (binned distance) |
| novelty_flag | nearest_distance > threshold |

`region_novelty_summary.tsv` (per region):

| column | meaning |
|---|---|
| region_id | antiSMASH region |
| n_core_proteins | core genes considered |
| frac_with_close_hit | fraction with a close MIBiG neighbor (FP control) |
| region_novelty_score | 1 − frac_with_close_hit |
| flagged_novel | boolean, thresholded |

Everything is inspectable; nothing overwrites an antiSMASH field.

---

## 8. Scope boundaries & non-goals (anti-overclaim)

- **Day 1 delivers plumbing, not validated science.** Whether the novelty flag tracks *real* novelty is a separate evaluation against held-out known BGCs, run afterward. Ship the wiring; do not claim the biology yet.
- **The database dependency does not disappear.** ESM2 was itself trained on UniProt-scale data, so truly alien variants that resemble nothing in training remain hard. The PLM widens the radius; it does not make it infinite.
- **No replacement of any antiSMASH call.** The layer is additive only.
- **No phylogenetics, no de novo detection, no product-structure prediction** in this scope.

---

## 9. Related work grounding

- **Replacement direction (crowded — avoid competing here):** CoreFinder / BGC-Finder (ESM + genomic context, alignment-free detection, ~100–240× faster than antiSMASH, remote-homolog retrieval, one wet-lab-validated hit); BGC-Prophet (transformer + ESM embeddings); BiGCARP (ESM-1b Pfam-domain masked LM); DeepBGC; GECCO.
- **Augmentation direction (matches this plan):** npj 2026 — PLM embedding retrieval against MIBiG rescued a neomycin cluster antiSMASH missed in a fragmented genome, and controlled false positives via strain-level fraction ranking. CHAMOIS (chemical-fingerprint screen over antiSMASH/GECCO output). BGCat (ESM embedding of antiSMASH regions, core-gene-only).
- **Interpretability thread (validates the transparency crux):** BGC-MAC/BGC-MAP (cross-attention as explainable AI, surfaces the domains driving a call); PLM-eXplain (partitions the ESM embedding space so feature attribution becomes biologically meaningful).

**White space this plan occupies:** the above are standalone models or one-off analyses. None ships as a reproducible, provenance-carrying component integrated into an orchestration people actually run. A BGCFlow rule that emits an auditable novelty side-channel is the underoccupied contribution.

---

## 10. After day 1

1. Add the comparative-genomics clustering output (Insertion point B) from the cached embeddings.
2. Run the honest evaluation: does `region_novelty_score` separate held-out known BGCs from genuinely divergent ones?
3. Calibrate the distance threshold and confidence bands against that evaluation before trusting the flag.
4. Only then surface the novelty axis in any downstream interpretation layer.

---

## 11. Implementation notes (as-built)

**antiSMASH version requirement:** the lane is v7-only. `plm_novelty.smk` checks `antismash_major_version` at include time and registers no rules (printing a warning) on v6. This avoids a hard failure when users run older antiSMASH environments.

**MIBiG source:** rather than a separate MIBiG download, `build_mibig_esm2_index` reuses the reference GBK dump that `antismash_db_setup` already downloads to `resources/antismash_db/knownclusterblast/4.0/`. This avoids duplicating ~1 GB and keeps the build rule trivially connected to the existing dependency graph.

**Default model:** `esm2_t30_150M_UR50D` (dim 480, CPU-tractable, ~600 MB). To use the GPU model from the design spec, set `rule_parameters.plm.model: esm2_t33_650M_UR50D` in `config.yaml`. The FAISS index must then be rebuilt (delete `resources/mibig_esm2/` and re-run).

**Enabling the lane:** add `plm-novelty: TRUE` to a project's `rules:` block in its `project_config.yaml` (PEP format). The global `rules.yaml` entry wires the final output automatically through `get_project_outputs` in `common.smk` — no other changes needed.

**Model consistency:** `plm_novelty_query` reads `resources/mibig_esm2/model_version.txt` at runtime and logs a warning if the model used to build the index differs from the one specified in `config.yaml`. This prevents silent embedding-space mismatch if the config is changed after the index is cached.

**Dim mismatch guard:** if query embedding vectors have a different dimension than the FAISS index (e.g. because index and embeddings were built with different models), `plm_novelty_query` raises `ValueError` immediately rather than producing silently wrong cosine distances.
