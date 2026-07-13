# Testing

Scope: index of every test suite in this repo — what it covers, how to run it, where a new test belongs.
Summary: three independent pytest suites (`workflow/bgcflow/tests`, `tests/plm_novelty`, `.tests/unit`), each covering a different layer of the pipeline.

---

## `workflow/bgcflow/tests/`

Unit tests for the `bgcflow_wrapper` pip package (`workflow/bgcflow/bgcflow/`) — database aggregation, feature extraction, visualization helpers, and the PLM novelty query module (`test_plm_novelty.py`).

```bash
cd workflow/bgcflow
python -m pytest tests/
```

New test for a change under `workflow/bgcflow/bgcflow/` → goes here, named `test_<module>.py`.

## `tests/plm_novelty/`

Suite for the PLM novelty-retrieval side-channel specifically: index building (`test_build_mibig_esm2_index.py`), embedding (`test_embed_esm2.py`), core-protein extraction (`test_extract_bgc_core_proteins.py`), the query path (`test_plm_novelty_query.py`), config schema validation (`test_config_schema.py`), and a guard test (`test_no_upstream_dag_edit.py`) asserting the side-channel never touches antiSMASH's own output directories — keep that guard test passing for any change to `plm_novelty.smk`.

```bash
pytest tests/plm_novelty
```

Fixtures live in `tests/plm_novelty/fixtures/`. New test for a PLM-novelty change → goes here.

## `.tests/unit/`

Unit tests for individual Snakemake rule scripts (ARTS combine/extract steps, FastANI conversion, GTDB taxonomy fixing, antiSMASH summary/overview, deepTFactor JSON conversion, NCBI info extraction). Each script has a matching `test_<script>.py` plus a same-named fixture directory holding expected inputs/outputs.

```bash
pytest .tests/unit
```

New test for a new or changed rule script under `workflow/scripts/` → goes here, following the existing `test_<script>.py` + `<script>/` fixture-dir pattern.

## `.tests/` (integration)

`.tests/config/` holds several full PEP project configs (`test1`, `test2`, `test3`, `test_pep`, `lanthipeptide_lactobacillus`) used for end-to-end dry-run / smoke-run validation of the whole Snakemake DAG, not just unit-level scripts. Use these when a change could affect DAG construction or cross-rule wiring, not just a single script's logic.

## Running everything

There's no single top-level `pytest` invocation that covers all three suites at once (different working directories, different conda environments per the package boundary) — run each suite from its own location as shown above.
