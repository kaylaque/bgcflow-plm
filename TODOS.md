# TODOS

Flat, current task queue. Not a durable record — once a task is done and has a worksheet (`docs/worksheets/`), drop the line here and let the worksheet be the record.

## In progress

- (none)

## Next

- Continue PLM novelty retrieval work per `docs/design-plm-retrieval-bgcflow.md` / `docs/systems/plm-novelty.md`
- Consider a scoped follow-up spec for the deferred items in `docs/superpowers/specs/2026-07-13-agentic-infra-design.md` (cross-agent review tooling, pre-commit auto-fix, `tools/bin/` scripts) if they prove necessary
- Fix `workflow/bgcflow/tests/test_plm_novelty.py`'s off-by-one `FEATURES_DIR` path (points at `workflow/bgcflow/features/`, should be `workflow/bgcflow/bgcflow/features/`) — 12 tests currently fail to even collect
- Get `test_config_schema.py` runnable in the `plm_novelty` conda env (needs `jsonschema`), or move it to run under the `bgcflow` env instead
- Actually run the (now resumable) MIBiG index build in tmux: `snakemake --use-conda --cores 8 resources/mibig_esm2/index.faiss`

## Done (recent)

- Scoped and built docs/workflow agentic infra tier (AGENTS.md, docs/systems/, docs/coding-conventions.md, docs/testing.md, docs/worksheets/, docs/agent-feedback.md, TODOS.md) — see `docs/worksheets/2026-07-13-agentic-infra-setup.md`
- Added checkpoint/resume to MIBiG ESM2 embedding + fixed 3 pre-existing broken PLM test mocks found along the way — see `docs/worksheets/2026-07-13-mibig-embed-checkpoint.md`
