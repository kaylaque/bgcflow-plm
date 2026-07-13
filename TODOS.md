# TODOS

Flat, current task queue. Not a durable record — once a task is done and has a worksheet (`docs/worksheets/`), drop the line here and let the worksheet be the record.

## In progress

- (none)

## Next

- Continue PLM novelty retrieval work per `docs/design-plm-retrieval-bgcflow.md` / `docs/systems/plm-novelty.md`
- Consider a scoped follow-up spec for the deferred items in `docs/superpowers/specs/2026-07-13-agentic-infra-design.md` (cross-agent review tooling, pre-commit auto-fix, `tools/bin/` scripts) if they prove necessary
- Fix `workflow/bgcflow/tests/test_plm_novelty.py`'s off-by-one `FEATURES_DIR` path (points at `workflow/bgcflow/features/`, should be `workflow/bgcflow/bgcflow/features/`) — 12 tests currently fail to even collect
- Fix `test_config_schema.py`: even in the `bgcflow` env (which has `jsonschema`), it fails with a `SchemaError` — `workflow/schemas/config.schema.yaml` has an invalid `'type': ['string', 'None']` on `projects.items.prokka-db` (`'None'` isn't a valid JSON Schema type, should be `'null'`). Pre-existing, confirmed via `git stash` against HEAD before today's changes.
- Snakemake's own `validate()` only warns (doesn't fail) on the schema's `$schema` draft version mismatch, so this hasn't blocked real runs — but it should still be fixed properly.

## Done (recent)

- Scoped and built docs/workflow agentic infra tier (AGENTS.md, docs/systems/, docs/coding-conventions.md, docs/testing.md, docs/worksheets/, docs/agent-feedback.md, TODOS.md) — see `docs/worksheets/2026-07-13-agentic-infra-setup.md`
- Added checkpoint/resume to MIBiG ESM2 embedding + fixed 3 pre-existing broken PLM test mocks found along the way — see `docs/worksheets/2026-07-13-mibig-embed-checkpoint.md`
