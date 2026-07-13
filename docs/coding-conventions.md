# Coding Conventions

Scope: repo-specific conventions for `bgcflow-plm`, not already enforced by tooling.
Summary: rules covering `.smk` rule files, config toggles, log paths, and additive side-channels — things a linter can't check for you.

---

Style mechanics (line length, import order, trailing whitespace, tabs/CRLF) are enforced by `.pre-commit-config.yaml` (black, isort, flake8 with `E501,W503` ignored). Don't restate those here; this doc only covers what tooling doesn't catch.

## `.smk` rule files

- One `.smk` file per tool, self-contained, under `workflow/rules/`. Don't split a single tool's rules across files or merge unrelated tools into one file.
- A rule file that adds a **side-channel** (reads existing outputs, doesn't sit on the main DAG path — e.g. `plm_novelty.smk`) must never write into another tool's output directory, and should say so in its module docstring, the way `plm_novelty.smk` does.
- Guard version-dependent rule registration at the top of the file (see the `antismash_major_version` check in `plm_novelty.smk`) rather than inside individual rules — keeps the "why no rule exists" reasoning in one place.
- Large external databases are installed as rule outputs under `resources/`, referenced via symlinks when the user supplies an external path. Follow `custom_resource_dir()` in `common.smk` for this pattern rather than hand-rolling a new one.

## Config toggles

- New optional features get a boolean under `rules:` in `config/config.yaml`, and any tunable parameters live under `rule_parameters.<tool>` (see `rule_parameters.plm` for the pattern: model name, batch size, device, thresholds — not paths, which belong in `resources_path` or rule outputs).
- Validate new config sections against `workflow/schemas/config.schema.yaml` — don't add a config key without a corresponding schema entry.

## Logs

- Every rule writes its log to `logs/<tool>/<rule>-<wildcards>.log`. Don't write logs elsewhere or omit them — downstream debugging assumes this path.

## Wildcards

- Reuse the existing wildcard vocabulary (`{name}`, `{strains}`, `{strains_fna}`, `{strains_genbank}`, `{ncbi}`, `{patric}`, `{custom_fna}`, `{custom_genbank}`, `{version}` — defined in `CLAUDE.md`) rather than inventing new wildcard names for the same concepts.

## Tests

See `docs/testing.md` for suite locations and conventions. In short: a new rule-script gets a unit test in `.tests/unit/`, a new `bgcflow_wrapper` feature gets a test in `workflow/bgcflow/tests/`, and a new PLM-novelty component gets a test in `tests/plm_novelty/`.
