# Agentic Workflow Infrastructure — Docs & Workflow Tier

**Scoped subset of an 18-item agentic-setup checklist, tailored to this repo's actual shape**

Status: Approved · Scope: docs/workflow infra only (no linter automation, no review/loop skills, no perf/visual-regression tooling)

---

## 1. Motivation

The user proposed a comprehensive 18-item agentic-coding setup (router doc, self-healing docs, always-run-the-app, e2e tests, auto-fixing pre-commit, cross-agent review, agent worksheets/traces, session feedback loop, tools/bin scripts, commit sweeps, coding conventions, night-shift loop, task queue, false-confidence test audits, visual regression, perf benchmarks/profiling, end-of-shift validation) and asked for an evaluation of fit followed by implementation of what's necessary.

This repo (`bgcflow-plm`) is a solo/small-team Snakemake bioinformatics pipeline with no UI. The user confirmed they also drive Codex and Cursor against this repo (not Claude Code alone), so docs need to be legible to any of those tools, not just Claude-specific conventions. Several items on the original list assume a team shipping a product with a browser-facing UI and heavy perf sensitivity — that doesn't describe this repo, so building visual regression tests or perf-benchmark infrastructure here would be dead weight.

The user selected the **"Docs & workflow infra"** tier: router doc, self-healing per-system docs with greppable headers, coding conventions doc, test index doc, a flat task queue, agent worksheets tied to git tags, and an end-of-session feedback log. Deferred: pre-commit auto-fix automation, a `tools/bin/` scripts folder, cross-agent review scripts/skill, night-shift autonomous loop, commit sweeps, false-confidence test audits, visual regression, perf benchmarking/profiling, end-of-shift full validation. These can be layered on later as separate scoped specs if they prove necessary.

---

## 2. Fit assessment (why each of the 18 items was in/out of this pass)

| # | Item | Fit here | Disposition |
|---|------|----------|-------------|
| 0 | AGENTS.md router | High — no router exists today | **Build now** |
| 1 | Standard workflow doc/skill | Low marginal value — `superpowers` skills (brainstorming, TDD, systematic-debugging, writing-plans) already cover process; duplicating into `AGENT_WORKFLOW.md` would drift out of sync | Folded into AGENTS.md as pointers, not duplicated |
| 2 | Self-healing docs, greppable headers | High — `docs/design-plm-retrieval-bgcflow.md` already shows the pattern works here | **Build now** (`docs/systems/`) |
| 3 | Agents always run the app | Partial fit — this is a Snakemake pipeline wrapping external bioinformatics tools (antiSMASH, Prokka, GTDB-Tk); "running the app" means `snakemake -n` / targeted rule runs / pytest, which is already covered by `run` and `verify` skills at the harness level | No new doc needed; noted in AGENTS.md as a pointer to existing skills |
| 4 | E2E tests + test docs | High — three test suites exist (`workflow/bgcflow/tests`, `tests/plm_novelty`, `.tests/unit`) with no index | **Build now** (`docs/testing.md`) |
| 5 | Auto-fixing pre-commit / LLM-fix hook | Deferred — real automation work, not a docs task; user explicitly deferred it | Not built this pass |
| 6 | Cross-agent review (Codex/Cursor/Claude) + personas | Deferred — user confirmed Codex/Cursor are used, so this is real, but it's a skill+script deliverable, not docs infra; user explicitly deferred it | Not built this pass |
| 7 | Agent worksheets/traces + git tags | High — directly requested, cheap, high handoff value across Claude/Codex/Cursor | **Build now** (`docs/worksheets/`) |
| 8 | End-of-session feedback log | High — cheap, user explicitly wants a periodic-ingest log | **Build now** (`docs/agent-feedback.md`) |
| 9 | `tools/bin/` helper scripts | Deferred — automation work, user deferred it | Not built this pass |
| 10 | Periodic commit sweeps | Deferred — this is a recurring *skill/loop* deliverable, not a doc; user deferred it | Not built this pass |
| 11 | Coding conventions doc | High — cheap, distinct from what black/isort/flake8 already enforce | **Build now** (`docs/coding-conventions.md`) |
| 12 | Night-shift/agent-loop skill | Deferred — user deferred it; no evidence of unattended runs today | Not built this pass |
| 13 | Task queue | High — user's own reference implementation is a flat `TODOS.md`; matches this repo's scale | **Build now** (`TODOS.md`) |
| 14 | False-confidence test audit | Deferred — a periodic skill, not docs infra; no evidence yet of misleading tests | Not built this pass |
| 15 | Visual regression tests | No fit — no UI in this repo (Snakemake pipeline + CLI wrapper) | Not built |
| 16 | Perf benchmark tests | Low fit — wall-clock is dominated by external bioinformatics tools (antiSMASH, GTDB-Tk), not code the agent tunes; revisit if the PLM embedding/query path becomes a hot path | Not built this pass |
| 17 | Perf profiling tools | Same as #16 | Not built this pass |
| 18 | End-of-shift full validation | Deferred — a skill/checklist for autonomous sessions; no evidence of unattended multi-hour sessions today | Not built this pass |

---

## 3. Components being built

### 3.1 `AGENTS.md` (repo root)
Router doc. Any agent (Claude, Codex, Cursor) reading the repo root should land here first (alongside `CLAUDE.md`, which stays as-is and is linked from AGENTS.md rather than duplicated). Contents:
- One-paragraph repo description + pointer to `CLAUDE.md` for architecture/workflow details
- Table of where to find things: coding conventions, testing, per-system docs, worksheets, feedback log, task queue
- The greppable-header convention: every doc under `docs/systems/`, plus `docs/coding-conventions.md` and `docs/testing.md`, opens with a short frontmatter-style block (topic, scope, one-line summary) in its first ~7 lines so `grep -A5 "^# "` across `docs/` surfaces the right file without opening it
- The worksheet workflow: when to create one, naming (`docs/worksheets/YYYY-MM-DD-<topic>.md`), the `worksheet/<topic>` git tag convention, and that tags are local-only unless the user asks to push
- Pointer to `docs/agent-feedback.md` and the append-only convention

### 3.2 `docs/systems/plm-novelty.md`
First system doc, seeded now because it's the live feature. Greppable header, then: what the PLM side-channel does, where its rule lives (`workflow/rules/plm_novelty.smk`), its scripts (`build_mibig_esm2_index.py`, `embed_esm2.py`, `plm_novelty_query.py`), config surface (`rule_parameters.plm` in `config/config.yaml`), outputs, and the v6/v7 antiSMASH gating behavior. This doc is explicitly documented (in AGENTS.md) as self-healing: the next agent that touches this system updates this file rather than letting it drift, but the doc doesn't try to duplicate `docs/design-plm-retrieval-bgcflow.md`'s rationale — it's the living reference, the design doc is the historical rationale.

New system docs are **not** pre-created for every subsystem — AGENTS.md documents that they're created the first time an agent does non-trivial work in an area lacking one, to avoid stale unreferenced docs.

### 3.3 `docs/coding-conventions.md`
Repo-specific conventions not already enforced by black/isort/flake8 (line length ignored via `E501`, tabs/CRLF forbidden by pre-commit). Covers things specific to this codebase: `.smk` rule file conventions (self-contained per tool, additive-only for side-channel rules like PLM), wildcard semantics reference pointer (already in CLAUDE.md, not duplicated), config toggle conventions (`rule_parameters.<tool>`), and log path conventions (`logs/<tool>/<rule>-<wildcards>.log`).

### 3.4 `docs/testing.md`
Index of the three existing suites:
- `workflow/bgcflow/tests/` — package unit tests for `bgcflow_wrapper`, run via `cd workflow/bgcflow && python -m pytest tests/`
- `tests/plm_novelty/` — PLM side-channel test suite (index build, embedding, query, config schema, no-upstream-DAG-edit guard), run via `pytest tests/plm_novelty`
- `.tests/unit/` — Snakemake rule-script unit tests (arts, fastani, gtdb, antismash summary, etc.), run via `pytest .tests/unit`
Each entry: what it covers, how to run just that suite, and where a new test for that area should go.

### 3.5 `TODOS.md` (repo root)
Flat checkbox task queue, sections: `## In progress`, `## Next`, `## Done` (recent only — old entries get dropped once a worksheet references them, since the worksheet is the durable record).

### 3.6 `docs/worksheets/`
Session/task trace docs. One file per unit of work: `docs/worksheets/YYYY-MM-DD-<topic>.md` with sections **Goal**, **Context** (links to design docs / related worksheets), **Steps taken**, **Current state**, **Next step if handed off**. Updated as work progresses. Committed with the related code changes. At a natural checkpoint, the agent tags the commit `worksheet/<topic>` (local tag only, per this repo's git safety norms — never auto-pushed).

### 3.7 `docs/agent-feedback.md`
Append-only log, one short entry per session: what worked, what didn't, a suggested workflow tweak. Header documents that this file is periodically read by the user in an interactive session to refine workflows, and that agents should append, never rewrite, past entries.

---

## 4. Non-goals for this pass

No automation is being built (no auto-fixing hooks, no `tools/bin/` scripts, no review/loop skills). This pass is pure documentation/process scaffolding. If any of the deferred items later prove necessary, they get their own scoped spec rather than being retrofitted here.
