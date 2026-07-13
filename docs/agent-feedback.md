# Agent Feedback Log

Append-only. One short entry per session: what worked, what didn't, a suggested workflow tweak. Agents append new entries at the bottom — never edit or delete past entries. The user periodically reads this whole file in an interactive session to refine how agents should work in this repo, so keep entries short and concrete rather than restating context that's already in a worksheet.

---

## 2026-07-13 — Agentic infra setup

Scoping an 18-item checklist down to what actually fits the repo (via brainstorming's clarifying questions) worked well and avoided building dead-weight infra (visual regression, perf benchmarking) that had no fit in a UI-less Snakemake pipeline. Worth keeping: asking about *other* agent tools in use (Codex/Cursor here) before designing docs, since it changes whether docs can assume Claude-specific conventions. No process problems this session.

## 2026-07-13 — MIBiG embedding checkpoint

The PLM test suite (`tests/plm_novelty/`) had never actually been run successfully since a refactor a few commits back — three independent mock/fixture bugs, invisible until someone tried to run it in the *correct* conda env. Worth remembering: this repo's PLM tests need the Snakemake-managed `plm_novelty` conda env (`.snakemake/conda/<hash>_`, matched by grepping `conda-meta/` for `fair-esm`/`faiss`/`pytorch`), not the `bgcflow` env — `bgcflow` lacks torch/faiss entirely and silently gives misleading errors. Before trusting a "tests pass" claim in this repo, check which env they actually ran in. Also worth a standing suspicion: a test suite added in one commit and never touched again after the code it tests was refactored is a strong "actually run this before believing it" signal — grep the git log for the test file's last-touched commit vs. the code file's, and if the code moved on without it, run the tests before relying on them.
