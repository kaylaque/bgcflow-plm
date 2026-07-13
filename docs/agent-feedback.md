# Agent Feedback Log

Append-only. One short entry per session: what worked, what didn't, a suggested workflow tweak. Agents append new entries at the bottom — never edit or delete past entries. The user periodically reads this whole file in an interactive session to refine how agents should work in this repo, so keep entries short and concrete rather than restating context that's already in a worksheet.

---

## 2026-07-13 — Agentic infra setup

Scoping an 18-item checklist down to what actually fits the repo (via brainstorming's clarifying questions) worked well and avoided building dead-weight infra (visual regression, perf benchmarking) that had no fit in a UI-less Snakemake pipeline. Worth keeping: asking about *other* agent tools in use (Codex/Cursor here) before designing docs, since it changes whether docs can assume Claude-specific conventions. No process problems this session.
