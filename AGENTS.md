# AGENTS.md

Router for any agent (Claude, Codex, Cursor, or otherwise) working in this repo. Start here, then go to `CLAUDE.md` for architecture and workflow details — this file doesn't duplicate that content, it points to it.

## What this is

`bgcflow-plm` — BGCFlow, a Snakemake workflow for BGC analysis across prokaryotic pangenomes, plus an in-progress PLM (protein language model) novelty-retrieval side-channel. Full description: `CLAUDE.md`.

## Where to find things

| Need | Location |
|---|---|
| Architecture, config, running the workflow | `CLAUDE.md` |
| Design rationale for a feature (historical, why-focused) | `docs/*-design.md`, `docs/superpowers/specs/*-design.md` |
| Living reference for a specific system (current, what-focused) | `docs/systems/*.md` |
| Repo-specific coding conventions | `docs/coding-conventions.md` |
| Test suite index — what exists, how to run it, where to add tests | `docs/testing.md` |
| Current task queue | `TODOS.md` |
| Session/task work logs, handoff points | `docs/worksheets/` |
| Cross-session feedback on how agents should work here | `docs/agent-feedback.md` |

## Greppable-header convention

Every doc under `docs/systems/`, plus `docs/coding-conventions.md` and `docs/testing.md`, opens with a short header in its first ~7 lines: title, one-line scope, one-line summary. This lets an agent run `grep -A5 "^# " docs/systems/*.md docs/coding-conventions.md docs/testing.md` to find the right doc without opening every file. Keep headers accurate — if a doc's scope changes, update the header in the same edit.

## System docs are self-healing, not pre-built

`docs/systems/` is not meant to have one file per subsystem from day one. A system doc gets created the first time an agent does non-trivial work in an area that doesn't have one yet, and gets updated (not left stale) by whichever agent next touches that system. If you make a meaningful change to a system that has a doc, update the doc as part of the same change — don't leave it to drift.

## Worksheets: session traces for handoff

For any non-trivial task, create `docs/worksheets/YYYY-MM-DD-<topic>.md` at the start, with sections: **Goal**, **Context** (links to design docs / related worksheets), **Steps taken**, **Current state**, **Next step if handed off**. Update it as you go — the point is that a partial worksheet is enough for a different agent (or a different tool entirely — Codex, Cursor, Claude) to pick up the work cold. Commit the worksheet together with the related code changes.

When a worksheet reaches a checkpoint (task done, or handed off mid-stream), tag the commit `worksheet/<topic>` locally:

```bash
git tag worksheet/2026-07-13-agentic-infra-setup
```

Tags are **local only** — do not push tags unless explicitly asked. Find prior work with `git tag -l 'worksheet/*'`.

## Feedback log

`docs/agent-feedback.md` is an append-only log of what worked and what didn't, one short entry per session. Append to it; don't rewrite past entries. The user periodically reads this to refine how agents should work in this repo.

## Task queue

`TODOS.md` is the flat, current task list. It's not a durable record — once a task is done and has a worksheet, the worksheet is the durable record and the `TODOS.md` line can be dropped.

## Running and testing the app

This is a Snakemake pipeline wrapping external bioinformatics tools — there's no long-running "app" to launch. "Running it" means a targeted `snakemake -n` dry run, a scoped rule run, or the relevant pytest suite (see `docs/testing.md`). Prefer this project's `run` and `verify` skills (if available in your harness) over inventing ad hoc invocations.
