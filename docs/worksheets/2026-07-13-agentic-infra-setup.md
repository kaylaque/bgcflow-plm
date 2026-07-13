# Worksheet: Agentic workflow infra setup

## Goal

Evaluate the user's 18-item agentic-coding-setup checklist for fit in this repo, then implement the subset that's actually beneficial here (solo/small-team Snakemake bioinformatics pipeline, no UI, Codex/Cursor also used alongside Claude Code).

## Context

- Design doc: `docs/superpowers/specs/2026-07-13-agentic-infra-design.md` (full fit assessment for all 18 items, and exact spec for each file built)
- User picked the "Docs & workflow infra" tier via brainstorming skill: `AGENTS.md`, `docs/systems/`, `docs/coding-conventions.md`, `docs/testing.md`, `TODOS.md`, `docs/worksheets/`, `docs/agent-feedback.md`
- User confirmed Codex/Cursor are also used against this repo — deferred cross-agent review tooling to a later scoped pass rather than building it now

## Steps taken

1. Explored repo state: no `AGENTS.md`, empty `.claude/`, three separate test suites, existing `docs/design-plm-retrieval-bgcflow.md` as the one design-doc precedent
2. Ran brainstorming skill: scoped down from 18 items to a "docs & workflow infra" tier via two rounds of clarifying questions
3. Wrote and committed the design spec (`docs/superpowers/specs/2026-07-13-agentic-infra-design.md`)
4. Built `AGENTS.md` router
5. Built `docs/coding-conventions.md`, `docs/testing.md`, `docs/systems/plm-novelty.md`
6. Built this worksheet (dogfooding the convention it documents)
7. Building `docs/agent-feedback.md` and `TODOS.md` next
8. Will add a pointer from `CLAUDE.md` to `AGENTS.md`, tag the commit `worksheet/2026-07-13-agentic-infra-setup`, commit, and push

## Current state

All planned files built except the final CLAUDE.md pointer + commit/tag/push step, which is in progress in the same session.

## Next step if handed off

If picked up cold: check `git log` for whether the final commit (CLAUDE.md pointer + tag `worksheet/2026-07-13-agentic-infra-setup`) landed. If not, finish step 8 above — everything else is done.
