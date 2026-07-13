# Worksheet: MIBiG embedding resilience (threads, memory cap, retry wrapper)

## Goal

Make the MIBiG ESM2 index build ("`resources/mibig_esm2/index.faiss`") actually seamless to run on a shared multi-user host: survive OOM, survive disconnects, resume automatically, and not endanger other users' processes.

## Context

- `docs/worksheets/2026-07-13-mibig-embed-checkpoint.md` — prior worksheet, added checkpoint/resume to `embed_sequences()`
- `docs/systems/plm-novelty.md` — living reference, "Running it" section rewritten by this worksheet's work

## Steps taken (chronological, across live debugging in chat)

1. First real run (no `threads:` on the rule) sat at 99.9% CPU (~1 core) for 25+ minutes with zero checkpoint progress. Root-caused: Snakemake defaults a rule with no `threads:` directive to 1, setting `OMP_NUM_THREADS=1`/`MKL_NUM_THREADS=1`/etc. for the job regardless of `--cores` or host core count. Fixed: added `threads: 8` to `build_mibig_esm2_index` and `threads: 4` to `embed_bgc_core_esm2` in `workflow/rules/plm_novelty.smk`.
2. Restarted with `--rerun-triggers mtime` (a plain restart would have also triggered an unrelated, expensive `antismash_db_setup` rebuild due to stale software-environment provenance metadata — unrelated to our change).
3. Confirmed via `/proc/<pid>/environ` that the fix worked (`OMP_NUM_THREADS=8` etc.), and CPU usage jumped to ~658-691%.
4. ~1h45m later, the tmux server itself was gone. `journalctl` showed the kernel OOM killer had SIGKILLed the tmux scope at 12:51:43 — RSS had climbed to ~28GB after only ~500-800 of 39,992 proteins, a growth pattern disproportionate to the ~100MB the finished vectors alone would need. Root cause not fully diagnosed (suspected CPU allocator fragmentation from many variable-length-sequence batches over thousands of iterations) — mitigated instead of fixed, since fully diagnosing a PyTorch/glibc allocator fragmentation issue was out of proportion to the actual goal (get the one-time index built).
5. Confirmed no lingering processes/zombies after the OOM event; found swap left fully exhausted (8Gi/8Gi) even after RAM recovered — user cleared it manually with `sudo swapoff -a && sudo swapon -a` (needed a real password, out of my reach).
6. Built the actual resilience layer:
   - `checkpoint_every` exposed as a CLI flag (`--checkpoint-every`) and `rule_parameters.plm.checkpoint_every` config option (default 10 batches ≈ 80 proteins, down from the previous hardcoded 50/400 — smaller interval means less rework lost per kill, since kills are expected to recur until the underlying memory-growth issue is properly root-caused). TDD'd in `test_build_mibig_esm2_index.py`.
   - Added `checkpoint_every` to `workflow/schemas/config.schema.yaml` (`rule_parameters.plm` has `additionalProperties: false` — would have failed schema validation otherwise).
   - `tools/run_mibig_embed.sh` — new script (first in a new `tools/` dir, previously deferred in the agentic-infra design). Loops running the Snakemake target inside `systemd-run --user --scope -p MemoryMax=16G -p MemorySwapMax=0` (no sudo needed — confirmed via `systemctl --user status` and a live test scope), retrying with a short delay on any failure until the target file exists or retries are exhausted. Smoke-tested the retry/resume loop logic in isolation (`set -e` + `if cmd; then` semantics) before trusting it on the real job.
7. Discovered `test_config_schema.py` fails even in the `bgcflow` env (which has `jsonschema`) with a real `SchemaError` — `'type': ['string', 'None']` isn't valid JSON Schema (`'None'` should be `'null'`), on an unrelated field (`projects.items.prokka-db`). Confirmed pre-existing via `git stash`. Logged in `TODOS.md`, not fixed (out of scope, unrelated field).
8. Rewrote `docs/systems/plm-novelty.md`'s "Running it" section to document all of the above as the living reference — the wrapper script, the threads fix, the memory situation, and the CPU-only correction below.
9. Corrected an earlier wrong claim in this doc: `device: auto` will **never** use a GPU in the `plm_novelty` conda env, because it deliberately pins `pytorch-cpu`/`faiss-cpu` — confirmed via `torch.version.cuda is None` in that env. This matters even though the host does have a GPU, because that GPU is also ~95% used by another user's `llama-server` process (46GB total, 43.9GB used) — even if the env were CUDA-enabled, the ~1.5GB free wouldn't reliably fit `esm2_t30_150M_UR50D`'s weights + CUDA context + activation memory.

10. First real launch of the fixed+wrapped job (attempt 1) hit RSS 15.2GB in 2m41s — an order of magnitude *faster* growth than the original single-threaded run (28GB over 1h41m), almost hitting the 16G cap within 3 minutes. Root-caused: glibc's malloc allocates one arena per thread by default, and `threads: 8` (step 1) meant 8x the arenas, each accumulating memory that isn't returned to the OS across many small alloc/free cycles per forward pass — the same underlying growth pattern, just multiplied by thread count. Fixed by exporting `MALLOC_ARENA_MAX=1` in `tools/run_mibig_embed.sh` before invoking `systemd-run`/`snakemake`, capping glibc to a single arena regardless of thread count.
11. Verified live: killed the still-growing attempt (had to kill both the outer wrapper loop and its orphaned `snakemake`/`python` children — a mid-run script edit doesn't affect a bash process already executing, confirmed the hard way), relaunched, confirmed `MALLOC_ARENA_MAX=1` reached the actual process via `/proc/<pid>/environ`. New trajectory: 1.9GB early, up to 9.45GB at 3m36s, then *down* to 6.76GB at 4m02s, settling ~6.3-6.8GB by 4m21s — a bounded sawtooth (rise then reclaim) well under the 16GB cap, versus the prior unbounded climb. Not a full root-cause fix (the sawtooth peaks aren't fully explained — likely protein-length variance across batches), but it took the job from "guaranteed to hit the cap every ~3 minutes" to "plateaus safely."

## Current state

Code, config, schema, tests, and docs changes are complete and passing (29/29 in `tests/plm_novelty` excluding the pre-existing `test_config_schema.py` issue). The `MALLOC_ARENA_MAX=1` fix (steps 10-11) landed live, during monitoring, after the first commit/tag for this worksheet — it gets its own follow-up commit/tag (`worksheet/2026-07-13-mibig-embed-resilience-arena-fix`) rather than amending the first one. Job is running in `tmux -t plm_index`, detached, memory plateaued safely as of the last check.

## Next step if handed off

Check `git tag -l 'worksheet/*'` for both `worksheet/2026-07-13-mibig-embed-resilience` and `worksheet/2026-07-13-mibig-embed-resilience-arena-fix` — if both present, all code/docs work in this worksheet is committed. Then check `tmux ls` for a `plm_index` session and `ls resources/mibig_esm2/.embed_checkpoint/*.npz | wc -l` for live progress — if the tmux session is gone AND `resources/mibig_esm2/index.faiss` doesn't exist, the wrapper's retry loop (up to 100 attempts by default) may have been interrupted by something outside its own retry logic (e.g. the whole tmux server dying again, or a `MAX_RETRIES` exhaustion) — restart with `tmux new -s plm_index` → `tools/run_mibig_embed.sh` → detach.

Still-unsolved and explicitly out of scope for this worksheet: the actual root cause of the remaining sawtooth memory pattern (worth a real profiling pass — `tracemalloc` or `memory_profiler` across a few hundred proteins — if the 16GB cap still gets hit periodically over a full multi-hour run).
