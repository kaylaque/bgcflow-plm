# PLM Novelty Retrieval

Scope: living reference for the ESM2/FAISS novelty-retrieval side-channel.
Summary: additive Snakemake lane that embeds antiSMASH biosynthetic-core proteins with frozen ESM2 and retrieves nearest MIBiG neighbors via FAISS, without ever modifying antiSMASH's own outputs. For the original rationale/design tradeoffs, see `docs/design-plm-retrieval-bgcflow.md`.

---

## What it does

Reads antiSMASH region GBKs and the MIBiG `knownclusterblast` directory (already on disk after `antismash_db_setup`), extracts biosynthetic-core proteins, embeds them with a frozen ESM2 model, and queries a prebuilt FAISS index built from MIBiG 4.0 core proteins. Output is a novelty/similarity ranking per BGC region — it never overwrites or reinterprets any antiSMASH call. It's a retrieval side-channel, not a replacement classifier.

## Where it lives

- Rule file: `workflow/rules/plm_novelty.smk`
- Scripts: `workflow/bgcflow/bgcflow/features/build_mibig_esm2_index.py`, `workflow/bgcflow/bgcflow/features/embed_esm2.py`, `workflow/bgcflow/bgcflow/features/plm_novelty_query.py`
- Conda env: `workflow/envs/plm_novelty.yaml`
- Config: `rule_parameters.plm` in `config/config.yaml` (`model`, `batch_size`, `device`, `checkpoint_every`, `max_tokens_per_batch`, `top_k`, `distance_threshold`); toggled on/off via `rules.plm-novelty` (currently `FALSE` by default)

## Outputs

- `resources/mibig_esm2/{index.faiss, metadata.parquet, model_version.txt}` — the MIBiG index, built once and cached; only rebuilds when `antismash_db_setup`'s `knownclusterblast` output refreshes
- `data/interim/plm_novelty/` and `data/processed/{name}/plm_novelty/` — per-run outputs (never written into antiSMASH's own output directories)

## Version gating

`plm_novelty.smk` requires antiSMASH v7+. If `antismash_major_version < 7`, the file registers no rules and emits a stderr warning instead — antiSMASH v6's `antismash_db_setup` doesn't expose a `knownclusterblast/4.0` directory to build the index from. This check lives at the top of the rule file, not inside individual rules.

## Running it: use `tools/run_mibig_embed.sh`, not a bare `snakemake` call

The recommended way to (re)build `resources/mibig_esm2/index.faiss` is:

```bash
tmux new -s plm_index
tools/run_mibig_embed.sh
# Ctrl-b d to detach; tmux attach -t plm_index to check on it later
```

Don't just run `snakemake ... resources/mibig_esm2/index.faiss` directly for this target — see why below. `tools/run_mibig_embed.sh` wraps that same command with the resilience this job actually needs on a shared multi-user host.

**Why it's not a bare `snakemake` call.** Building the full index means embedding all 39,992 MIBiG core proteins on CPU (see the CPU-only note below), which takes hours even when healthy. Four distinct failure modes showed up across real attempts, in the order they were found:

1. **No progress checkpointing** (fixed): `embed_sequences()` in `embed_esm2.py` now checkpoints — pass `checkpoint_dir` (which `build_index()` in `build_mibig_esm2_index.py` sets automatically, defaulting to `<output_dir>/.embed_checkpoint`) and it atomically writes a shard every `checkpoint_every` batches (config: `rule_parameters.plm.checkpoint_every`, default 10 batches ≈ 80 proteins — deliberately small, since kills were expected to recur). Re-running the same command resumes from the last shard rather than starting over. The checkpoint is scratch state, deleted once the build finishes successfully. A checkpoint built with a different `model` is detected via `manifest.json` and discarded rather than silently mixed in.
2. **No `threads:` directive** (fixed): the rule had none, so Snakemake defaulted it to 1 and set `OMP_NUM_THREADS=1`/`MKL_NUM_THREADS=1`/etc. for the job — throttling PyTorch's CPU matmuls to a single core regardless of `--cores` or how many cores the host has. It now declares `threads: 8` (`embed_bgc_core_esm2`, the per-genome rule, got `threads: 4`). If you ever see a PLM rule pegged at ~100% CPU on a many-core box, check for this first.
3. **glibc malloc-arena bloat from more threads** (mitigated): raising to 8 threads made memory growth *worse*, not better — glibc allocates one malloc arena per thread by default, and none of them return memory to the OS across the many small alloc/free cycles a batched transformer forward pass does. RSS hit 15GB in under 3 minutes at one point. `tools/run_mibig_embed.sh` exports `MALLOC_ARENA_MAX=1`, capping glibc to one arena regardless of thread count — turned an unbounded climb into a bounded sawtooth.
4. **Outlier-length proteins forcing catastrophic batch padding** (fixed, root cause): MIBiG contains megasynthase (PKS/NRPS) proteins up to 5,145 aa against a dataset mean of 628 aa (max overall: 18,447 aa). `embed_sequences()` originally batched by *count* (`batch_size`), and ESM pads every sequence in a batch to the batch's longest member — a batch containing one such outlier forced all `batch_size` sequences through a multi-thousand-token forward pass, and since self-attention memory scales with the *square* of sequence length, this could spike memory by many GB in a single batch. Because checkpointed resumes always restart at the same "remaining" position, this reproduced *deterministically* every retry (~57-60s to OOM, every time) once the resume point reached one of these clusters — not a slow leak, a wall. Fixed with `_length_bucketed_batches()`: sequences are now grouped by length, bounded by both `batch_size` and a `max_tokens_per_batch` token budget (config: `rule_parameters.plm.max_tokens_per_batch`, default `batch_size * 1024`) — short sequences still batch together at full `batch_size`, while outliers get progressively smaller (down to solo) batches. A single sequence longer than the whole budget still gets processed alone, never skipped.

On a shared multi-user machine, items 2-3 also matter beyond just "our job is slow": an uncontained OOM event risks the kernel picking a victim process belonging to a *different* user, not just ours. `tools/run_mibig_embed.sh` runs the job inside a `systemd-run --user --scope -p MemoryMax=... -p MemorySwapMax=0` cgroup (no sudo needed) so a runaway job only kills *itself*, cleanly, and the wrapper's retry loop resumes it from the last checkpoint. Tune `MEMORY_MAX` (default `16G`) and `MAX_RETRIES` (default `100`) via env vars if needed — with item 4 fixed, this cap should now rarely if ever trigger, but it stays as defense-in-depth.

**Restarting manually** (bypassing the wrapper, e.g. for debugging) must include `--rerun-triggers mtime`:
```bash
snakemake --use-conda --cores 8 --rerun-triggers mtime -- resources/mibig_esm2/index.faiss
```
Without it, Snakemake may also decide to rebuild `antismash_db_setup` (a large, expensive external database) based on stale software-environment provenance metadata — unrelated to anything you actually changed. This targets the file directly, which works even when `rules.plm-novelty` is `FALSE` — that toggle only controls `rule all`'s default target set, not what you can request explicitly.

**CPU vs GPU:** this repo's `plm_novelty` conda env pins `pytorch-cpu`/`faiss-cpu` deliberately (`workflow/envs/plm_novelty.yaml`) — a genuinely CPU-only PyTorch build (`torch.version.cuda` is `None`), so `device: auto` will **never** use a GPU in this env regardless of what hardware is available on the host, even though the code's device-resolution logic itself isn't hardcoded. If you need GPU throughput, that requires a separate CUDA-enabled env, not just a config toggle — not built as of this writing.

The ESM2 checkpoint weights (`esm2_t30_150M_UR50D`) are a one-time download cached by `torch.hub` outside the repo (`~/.cache/torch/hub/checkpoints/`) — check there before assuming a slow first run is a re-download.

## Tests

`tests/plm_novelty/` (see `docs/testing.md`) — includes a guard test (`test_no_upstream_dag_edit.py`) that must keep passing for any change here, since the side-channel's core invariant is that it never touches antiSMASH's own outputs.

## Status

Prototype plumbing per the design doc; not yet validated science. Default-off (`plm-novelty: FALSE`).
