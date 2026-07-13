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
- Config: `rule_parameters.plm` in `config/config.yaml` (`model`, `batch_size`, `device`, `checkpoint_every`, `top_k`, `distance_threshold`); toggled on/off via `rules.plm-novelty` (currently `FALSE` by default)

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

**Why it's not a bare `snakemake` call.** Building the full index means embedding all 39,992 MIBiG core proteins on CPU (see the CPU-only note below), which takes many hours. Two failure modes showed up on the very first real attempt:

1. **No progress checkpointing** (fixed): `embed_sequences()` in `embed_esm2.py` now checkpoints — pass `checkpoint_dir` (which `build_index()` in `build_mibig_esm2_index.py` sets automatically, defaulting to `<output_dir>/.embed_checkpoint`) and it atomically writes a shard every `checkpoint_every` batches (config: `rule_parameters.plm.checkpoint_every`, default 10 batches ≈ 80 proteins — deliberately small, since attempt #2 below showed kills can recur). Re-running the same command resumes from the last shard rather than starting over. The checkpoint is scratch state, deleted once the build finishes successfully. A checkpoint built with a different `model` is detected via `manifest.json` and discarded rather than silently mixed in.
2. **Unbounded memory growth → OOM killer** (not fully root-caused, mitigated instead): a real run's RSS climbed to ~28GB after only ~500 of 39,992 proteins, eventually getting the whole job SIGKILLed by the kernel OOM killer. On a *shared* machine, an uncontained OOM event risks the kernel picking a victim process belonging to a different user, not just ours. `tools/run_mibig_embed.sh` runs the job inside a `systemd-run --user --scope -p MemoryMax=... -p MemorySwapMax=0` cgroup (no sudo needed) so a runaway job only kills *itself*, cleanly, and the wrapper's retry loop just resumes it from the last checkpoint. Tune `MEMORY_MAX` (default `16G`) and `MAX_RETRIES` (default `100`) via env vars if needed.

Also fixed: the `build_mibig_esm2_index` rule had no `threads:` directive, so Snakemake defaulted it to 1 and set `OMP_NUM_THREADS=1`/`MKL_NUM_THREADS=1`/etc. for the job — throttling PyTorch's CPU matmuls to a single core regardless of `--cores` or how many cores the host has. It now declares `threads: 8` (`embed_bgc_core_esm2`, the per-genome rule, got `threads: 4`). If you ever see a PLM rule pegged at ~100% CPU on a many-core box, check for this first.

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
