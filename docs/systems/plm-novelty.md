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
- Config: `rule_parameters.plm` in `config/config.yaml` (`model`, `batch_size`, `device`, `top_k`, `distance_threshold`); toggled on/off via `rules.plm-novelty` (currently `FALSE` by default)

## Outputs

- `resources/mibig_esm2/{index.faiss, metadata.parquet, model_version.txt}` — the MIBiG index, built once and cached; only rebuilds when `antismash_db_setup`'s `knownclusterblast` output refreshes
- `data/interim/plm_novelty/` and `data/processed/{name}/plm_novelty/` — per-run outputs (never written into antiSMASH's own output directories)

## Version gating

`plm_novelty.smk` requires antiSMASH v7+. If `antismash_major_version < 7`, the file registers no rules and emits a stderr warning instead — antiSMASH v6's `antismash_db_setup` doesn't expose a `knownclusterblast/4.0` directory to build the index from. This check lives at the top of the rule file, not inside individual rules.

## Running it (SSH gotcha)

`build_index()` in `build_mibig_esm2_index.py` has **no checkpointing** — it embeds all MIBiG core proteins (currently 39,992) in one pass and only writes `index.faiss`/`metadata.parquet`/`model_version.txt` at the very end. If the process is killed mid-run (e.g. an SSH session drops, which SIGHUPs a foreground job), **all progress is lost silently** — no partial output, no error, no Snakemake lock left behind, it just looks like the target was never built.

Always run it detached so an SSH drop can't kill it:

```bash
tmux new -s plm_index
snakemake --use-conda --cores 8 resources/mibig_esm2/index.faiss
# Ctrl-b d to detach; tmux attach -t plm_index to check on it later
```

This targets `resources/mibig_esm2/index.faiss` directly, which works even when `rules.plm-novelty` is `FALSE` — that toggle only controls `rule all`'s default target set, not what you can request explicitly.

The ESM2 checkpoint itself (`esm2_t30_150M_UR50D`) is a one-time download cached by `torch.hub` outside the repo (`~/.cache/torch/hub/checkpoints/`) — check there before assuming a slow first run is a re-download; `device: auto` in `rule_parameters.plm` will use a GPU if one is present on the host.

## Tests

`tests/plm_novelty/` (see `docs/testing.md`) — includes a guard test (`test_no_upstream_dag_edit.py`) that must keep passing for any change here, since the side-channel's core invariant is that it never touches antiSMASH's own outputs.

## Status

Prototype plumbing per the design doc; not yet validated science. Default-off (`plm-novelty: FALSE`).
