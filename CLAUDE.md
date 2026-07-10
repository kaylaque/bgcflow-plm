# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

**BGCFlow** is a Snakemake workflow for systematic analysis of biosynthetic gene clusters (BGCs) across prokaryotic pangenomes. It orchestrates 20+ bioinformatics tools (antiSMASH, Prokka, GTDB-Tk, BiG-SCAPE, etc.) against collections of genomes from NCBI, PATRIC, or custom sources.

The current branch adds a **PLM (protein language model) novelty-retrieval side-channel** — an additive ESM2-based layer that detects divergent BGCs antiSMASH would miss, without modifying the core DAG. See `docs/design-plm-retrieval-bgcflow.md` for the design.

## Environment setup

```bash
# Create the base conda environment
mamba env create -f envs.yaml
conda activate bgcflow

# Disable conda channel priority (required)
conda config --set channel_priority disabled

# Or install the wrapper CLI instead
pip install bgcflow_wrapper==0.6.2
```

Each pipeline step runs in its **own isolated conda env** defined in `workflow/envs/*.yaml`. Snakemake creates and manages these automatically on first run.

## Running the workflow

```bash
# Dry run (shows what would execute)
snakemake --use-conda -n

# Full run (uses conda envs per rule)
snakemake --use-conda --cores <N>

# Via the wrapper CLI
bgcflow run -n          # dry run
bgcflow run             # full run
bgcflow build report    # build interactive HTML report
bgcflow serve --project <ProjectName>
```

## Configuration

The workflow is configured through two things:

1. **`config/config.yaml`** — top-level config referencing projects, resource paths, and enabled rules. Validated against `workflow/schemas/config.schema.yaml`.
2. **PEP project files** — each project has a `project_config.yaml` (PEP format, loaded via `peppy`) with a sample table, optional `prokka-db`, and optional `gtdb-tax`. Samples can come from `ncbi`, `patric`, or `custom` sources with `fna` or `gbk` input types.

Rules to run are toggled as booleans inside the project config (or globally in `config/config.yaml` under `rules:`). The full list of available rule keywords and their final output paths is in `workflow/rules.yaml`, `workflow/rules_bgc.yaml`, and `workflow/rules_ppanggolin.yaml`.

## Workflow architecture

**Entry point:** `workflow/Snakefile` calls `common.smk`, extracts project/sample info, sets wildcards, then includes all rule modules.

**`workflow/rules/common.smk`** is the backbone — it contains:
- `extract_project_information()` — reads config + PEP projects into DataFrames
- `get_final_output()` / `get_project_outputs()` — builds the target file list from enabled rules
- All helper lambda functions for rule I/O (`get_prokka_outputs`, `get_antismash_regions`, etc.)
- `custom_resource_dir()` — manages symlinks for large external databases (antismash_db, eggnog_db, BiG-SCAPE, bigslice, checkm, gtdbtk)

**Data flow:**
```
config/config.yaml + PEP project files
  → common.smk: DataFrames (DF_PROJECTS, DF_SAMPLES) + wildcard constants
  → rule all: expand final outputs across {name} (project) and {strains} (genome_id) wildcards
  → individual .smk rules in workflow/rules/
```

**Intermediate data lives in `data/interim/`** (per-tool, per-genome subdirs). **Final processed outputs go to `data/processed/{project_name}/`**. Raw custom genomes go in `data/raw/fasta/`.

**Wildcard semantics:**
- `{name}` = project name
- `{strains}` = genome_id (any source)
- `{strains_fna}`, `{strains_genbank}` = genome subsets by input type
- `{ncbi}`, `{patric}`, `{custom_fna}`, `{custom_genbank}` = subsets by source
- `{version}` = antiSMASH version (auto-detected from `workflow/envs/antismash.yaml`)

## The bgcflow Python package

`workflow/bgcflow/` is a pip-installable package (`bgcflow_wrapper`) containing:
- `bgcflow/database/` — functions to gather/process tool outputs into CSV/parquet tables
- `bgcflow/features/` — extraction utilities (antiSMASH gene kinds, MMseqs2, clinker, BGC selection)
- `bgcflow/visualization/` — plotting helpers

Tests are in `workflow/bgcflow/tests/`. Run with:
```bash
cd workflow/bgcflow
python -m pytest tests/
```

Lint:
```bash
flake8 bgcflow tests
```

## PLM novelty retrieval (new work)

The design calls for a new Snakemake rule `plm_novelty_retrieval` in `workflow/rules/` that:
- Reads existing antiSMASH region GBKs (does **not** modify upstream DAG)
- Extracts biosynthetic-core proteins, embeds with frozen ESM2
- Queries against a prebuilt MIBiG FAISS index in `resources/`
- Outputs `data/processed/plm_novelty/{strain}/bgc_novelty_retrieval.tsv` and `region_novelty_summary.tsv`

The PLM layer is an **additive evidence/ranking side-channel only** — it never overwrites any antiSMASH call. Config toggle goes under `rule_parameters.plm` in `config/config.yaml`.

## Key conventions

- Each `.smk` rule file in `workflow/rules/` is self-contained for one tool. Large databases are installed as rule outputs in `resources/` and referenced via symlinks if user provides an external path.
- `workflow/schemas/` contains YAML schemas used by `snakemake.utils.validate` to check config and sample tables at startup.
- Log files: `logs/<tool>/<rule>-<wildcards>.log`
- The Snakefile version (`__version__` in `common.smk`) and all tool versions are pinned in `workflow/envs/*.yaml` files. `get_dependency_version()` in `common.smk` auto-reads versions from those env files to fill `{version}` wildcards.
- Only Linux is officially supported (conda/mamba required).
