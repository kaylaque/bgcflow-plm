"""
PLM novelty-retrieval side-channel for BGCFlow.

Additive lane: reads antiSMASH region GBKs and the MIBiG knownclusterblast
directory (already on disk after antismash_db_setup), embeds biosynthetic-core
proteins with frozen ESM2, and queries a prebuilt FAISS index.

Outputs are written to data/interim/plm_novelty/ and
data/processed/{name}/plm_novelty/ — the lane never writes into the
antiSMASH output directories.

v6 note: antismash_db_setup for v6 does not expose a knownclusterblast/4.0
directory.  If antismash_major_version < 7, this file registers no rules and
emits a warning so the pipeline continues cleanly.
"""

if antismash_major_version < 7:
    print(
        " - WARNING: plm-novelty requires antiSMASH v7+ (current: "
        f"v{antismash_major_version}). Rules not registered.",
        file=sys.stderr,
    )
else:
    rule build_mibig_esm2_index:
        """
        Embed MIBiG 4.0 biosynthetic-core proteins with frozen ESM2 and build
        a FAISS IndexFlatIP (cosine via L2-normalised inner product).

        Runs once; cached in resources/mibig_esm2/.  Re-runs only when the
        knownclusterblast/4.0 directory is refreshed by antismash_db_setup.
        """
        input:
            knownclusterblast=rules.antismash_db_setup.output.knownclusterblast,
        output:
            index=   "resources/mibig_esm2/index.faiss",
            metadata="resources/mibig_esm2/metadata.parquet",
            model_version="resources/mibig_esm2/model_version.txt",
        conda:
            "../envs/plm_novelty.yaml"
        threads: 8
        log:
            "logs/plm_novelty/build_mibig_esm2_index.log",
        params:
            model=config.get("rule_parameters", {}).get("plm", {}).get(
                "model", "esm2_t30_150M_UR50D"
            ),
            batch_size=config.get("rule_parameters", {}).get("plm", {}).get(
                "batch_size", 8
            ),
            device=config.get("rule_parameters", {}).get("plm", {}).get(
                "device", "auto"
            ),
            checkpoint_every=config.get("rule_parameters", {}).get("plm", {}).get(
                "checkpoint_every", 10
            ),
            output_dir="resources/mibig_esm2",
        shell:
            """
            python workflow/bgcflow/bgcflow/features/build_mibig_esm2_index.py \
                --mibig-dir {input.knownclusterblast} \
                --output-dir {params.output_dir} \
                --model {params.model} \
                --batch-size {params.batch_size} \
                --device {params.device} \
                --checkpoint-every {params.checkpoint_every} \
                2>> {log}
            """

    rule extract_bgc_core_proteins:
        """
        Extract biosynthetic-core CDS translations from antiSMASH region GBKs
        for one strain.  Produces a FASTA (may be empty for strains with no
        core proteins — downstream handles this gracefully).
        """
        input:
            gbk="data/interim/antismash/{version}/{strains}/{strains}.gbk",
        output:
            faa="data/interim/plm_novelty/{version}/{strains}/core_proteins.faa",
        params:
            antismash_dir="data/interim/antismash/{version}/{strains}/",
        log:
            "logs/plm_novelty/extract_core/{version}/{strains}.log",
        shell:
            """
            python workflow/bgcflow/bgcflow/features/extract_bgc_core_proteins.py \
                --input {params.antismash_dir} \
                --strains {wildcards.strains} \
                --output {output.faa} \
                2>> {log}
            """

    rule embed_bgc_core_esm2:
        """
        Embed per-strain biosynthetic-core proteins with frozen ESM2.
        Saves L2-normalised float32 vectors to an .npz (ids + vectors keys).
        """
        input:
            faa="data/interim/plm_novelty/{version}/{strains}/core_proteins.faa",
        output:
            npz="data/interim/plm_novelty/{version}/{strains}/embeddings.npz",
        conda:
            "../envs/plm_novelty.yaml"
        threads: 4
        log:
            "logs/plm_novelty/embed/{version}/{strains}.log",
        params:
            model=config.get("rule_parameters", {}).get("plm", {}).get(
                "model", "esm2_t30_150M_UR50D"
            ),
            batch_size=config.get("rule_parameters", {}).get("plm", {}).get(
                "batch_size", 8
            ),
            device=config.get("rule_parameters", {}).get("plm", {}).get(
                "device", "auto"
            ),
        shell:
            """
            python workflow/bgcflow/bgcflow/features/embed_esm2.py \
                --fasta {input.faa} \
                --output {output.npz} \
                --model {params.model} \
                --batch-size {params.batch_size} \
                --device {params.device} \
                2>> {log}
            """

    rule plm_novelty_query:
        """
        Query the MIBiG FAISS index with all per-strain embeddings for a
        project and emit two novelty TSVs:
          - bgc_novelty_retrieval.tsv  (per protein)
          - region_novelty_summary.tsv (per antiSMASH region)
        """
        input:
            index="resources/mibig_esm2/index.faiss",
            metadata="resources/mibig_esm2/metadata.parquet",
            embeddings=lambda wildcards: expand(
                "data/interim/plm_novelty/{version}/{strains}/embeddings.npz",
                version=dependency_version["antismash"],
                strains=PEP_PROJECTS[wildcards.name].sample_table.genome_id.unique(),
            ),
        output:
            per_protein="data/processed/{name}/plm_novelty/bgc_novelty_retrieval.tsv",
            per_region="data/processed/{name}/plm_novelty/region_novelty_summary.tsv",
        conda:
            "../envs/plm_novelty.yaml"
        log:
            "logs/plm_novelty/query/{name}.log",
        params:
            model=config.get("rule_parameters", {}).get("plm", {}).get(
                "model", "esm2_t30_150M_UR50D"
            ),
            top_k=config.get("rule_parameters", {}).get("plm", {}).get("top_k", 10),
            distance_threshold=config.get("rule_parameters", {}).get("plm", {}).get(
                "distance_threshold", 0.5
            ),
        shell:
            """
            python workflow/bgcflow/bgcflow/features/plm_novelty_query.py \
                --embeddings {input.embeddings} \
                --index {input.index} \
                --meta {input.metadata} \
                --top-k {params.top_k} \
                --distance-threshold {params.distance_threshold} \
                --expected-model {params.model} \
                --out-per-protein {output.per_protein} \
                --out-per-region {output.per_region} \
                2>> {log}
            """
