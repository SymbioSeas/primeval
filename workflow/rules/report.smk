rule aggregate_report:
    input:
        detections=expand(
            config["results_dir"] + "/amplicons/{accession}.csv",
            accession=glob_wildcards(config["assembly_dir"] + "/{accession}.fna").accession
        ),
        metadata=config["metadata"],
        assay_table=config["assay_table"],
    output:
        matrices=expand(config["results_dir"] + "/reports/{gcol}_detection_matrix.csv", gcol=GROUP_COLS),
        heatmaps=expand(config["results_dir"] + "/reports/figures/{gcol}_detection_heatmap.pdf", gcol=GROUP_COLS),
        detection_long=config["results_dir"] + "/reports/detection_summary_long.csv",
        assay_summary_xlsx=config["results_dir"] + "/reports/assay_summary.xlsx",
        detection_by_assembly=config["results_dir"] + "/reports/detection_by_assembly.csv",
        assay_performance=config["results_dir"] + "/reports/assay_performance.csv",
        manifest=config["results_dir"] + "/reports/run_manifest.txt",
    params:
        amplicons_dir=config["results_dir"] + "/amplicons",
        reports_dir=config["results_dir"] + "/reports",
        group_by=" ".join(GROUP_BY),
        max_primer_mismatches=config["max_primer_mismatches"],
        prime3_exact_nt=config["prime3_exact_nt"],
        max_probe_mismatches=config["max_probe_mismatches"],
        max_amplicon_size=config["max_amplicon_size"],
        store_amplicon_sequences=config["store_amplicon_sequences"],
        keep_blast=config.get("keep_blast", False),
        keep_logs=config.get("keep_logs", False),
    resources:
        mem_mb=16000,
    shell:
        """
        python {SCRIPTS}/summarize.py \
            --amplicons-dir {params.amplicons_dir} \
            --metadata "{input.metadata}" \
            --assay-table "{input.assay_table}" \
            --reports-dir {params.reports_dir} \
            --group-by {params.group_by} \
            --max-primer-mismatches {params.max_primer_mismatches} \
            --prime3-exact-nt {params.prime3_exact_nt} \
            --max-probe-mismatches {params.max_probe_mismatches} \
            --max-amplicon-size {params.max_amplicon_size} \
            --store-amplicon-sequences {params.store_amplicon_sequences} \
            --keep-blast {params.keep_blast} \
            --keep-logs {params.keep_logs}
        """
