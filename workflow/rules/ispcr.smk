rule run_ispcr:
    input:
        blast=config["results_dir"] + "/blast/{accession}.tsv",
        fna=config["assembly_dir"] + "/{accession}.fna",
        assay_table=config["assay_table"],
    output:
        detection=config["results_dir"] + "/amplicons/{accession}.csv",
    params:
        max_primer_mismatches=config["max_primer_mismatches"],
        prime3_exact_nt=config["prime3_exact_nt"],
        max_probe_mismatches=config["max_probe_mismatches"],
        max_amplicon_size=config["max_amplicon_size"],
        store_amplicon_sequences=config["store_amplicon_sequences"],
    resources:
        mem_mb=2000,
    shell:
        """
        python {SCRIPTS}/run_ispcr.py \
            --blast "{input.blast}" \
            --fna "{input.fna}" \
            --assay-table "{input.assay_table}" \
            --max-primer-mismatches {params.max_primer_mismatches} \
            --prime3-exact-nt {params.prime3_exact_nt} \
            --max-probe-mismatches {params.max_probe_mismatches} \
            --max-amplicon-size {params.max_amplicon_size} \
            --store-amplicon-sequences {params.store_amplicon_sequences} \
            --detection-out "{output.detection}"
        """
