rule prepare_oligos:
    input:
        assay_table=config["assay_table"],
    output:
        fasta="resources/oligos/all_oligos.fasta",
        log="resources/oligos/prep.log",
    shell:
        """
        python {SCRIPTS}/prepare_oligos.py \
            --assay-table "{input.assay_table}" \
            --fasta-out "{output.fasta}" \
            --log-out "{output.log}"
        """


rule make_blast_db:
    input:
        fna=config["assembly_dir"] + "/{accession}.fna",
    output:
        nhr="resources/blast_db/{accession}/{accession}.nhr",
        nin="resources/blast_db/{accession}/{accession}.nin",
        nsq="resources/blast_db/{accession}/{accession}.nsq",
    params:
        db_prefix="resources/blast_db/{accession}/{accession}",
    log:
        config["results_dir"] + "/logs/{accession}.makeblastdb.log",
    shell:
        """
        makeblastdb \
            -in "{input.fna}" \
            -dbtype nucl \
            -out "{params.db_prefix}" \
            -title {wildcards.accession} > "{log}" 2>&1
        """
