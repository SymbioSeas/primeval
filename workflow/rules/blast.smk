rule run_blast:
    input:
        oligos="resources/oligos/all_oligos.fasta",
        db_nhr="resources/blast_db/{accession}/{accession}.nhr",
    output:
        tsv=blast_output(config["results_dir"] + "/blast/{accession}.tsv"),
    params:
        db="resources/blast_db/{accession}/{accession}",
        evalue=config["blast_evalue"],
        perc_identity=config["blast_perc_identity"],
        word_size=config["blast_word_size"],
    resources:
        mem_mb=4000,
    log:
        config["results_dir"] + "/logs/{accession}.blast.log",
    shell:
        """
        blastn \
            -task blastn-short \
            -query "{input.oligos}" \
            -db "{params.db}" \
            -evalue {params.evalue} \
            -perc_identity {params.perc_identity} \
            -word_size {params.word_size} \
            -strand both \
            -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq" \
            -out "{output.tsv}" 2> "{log}"
        """
