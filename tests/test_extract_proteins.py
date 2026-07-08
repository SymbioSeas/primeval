"""Tests for Stage 2 (extract_proteins) failure signalling."""
import assay_design.extract_proteins as ep


def _setup(tmp_path, rep_value):
    orthologs = tmp_path / "orth"
    orthologs.mkdir()
    (orthologs / "G_conserved_orthologs.tsv").write_text(
        "Gene\trepA\ngeneX\tlocus1\n"
    )
    gene_data = tmp_path / "gene_data.csv"
    gene_data.write_text(
        "gff_file,scaffold_name,clustering_id,annotation_id,prot_sequence,dna_sequence,gene_name,description\n"
        "repA,s,c,locus1,MPROT,ATGC,geneX,desc\n"
    )
    reps = tmp_path / "reps.tsv"
    reps.write_text(f"group_stem\trepresentative_assembly\nG\t{rep_value}\n")
    return orthologs, gene_data, reps


def test_process_dataset_returns_count_written_on_success(tmp_path):
    orthologs, gene_data, reps = _setup(tmp_path, "repA")
    n = ep.process_dataset(orthologs, gene_data, reps)
    assert n == 1
    assert (orthologs / "G_conserved_proteins.faa").exists()


def test_process_dataset_returns_zero_when_representative_is_a_path(tmp_path):
    orthologs, gene_data, reps = _setup(tmp_path, "/some/path/repA.faa")
    n = ep.process_dataset(orthologs, gene_data, reps)
    assert n == 0
