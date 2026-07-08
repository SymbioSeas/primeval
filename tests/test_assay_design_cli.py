"""Unit test for the assay-design CLI entry point on the worked example."""
from pathlib import Path
import assay_design.cli as adc

REPO = Path(__file__).resolve().parents[1]
EXAMPLE = REPO / "assay-design" / "example"


def test_cli_main_runs_example_and_prints_header(tmp_path, capsys):
    out = tmp_path / "out"
    rc = adc.main([
        "--matrix", str(EXAMPLE / "example_gene_presence_absence.csv"),
        "--isolates-dir", str(EXAMPLE / "isolate_groups"),
        "--gene-data", str(EXAMPLE / "gene_data.csv"),
        "--representatives", str(EXAMPLE / "representatives.tsv"),
        "--output-dir", str(out),
    ])
    assert rc == 0
    assert list(out.glob("*_conserved_orthologs.tsv"))
    assert list(out.glob("*_specific_orthologs.tsv"))
    assert any(f.stat().st_size > 0 for f in out.glob("*_conserved_proteins.faa"))
    printed = capsys.readouterr().out
    assert "assay-design  |  dataset:" in printed
    assert "Done. Outputs in" in printed


def test_cli_returns_nonzero_when_no_sequences_extracted(tmp_path):
    out = tmp_path / "out"
    bad_reps = tmp_path / "bad_reps.tsv"
    bad_reps.write_text(
        "group_stem\trepresentative_assembly\n"
        "Vmed_pathogenic\t/some/path/McD53.faa\n"
        "Vmed_protective\t/some/path/PNB23_20_7.faa\n"
        "Vmed_intermediate-avirulent\t/some/path/Fyc26.faa\n"
    )
    rc = adc.main([
        "--matrix", str(EXAMPLE / "example_gene_presence_absence.csv"),
        "--isolates-dir", str(EXAMPLE / "isolate_groups"),
        "--gene-data", str(EXAMPLE / "gene_data.csv"),
        "--representatives", str(bad_reps),
        "--output-dir", str(out),
    ])
    assert rc != 0
