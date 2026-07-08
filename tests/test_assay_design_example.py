"""End-to-end smoke test for the assay-design worked example via the console command."""
import shutil
import subprocess
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
EXAMPLE = REPO / "assay-design" / "example"


def test_example_runs_and_produces_nonempty_outputs(tmp_path):
    exe = shutil.which("assay-design")
    assert exe, "assay-design command not on PATH (is the package installed?)"
    out = tmp_path / "output"
    subprocess.run(
        [exe,
         "--matrix", str(EXAMPLE / "example_gene_presence_absence.csv"),
         "--isolates-dir", str(EXAMPLE / "isolate_groups"),
         "--gene-data", str(EXAMPLE / "gene_data.csv"),
         "--representatives", str(EXAMPLE / "representatives.tsv"),
         "--output-dir", str(out)],
        check=True,
    )
    conserved = list(out.glob("*_conserved_orthologs.tsv"))
    specific = list(out.glob("*_specific_orthologs.tsv"))
    faa = list(out.glob("*_conserved_proteins.faa"))
    assert conserved, "no conserved orthologs written"
    assert specific, "no specific orthologs written"
    assert any(f.stat().st_size > 0 for f in faa), "no non-empty protein FASTA"
    assert any(len(f.read_text().splitlines()) > 1 for f in specific), "no specific ortholog rows"
