"""assay-design command: identify clade-conserved / clade-specific orthologs
from a Panaroo pangenome and extract representative sequences (Stage 1 + Stage 2)."""
import argparse
import pathlib
import sys

from .parse_group_specific_orthologs import process_dataset as run_stage1, derive_dataset_name
from .extract_proteins import process_dataset as run_stage2


def build_parser():
    p = argparse.ArgumentParser(
        prog="assay-design",
        description="Identify group-conserved / group-specific orthologs from a Panaroo "
                    "pangenome and extract representative sequences (Stage 1 + Stage 2).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--matrix", required=True, type=pathlib.Path, metavar="FILE",
                   help="Panaroo gene_presence_absence matrix (.csv or .tsv)")
    p.add_argument("--isolates-dir", required=True, type=pathlib.Path, metavar="DIR",
                   help="Directory of .txt inclusion-group files (one genome per line)")
    p.add_argument("--gene-data", type=pathlib.Path, metavar="FILE", default=None,
                   help="Panaroo gene_data.csv (required unless --skip-stage2)")
    p.add_argument("--representatives", type=pathlib.Path, metavar="FILE", default=None,
                   help="TSV of group_stem<TAB>representative_assembly (required unless --skip-stage2)")
    p.add_argument("--output-dir", type=pathlib.Path, default=None, metavar="DIR",
                   help="Stage 1 output dir. Default: <dataset>_orthologs/ beside --matrix")
    p.add_argument("--threshold", type=float, default=0.9, metavar="FLOAT",
                   help="Conserved: present in >= THRESHOLD of inclusion assemblies. "
                        "Specific: conserved AND present in <= (1 - THRESHOLD) of exclusion. Default: 0.9")
    p.add_argument("--min-match-rate", type=float, default=0.5, metavar="FLOAT",
                   help="Stage 2 skip guard: min fraction of locus IDs that must match "
                        "gene_data annotation_ids to write sequences for a group. Default: 0.5")
    p.add_argument("--skip-stage2", action="store_true",
                   help="Run Stage 1 only (ortholog identification); skip sequence extraction.")
    return p


def _print_header(dataset, matrix, n_groups, output_dir):
    print(f"assay-design  |  dataset: {dataset}")
    print(f"  matrix   : {matrix}")
    print(f"  groups   : {n_groups} inclusion group(s)")
    print(f"  output   : {output_dir}")
    print()


def _print_footer(output_dir, skip_stage2):
    print()
    print(f"Done. Outputs in {output_dir}/")
    print("  <group>_conserved_orthologs.tsv, <group>_specific_orthologs.tsv")
    if not skip_stage2:
        print("  <group>_conserved_proteins.faa, <group>_conserved_genes.fna")


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)

    if not (0.0 < args.threshold <= 1.0):
        parser.error("--threshold must be between 0.0 (exclusive) and 1.0 (inclusive)")
    if not args.skip_stage2 and (args.gene_data is None or args.representatives is None):
        parser.error("--gene-data and --representatives are required unless --skip-stage2 is set")

    dataset = derive_dataset_name(args.matrix)
    output_dir = args.output_dir or (args.matrix.parent / f"{dataset}_orthologs")
    n_groups = len(list(pathlib.Path(args.isolates_dir).glob("*.txt")))
    _print_header(dataset, args.matrix, n_groups, output_dir)

    print(f"{'='*60}\nSTAGE 1 — {dataset}\n{'='*60}")
    run_stage1(args.matrix, args.isolates_dir, output_dir=args.output_dir, threshold=args.threshold)

    n_written = None
    if not args.skip_stage2:
        print(f"\n{'='*60}\nSTAGE 2 — {dataset}\n{'='*60}")
        n_written = run_stage2(
            orthologs_dir=output_dir,
            gene_data_path=args.gene_data,
            representatives_path=args.representatives,
            min_match_rate=args.min_match_rate,
        )

    if n_written == 0:
        print("\nERROR: Stage 2 extracted no sequences. Check that each "
              "representative_assembly in the representatives file is a genome "
              "name matching a matrix column and a gene_data 'gff_file' value "
              "(not a file path or .faa filename).", file=sys.stderr)
        return 1

    _print_footer(output_dir, args.skip_stage2)
    return 0


if __name__ == "__main__":
    sys.exit(main())
