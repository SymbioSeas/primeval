"""
Run the assay-design ortholog pipeline on a single dataset.

Stage 1  parse_group_specific_orthologs.py
           Identify group-conserved and group-specific orthologs for each
           inclusion group defined by a .txt file. Writes
           <stem>_conserved_orthologs.tsv and <stem>_specific_orthologs.tsv.

Stage 2  extract_proteins.py
           For each group, extract representative protein (.faa) and nucleotide
           (.fna) sequences from gene_data.csv using the locus IDs in the
           representative-assembly column.

Usage
-----
    python run_pipeline.py \\
        --matrix        gene_presence_absence.csv \\
        --isolates-dir  isolate_groups \\
        --gene-data     gene_data.csv \\
        --representatives representatives.tsv

    python run_pipeline.py ... --threshold 1.0      # strict all-or-nothing
    python run_pipeline.py ... --skip-stage2        # ortholog identification only
"""

import argparse
import pathlib

from parse_group_specific_orthologs import process_dataset as run_stage1, derive_dataset_name
from extract_proteins import process_dataset as run_stage2


def main():
    parser = argparse.ArgumentParser(
        description="Run the assay-design ortholog pipeline (Stage 1 + Stage 2) on one dataset.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--matrix", required=True, type=pathlib.Path, metavar="FILE",
                        help="Panaroo gene_presence_absence matrix (.csv or .tsv)")
    parser.add_argument("--isolates-dir", required=True, type=pathlib.Path, metavar="DIR",
                        help="Directory of .txt inclusion-group files (one genome per line)")
    parser.add_argument("--gene-data", type=pathlib.Path, metavar="FILE", default=None,
                        help="Panaroo gene_data.csv (required unless --skip-stage2)")
    parser.add_argument("--representatives", type=pathlib.Path, metavar="FILE", default=None,
                        help="TSV of group_stem<TAB>representative_assembly (required unless --skip-stage2)")
    parser.add_argument("--output-dir", type=pathlib.Path, default=None, metavar="DIR",
                        help="Stage 1 output dir. Default: <dataset>_orthologs/ beside --matrix")
    parser.add_argument("--threshold", type=float, default=0.9, metavar="FLOAT",
                        help="Conserved: present in >= THRESHOLD of inclusion assemblies. "
                             "Specific: conserved AND present in <= (1 - THRESHOLD) of exclusion. Default: 0.9")
    parser.add_argument("--min-match-rate", type=float, default=0.5, metavar="FLOAT",
                        help="Stage 2 skip guard: min fraction of locus IDs that must match "
                             "gene_data annotation_ids to write sequences for a group. Default: 0.5")
    parser.add_argument("--skip-stage2", action="store_true",
                        help="Run Stage 1 only (ortholog identification); skip sequence extraction.")
    args = parser.parse_args()

    if not (0.0 < args.threshold <= 1.0):
        parser.error("--threshold must be between 0.0 (exclusive) and 1.0 (inclusive)")
    if not args.skip_stage2 and (args.gene_data is None or args.representatives is None):
        parser.error("--gene-data and --representatives are required unless --skip-stage2 is set")

    dataset_name = derive_dataset_name(args.matrix)

    print(f"\n{'='*60}\nSTAGE 1 — {dataset_name}\n{'='*60}")
    run_stage1(args.matrix, args.isolates_dir, output_dir=args.output_dir, threshold=args.threshold)

    if not args.skip_stage2:
        orthologs_dir = args.output_dir or (args.matrix.parent / f"{dataset_name}_orthologs")
        print(f"\n{'='*60}\nSTAGE 2 — {dataset_name}\n{'='*60}")
        run_stage2(
            orthologs_dir=orthologs_dir,
            gene_data_path=args.gene_data,
            representatives_path=args.representatives,
            min_match_rate=args.min_match_rate,
        )

    print("\nDone.")


if __name__ == "__main__":
    main()
