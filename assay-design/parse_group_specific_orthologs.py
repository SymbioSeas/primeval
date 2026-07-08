"""
Parse a Panaroo gene_presence_absence matrix to identify group-conserved and
group-specific orthologs for each inclusion group defined by .txt files.

Definitions (governed by --threshold, default 0.9)
---------------------------------------------------
Conserved   : orthogroup present in >= threshold fraction of inclusion assemblies.
Specific    : conserved AND present in <= (1 - threshold) fraction of exclusion
              assemblies (i.e., effectively absent from the exclusion set).

At the default threshold of 0.9:
  conserved  = present in >= 90% of inclusion assemblies
  specific   = conserved AND present in <= 10% of exclusion assemblies

At threshold 1.0 these collapse to the strict all-or-nothing definitions.

For every inclusion .txt file, two output TSVs are written:

  <stem>_conserved_orthologs.tsv
      All conserved orthogroups. Specificity metric columns are inserted between
      the metadata block and the inclusion genome locus-ID columns:
        n_inclusion_present      – count of inclusion assemblies where present
        frac_inclusion_present   – n_inclusion_present / total inclusion assemblies
        n_exclusion_present      – count of exclusion assemblies where present
        frac_exclusion_present   – n_exclusion_present / total exclusion assemblies
        exclusion_assemblies_present – semicolon-delimited list of those assemblies

  <stem>_specific_orthologs.tsv
      Subset where frac_exclusion_present <= (1 - threshold).
      Specificity columns omitted (values are at or below the threshold for all rows).

Handles both TSV and CSV Panaroo output. Metadata columns (Gene, Annotation,
etc.) are detected automatically and carried through to output.

Usage
-----
    python parse_group_specific_orthologs.py \\
        --matrix mSystems2025_gene_presence_absence_roary.tsv \\
        --isolates-dir mSystems2025_isolates_txt

    python parse_group_specific_orthologs.py \\
        --matrix FMS2026_gene_presence_absence_roary.csv \\
        --isolates-dir FMS2026_isolates_txt \\
        --threshold 0.9

    python parse_group_specific_orthologs.py \\
        --matrix FMS2026_gene_presence_absence_roary.csv \\
        --isolates-dir FMS2026_isolates_txt \\
        --threshold 1.0   # strict: all inclusion present, all exclusion absent
"""

import argparse
import pathlib
import re
import sys

import numpy as np
import pandas as pd

# All column names used by Panaroo/Roary as non-genome metadata
PANAROO_METADATA_COLS = {
    "Gene",
    "Non-unique Gene name",
    "Annotation",
    "No. isolates",
    "No. sequences",
    "Avg sequences per isolate",
    "Genome Fragment",
    "Order within Fragment",
    "Accessory Fragment",
    "Accessory Order with Fragment",
    "QC",
    "Min group size nuc",
    "Max group size nuc",
    "Avg group size nuc",
}

METRIC_COLS = [
    "n_inclusion_present",
    "frac_inclusion_present",
    "n_exclusion_present",
    "frac_exclusion_present",
    "exclusion_assemblies_present",
]


def read_matrix(matrix_path: pathlib.Path) -> pd.DataFrame:
    sep = "\t" if matrix_path.suffix.lower() == ".tsv" else ","
    # encoding='utf-8-sig' strips a leading BOM if present. Panaroo/Roary CSVs
    # exported on Windows/Excel often carry one, which would otherwise turn the
    # first metadata column ("Gene") into a BOM-prefixed name and misclassify it
    # as a genome column — corrupting exclusion counts and dropping the annotation.
    return pd.read_csv(matrix_path, sep=sep, low_memory=False, encoding="utf-8-sig")


def split_columns(df: pd.DataFrame):
    """Return (metadata_cols, genome_cols) by matching against known Panaroo names."""
    metadata_cols = [col for col in df.columns if col in PANAROO_METADATA_COLS]
    genome_cols = [col for col in df.columns if col not in PANAROO_METADATA_COLS]
    return metadata_cols, genome_cols


def read_isolate_list(txt_path: pathlib.Path):
    # utf-8-sig tolerates a BOM on isolate-list files saved from Excel/Windows.
    return [line.strip() for line in txt_path.read_text(encoding="utf-8-sig").splitlines() if line.strip()]


def derive_dataset_name(matrix_path: pathlib.Path) -> str:
    """Extract leading dataset label from matrix filename (e.g. 'mSystems2025')."""
    match = re.match(r"^([^_]+)", matrix_path.stem)
    return match.group(1) if match else matrix_path.stem


def add_specificity_metrics(
    conserved_df: pd.DataFrame, inclusion: list, exclusion: list
) -> pd.DataFrame:
    """
    Append inclusion and exclusion coverage metric columns to conserved_df.
    Returns a copy; does not modify in place.
    """
    conserved_df = conserved_df.copy()

    # Inclusion metrics
    inc_notna = conserved_df[inclusion].notna()
    n_inc = inc_notna.sum(axis=1).astype(int)
    frac_inc = (n_inc / len(inclusion)).round(4)
    conserved_df["n_inclusion_present"] = n_inc
    conserved_df["frac_inclusion_present"] = frac_inc

    # Exclusion metrics
    if exclusion:
        excl_notna = conserved_df[exclusion].notna()
        n_excl = excl_notna.sum(axis=1).astype(int)
        frac_excl = (n_excl / len(exclusion)).round(4)
        excl_names = np.array(exclusion)
        excl_bool_arr = excl_notna.to_numpy()
        excl_list_series = pd.Series(
            [";".join(excl_names[row]) for row in excl_bool_arr],
            index=conserved_df.index,
        )
    else:
        n_excl = pd.Series(0, index=conserved_df.index)
        frac_excl = pd.Series(0.0, index=conserved_df.index)
        excl_list_series = pd.Series("", index=conserved_df.index)

    conserved_df["n_exclusion_present"] = n_excl
    conserved_df["frac_exclusion_present"] = frac_excl
    conserved_df["exclusion_assemblies_present"] = excl_list_series

    return conserved_df


def process_dataset(matrix_path, isolates_dir, output_dir=None, threshold=0.9):
    """
    Parameters
    ----------
    threshold : float, 0.0–1.0
        Minimum fraction of inclusion assemblies that must carry an orthogroup
        for it to be called conserved. Specificity cutoff is (1 - threshold):
        an orthogroup is called specific only when <= (1 - threshold) of
        exclusion assemblies carry it.
    """
    matrix_path = pathlib.Path(matrix_path)
    isolates_dir = pathlib.Path(isolates_dir)

    dataset_name = derive_dataset_name(matrix_path)

    if output_dir is None:
        output_dir = matrix_path.parent / f"{dataset_name}_orthologs"
    else:
        output_dir = pathlib.Path(output_dir)
    output_dir.mkdir(exist_ok=True)

    incl_pct = f"{threshold * 100:.0f}"
    excl_pct = f"{(1 - threshold) * 100:.0f}"

    print(f"\nDataset: {dataset_name}")
    print(f"  Matrix:       {matrix_path}")
    print(f"  Isolates dir: {isolates_dir}")
    print(f"  Output dir:   {output_dir}")
    print(f"  Threshold:    conserved >= {incl_pct}% inclusion  |  specific <= {excl_pct}% exclusion")

    df = read_matrix(matrix_path)
    metadata_cols, genome_cols = split_columns(df)

    print(f"  Metadata columns ({len(metadata_cols)}): {metadata_cols}")
    print(f"  Genome columns: {len(genome_cols)}")
    print(f"  Orthogroups: {len(df)}")

    txt_files = sorted(isolates_dir.glob("*.txt"))
    if not txt_files:
        print(f"  ERROR: no .txt files found in {isolates_dir}", file=sys.stderr)
        return

    genome_col_set = set(genome_cols)

    for txt_path in txt_files:
        inclusion = read_isolate_list(txt_path)

        missing = [g for g in inclusion if g not in genome_col_set]
        if missing:
            print(
                f"  WARNING [{txt_path.name}]: {len(missing)} genome(s) not found in "
                f"matrix columns — skipping: {missing}",
                file=sys.stderr,
            )
            inclusion = [g for g in inclusion if g in genome_col_set]

        if not inclusion:
            print(f"  SKIP [{txt_path.name}]: no valid genomes after validation.")
            continue

        exclusion = [g for g in genome_cols if g not in set(inclusion)]

        # --- Conserved: >= threshold fraction of inclusion assemblies present ---
        conserved_mask = df[inclusion].notna().mean(axis=1) >= threshold
        conserved_df = add_specificity_metrics(df[conserved_mask], inclusion, exclusion)

        conserved_out_cols = metadata_cols + METRIC_COLS + inclusion
        conserved_file = output_dir / f"{txt_path.stem}_conserved_orthologs.tsv"
        conserved_df[conserved_out_cols].to_csv(conserved_file, sep="\t", index=False)

        # --- Specific: conserved AND <= (1 - threshold) exclusion presence ---
        specific_mask = conserved_df["frac_exclusion_present"] <= (1 - threshold)
        specific_df = conserved_df[specific_mask]
        specific_out_cols = metadata_cols + inclusion
        specific_file = output_dir / f"{txt_path.stem}_specific_orthologs.tsv"
        specific_df[specific_out_cols].to_csv(specific_file, sep="\t", index=False)

        print(
            f"  [{txt_path.stem}]: {len(conserved_df)} conserved, "
            f"{len(specific_df)} specific "
            f"({len(inclusion)} incl. / {len(exclusion)} excl. genomes)"
        )


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Identify group-conserved and group-specific orthologs from a Panaroo "
            "gene_presence_absence matrix."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--matrix",
        required=True,
        type=pathlib.Path,
        metavar="FILE",
        help="Path to gene_presence_absence_roary.tsv or .csv (Panaroo output)",
    )
    parser.add_argument(
        "--isolates-dir",
        required=True,
        type=pathlib.Path,
        metavar="DIR",
        help="Directory containing .txt inclusion-group files (one genome name per line)",
    )
    parser.add_argument(
        "--output-dir",
        type=pathlib.Path,
        default=None,
        metavar="DIR",
        help=(
            "Output directory for results. "
            "Default: <dataset>_orthologs/ in the same directory as --matrix"
        ),
    )
    parser.add_argument(
        "--threshold",
        type=float,
        default=0.9,
        metavar="FLOAT",
        help=(
            "Presence/absence threshold as a fraction (0.0–1.0). "
            "Conserved: orthogroup present in >= THRESHOLD of inclusion assemblies. "
            "Specific: conserved AND present in <= (1 - THRESHOLD) of exclusion assemblies. "
            "Default: 0.9 (90%% inclusion, <=10%% exclusion)."
        ),
    )
    args = parser.parse_args()

    if not (0.0 < args.threshold <= 1.0):
        parser.error("--threshold must be between 0.0 (exclusive) and 1.0 (inclusive)")

    process_dataset(args.matrix, args.isolates_dir, args.output_dir, args.threshold)


if __name__ == "__main__":
    main()
