"""
Extract representative protein and nucleotide sequences for each inclusion
group's conserved orthologs.

This is the second step of the assay-design workflow, run after
parse_group_specific_orthologs.py has identified conserved orthologs.

For each inclusion group:
  1. Read the group's *_conserved_orthologs.tsv (from the ortholog step)
  2. Identify the representative assembly from the representatives config
  3. Pull locus IDs from the representative assembly column (skipping NaN —
     gene absent in that assembly at the chosen threshold)
  4. Look up sequences in gene_data by matching annotation_id
  5. Write one .faa and one .fna per group

If fewer than --min-match-rate of locus IDs are found in gene_data (e.g. due
to annotation-ID format mismatch between gene_data and the PA matrix), the
group is skipped with a warning rather than erroring. This allows mixed
datasets where gene_data format is consistent for some assemblies but not
others.

Outputs (written into orthologs_dir):
  <group_stem>_conserved_proteins.faa   — protein FASTA
  <group_stem>_conserved_genes.fna      — nucleotide FASTA

Representatives config format (TSV, one row per group):
  group_stem<TAB>representative_assembly
  Vmed_protective<TAB>PNB23_20_7
  Vmed_pathogenic<TAB>McD53
  ...

  group_stem must exactly match the *_conserved_orthologs.tsv filename stem.
  representative_assembly must exactly match a column header in the
  presence/absence matrix (and therefore in the conserved orthologs TSV).

Usage
-----
    python extract_proteins.py \\
        --orthologs-dir FMS2026_orthologs \\
        --gene-data FMS2026_gene_data.csv \\
        --representatives FMS2026_representatives.tsv

    python extract_proteins.py \\
        --orthologs-dir mSystems2025_orthologs \\
        --gene-data mSystems2025_gene_data.csv \\
        --representatives mSystems2025_representatives.tsv \\
        --min-match-rate 0.5
"""

import argparse
import pathlib
import sys

import pandas as pd


def load_representatives(rep_path: pathlib.Path) -> dict:
    """Return {group_stem: representative_assembly} from a 2-column TSV."""
    df = pd.read_csv(rep_path, sep="\t", dtype=str)
    return dict(zip(df["group_stem"].str.strip(), df["representative_assembly"].str.strip()))


def build_header_desc(gene_name, description) -> str:
    parts = []
    if pd.notna(gene_name) and str(gene_name).strip() not in ("", "nan"):
        parts.append(f"gene={str(gene_name).strip()}")
    desc = str(description).strip() if pd.notna(description) else ""
    parts.append(f"product={desc}" if desc and desc != "nan" else "product=hypothetical protein")
    return " ".join(parts)


def write_fasta(records: list, out_path: pathlib.Path) -> None:
    with open(out_path, "w") as f:
        for seq_id, desc, seq in records:
            header = f">{seq_id} {desc}" if desc else f">{seq_id}"
            f.write(header + "\n")
            for i in range(0, len(seq), 60):
                f.write(seq[i : i + 60] + "\n")


def _extract_sequences(df, rep_assembly, gene_lookup, orthologs_dir, group_stem,
                       category, min_match_rate):
    """Write representative protein/nucleotide FASTAs for one group and category
    ("conserved" or "specific"). Returns True if files were written.

    `df` is the already-loaded orthologs table (conserved or specific); the
    representative column is assumed validated by the caller.
    """
    if rep_assembly not in df.columns:
        return False

    rep_loci = df[rep_assembly].dropna().tolist()
    if not rep_loci:
        print(f"  SKIP [{group_stem} / {category}]: representative '{rep_assembly}' "
              f"has no non-empty locus IDs")
        return False

    matched_ids = [lid for lid in rep_loci if lid in gene_lookup.index]
    match_rate = len(matched_ids) / len(rep_loci)
    if match_rate < min_match_rate:
        print(
            f"  SKIP [{group_stem} / {category}]: {len(matched_ids)}/{len(rep_loci)} locus IDs "
            f"({match_rate:.0%}) matched in gene_data for '{rep_assembly}' — "
            f"annotation ID format likely differs from PA matrix",
            file=sys.stderr,
        )
        return False

    n_missing = len(rep_loci) - len(matched_ids)
    if n_missing:
        print(f"  NOTE [{group_stem} / {category}]: {n_missing}/{len(rep_loci)} locus IDs "
              f"not found in gene_data (likely Panaroo refound genes — excluded from output)")

    prot_records, dna_records = [], []
    for locus_id in matched_ids:
        row = gene_lookup.loc[locus_id]
        # gene_lookup.loc may return a DataFrame if annotation_id is duplicated
        if isinstance(row, pd.DataFrame):
            row = row.iloc[0]
        desc = build_header_desc(row.get("gene_name"), row.get("description"))
        prot_seq = str(row.get("prot_sequence", ""))
        dna_seq = str(row.get("dna_sequence", ""))
        if prot_seq and prot_seq != "nan":
            prot_records.append((locus_id, desc, prot_seq))
        if dna_seq and dna_seq != "nan":
            dna_records.append((locus_id, desc, dna_seq))

    faa_file = orthologs_dir / f"{group_stem}_{category}_proteins.faa"
    fna_file = orthologs_dir / f"{group_stem}_{category}_genes.fna"
    write_fasta(prot_records, faa_file)
    write_fasta(dna_records, fna_file)
    print(
        f"  [{group_stem} / {category}]: {len(prot_records)} proteins / {len(dna_records)} genes "
        f"from '{rep_assembly}' ({match_rate:.0%} match) → {faa_file.name}, {fna_file.name}"
    )
    return True


def process_dataset(
    orthologs_dir,
    gene_data_path,
    representatives_path,
    min_match_rate: float = 0.5,
):
    orthologs_dir = pathlib.Path(orthologs_dir)
    gene_data_path = pathlib.Path(gene_data_path)
    representatives_path = pathlib.Path(representatives_path)

    print(f"\n  Orthologs dir: {orthologs_dir}")
    print(f"  Gene data:     {gene_data_path}")
    print(f"  Representatives: {representatives_path}")

    for path, label in [
        (gene_data_path, "gene_data"),
        (representatives_path, "representatives"),
        (orthologs_dir, "orthologs_dir"),
    ]:
        if not path.exists():
            print(f"  ERROR: {label} not found: {path}", file=sys.stderr)
            return 0

    representatives = load_representatives(representatives_path)
    if not representatives:
        print(f"  ERROR: no entries in {representatives_path}", file=sys.stderr)
        return 0

    conserved_files = sorted(orthologs_dir.glob("*_conserved_orthologs.tsv"))
    if not conserved_files:
        print(
            f"  ERROR: no *_conserved_orthologs.tsv files in {orthologs_dir}\n"
            f"  Run Stage 1 (parse_group_specific_orthologs.py) first.",
            file=sys.stderr,
        )
        return 0

    # Pre-filter gene_data to only the representative assemblies needed,
    # avoiding loading the full (potentially multi-GB) file into memory.
    needed_assemblies = set(representatives.values())
    print(f"  Loading gene_data (filtering to {len(needed_assemblies)} representative assemblies)...")

    gene_data = pd.read_csv(
        gene_data_path,
        usecols=["gff_file", "annotation_id", "prot_sequence", "dna_sequence", "gene_name", "description"],
        dtype=str,
        low_memory=False,
    )
    available_gff = sorted(gene_data["gff_file"].dropna().unique())
    gene_data = gene_data[gene_data["gff_file"].isin(needed_assemblies)].copy()
    print(f"  {len(gene_data)} gene_data rows for {len(needed_assemblies)} representative assemblies")

    if gene_data.empty:
        sample = ", ".join(available_gff[:5]) if available_gff else "(none)"
        print(
            f"  ERROR: none of the representative names {sorted(needed_assemblies)} match a "
            f"'gff_file' value in {gene_data_path.name}.\n"
            f"  representative_assembly must be a genome name (a 'gff_file' in gene_data, "
            f"which is also a matrix column) — not a file path or .faa filename.\n"
            f"  Valid genome names include: {sample}",
            file=sys.stderr,
        )

    # Build lookup: annotation_id → row (for fast access)
    gene_lookup = gene_data.set_index("annotation_id")

    n_written = 0

    for conserved_file in conserved_files:
        group_stem = conserved_file.name.replace("_conserved_orthologs.tsv", "")

        if group_stem not in representatives:
            print(f"  SKIP [{group_stem}]: not listed in {representatives_path.name}")
            continue

        rep_assembly = representatives[group_stem]

        conserved_df = pd.read_csv(conserved_file, sep="\t", dtype=str)

        if rep_assembly not in conserved_df.columns:
            hint = ""
            if "/" in rep_assembly or rep_assembly.endswith((".faa", ".fna", ".fasta")):
                hint = " (this looks like a file path or filename — use the genome name instead)"
            sample = ", ".join(list(conserved_df.columns)[-3:])
            print(
                f"  SKIP [{group_stem}]: representative '{rep_assembly}' is not a genome "
                f"column in {conserved_file.name}{hint}. It must be a genome name matching a "
                f"matrix column, e.g. one of: {sample}",
                file=sys.stderr,
            )
            continue

        # Extract representative sequences for the conserved orthologs …
        if _extract_sequences(conserved_df, rep_assembly, gene_lookup, orthologs_dir,
                              group_stem, "conserved", min_match_rate):
            n_written += 1

        # … and for the group-specific subset (the clade-specific candidate targets).
        specific_file = orthologs_dir / f"{group_stem}_specific_orthologs.tsv"
        if specific_file.exists():
            specific_df = pd.read_csv(specific_file, sep="\t", dtype=str)
            if _extract_sequences(specific_df, rep_assembly, gene_lookup, orthologs_dir,
                                  group_stem, "specific", min_match_rate):
                n_written += 1
        else:
            print(f"  NOTE [{group_stem}]: no {specific_file.name} found; "
                  f"skipping specific-ortholog extraction")

    if n_written == 0:
        print(
            "  No output files written. If all groups were skipped due to low match rate,\n"
            "  check that gene_data annotation_ids match locus IDs in the PA matrix."
        )

    return n_written


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Extract representative protein (.faa) and nucleotide (.fna) sequences "
            "for each inclusion group's conserved orthologs."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--orthologs-dir",
        required=True,
        type=pathlib.Path,
        metavar="DIR",
        help="Directory containing *_conserved_orthologs.tsv files (Stage 1 output)",
    )
    parser.add_argument(
        "--gene-data",
        required=True,
        type=pathlib.Path,
        metavar="FILE",
        help="gene_data CSV with columns: gff_file, annotation_id, prot_sequence, dna_sequence, gene_name, description",
    )
    parser.add_argument(
        "--representatives",
        required=True,
        type=pathlib.Path,
        metavar="FILE",
        help="TSV with columns group_stem and representative_assembly",
    )
    parser.add_argument(
        "--min-match-rate",
        type=float,
        default=0.5,
        metavar="FLOAT",
        help=(
            "Minimum fraction of locus IDs that must match gene_data annotation_ids "
            "for a group to be processed. Groups below this rate are skipped with a "
            "warning. Default: 0.5 (50%%)."
        ),
    )
    args = parser.parse_args()
    process_dataset(
        args.orthologs_dir,
        args.gene_data,
        args.representatives,
        args.min_match_rate,
    )


if __name__ == "__main__":
    main()
