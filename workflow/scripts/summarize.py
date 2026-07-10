import csv
import hashlib
import subprocess
import argparse
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from datetime import datetime


METADATA_COLS = [
    'accession', 'ani_best_match_organism', 'ani_match_status', 'ani_taxonomy_check',
    'infraspecific_strain', 'bs_strain', 'bs_host', 'bs_isolation_source', 'assembly_level',
]

DETECTION_COLS = [
    'accession', 'assay', 'detection_call', 'n_amplicons', 'multi_amplicon_flag',
    'amplicon_sizes', 'contig_ids', 'fwd_mismatches', 'rev_mismatches',
    'probe_mismatches', 'probe_strand',
]

PER_ASSAY_COLS = [
    'accession', 'infraspecific_strain', 'bs_strain', 'bs_host', 'bs_isolation_source',
    'assembly_level', 'ani_best_match_organism', 'ani_confidence', 'detection_call',
    'n_amplicons', 'multi_amplicon_flag', 'amplicon_sizes', 'contig_ids',
    'fwd_mismatches', 'rev_mismatches', 'probe_mismatches', 'probe_strand',
]

DETECTION_CALL_ORDER = ['Detected', 'Primer Only', 'Not Detected']

SPECIES_SUMMARY_COLS = [
    'species_group', 'assay', 'n_assemblies', 'n_detected',
    'n_primer_only', 'n_not_detected', 'pct_detected', 'pct_detected_or_primer',
    'n_multi_amplicon', 'max_amplicons',
]

ASSAY_SUMMARY_LONG_COLS = [
    'assay', 'species_group', 'n_assemblies', 'n_detected', 'n_primer_only',
    'n_not_detected', 'pct_detected', 'pct_detected_or_primer',
    'n_multi_amplicon', 'max_amplicons',
]


_HIGH_ANI_STATUSES = {'species_match', 'subspecies_match', 'derived_species_match'}


def _compute_ani_confidence(row) -> str:
    if row['ani_match_status'] in _HIGH_ANI_STATUSES and row['ani_taxonomy_check'] == 'OK':
        return 'High'
    if row['ani_match_status'] == 'genus_match' and row['ani_taxonomy_check'] == 'OK':
        return 'Genus'
    return 'Low'


def _effective_species(row) -> str:
    if row['ani_confidence'] == 'High':
        return row['ani_best_match_organism']
    if row['ani_confidence'] == 'Genus':
        return f"{row['ani_best_match_organism']} (genus match)"
    return 'Unclassified (low confidence ANI)'


def load_detection_results(amplicons_dir: str) -> pd.DataFrame:
    """Load per-assembly detection CSVs from a directory, skipping amplicons and empty files."""
    paths = [p for p in Path(amplicons_dir).glob('*.csv') if p.stat().st_size > 0]
    if not paths:
        return pd.DataFrame(columns=DETECTION_COLS)
    return pd.concat([pd.read_csv(p) for p in paths], ignore_index=True)


def join_metadata(det: pd.DataFrame, meta: pd.DataFrame) -> pd.DataFrame:
    """Left-join detection results with assembly metadata.
    Unmatched assemblies get ani_best_match_organism='Unclassified'; other metadata fields='Unknown'.
    Computes ani_confidence (High/Genus/Low) from ani_match_status and ani_taxonomy_check.
    """
    joined = det.merge(meta[METADATA_COLS], on='accession', how='left')
    joined['ani_best_match_organism'] = joined['ani_best_match_organism'].fillna('Unclassified')
    joined['ani_match_status'] = joined['ani_match_status'].fillna('no_match')
    joined['ani_taxonomy_check'] = joined['ani_taxonomy_check'].fillna('unknown')
    for col in ['infraspecific_strain', 'bs_strain', 'bs_host', 'bs_isolation_source', 'assembly_level']:
        joined[col] = joined[col].fillna('Unknown')
    joined['ani_confidence'] = joined.apply(_compute_ani_confidence, axis=1)
    return joined


def build_species_summary(joined: pd.DataFrame) -> pd.DataFrame:
    """Compute per-species-group × per-assay detection counts and percentages.
    Uses effective_species grouping (High-confidence: raw species; Genus: species + suffix; Low: Unclassified).
    """
    df = joined.copy()
    df['species_group'] = df.apply(_effective_species, axis=1)
    rows = []
    for (species, assay), grp in df.groupby(['species_group', 'assay']):
        n = len(grp)
        n_det = (grp['detection_call'] == 'Detected').sum()
        n_po = (grp['detection_call'] == 'Primer Only').sum()
        n_nd = (grp['detection_call'] == 'Not Detected').sum()
        namp = pd.to_numeric(grp['n_amplicons'], errors='coerce').fillna(0)
        rows.append({
            'species_group': species,
            'assay': assay,
            'n_assemblies': n,
            'n_detected': int(n_det),
            'n_primer_only': int(n_po),
            'n_not_detected': int(n_nd),
            'pct_detected': round(100 * n_det / n, 2) if n > 0 else 0.0,
            'pct_detected_or_primer': round(100 * (n_det + n_po) / n, 2) if n > 0 else 0.0,
            'n_multi_amplicon': int((namp > 1).sum()),
            'max_amplicons': int(namp.max()) if len(namp) else 0,
        })
    return pd.DataFrame(rows, columns=SPECIES_SUMMARY_COLS)


def build_assay_summary_long(species_summary: pd.DataFrame) -> pd.DataFrame:
    """Per assay × species_group, ALL species groups (including pct_detected=0),
    sorted by assay then pct_detected descending."""
    if species_summary.empty:
        return pd.DataFrame(columns=ASSAY_SUMMARY_LONG_COLS)
    df = species_summary[ASSAY_SUMMARY_LONG_COLS].copy()
    return df.sort_values(['assay', 'pct_detected'], ascending=[True, False]).reset_index(drop=True)


def build_species_matrix(species_summary: pd.DataFrame) -> pd.DataFrame:
    """Wide species × assay matrix of pct_detected with a leading n_assemblies column.
    Rows = species_group (as a column after reset), columns = assays."""
    if species_summary.empty:
        return pd.DataFrame(columns=['species_group', 'n_assemblies'])
    n_assemblies = species_summary.groupby('species_group')['n_assemblies'].first()
    matrix = species_summary.pivot_table(
        index='species_group', columns='assay', values='pct_detected', aggfunc='first'
    )
    matrix.insert(0, 'n_assemblies', n_assemblies)
    return matrix.reset_index()


def build_detection_matrix(joined: pd.DataFrame) -> pd.DataFrame:
    """Binary detection matrix: rows=accession, cols=assay, 1=Detected else 0."""
    dups = joined.duplicated(subset=['accession', 'assay'])
    if dups.any():
        raise ValueError(
            f"Duplicate accession×assay rows before pivot: "
            f"{joined[dups][['accession', 'assay']].values.tolist()}"
        )
    verbose = joined.pivot_table(
        index='accession', columns='assay', values='detection_call', aggfunc='first'
    )
    return verbose.map(lambda x: 1 if x == 'Detected' else 0)


def build_detection_by_assembly(binary: pd.DataFrame, meta: pd.DataFrame) -> pd.DataFrame:
    """One row per analyzed accession: all metadata columns + one 0/1 column per assay."""
    return meta.merge(binary.reset_index(), on='accession', how='right')


def write_species_heatmap(species_matrix: pd.DataFrame, figures_dir: str) -> None:
    """Species × assay heatmap colored by pct_detected (0-100)."""
    Path(figures_dir).mkdir(parents=True, exist_ok=True)
    if species_matrix.empty:
        return
    mat = species_matrix.set_index('species_group').drop(columns=['n_assemblies'], errors='ignore')
    if mat.empty:
        return
    annotate = mat.shape[0] * mat.shape[1] <= 300
    fig, ax = plt.subplots(figsize=(max(6, len(mat.columns)), max(4, len(mat) * 0.4)))
    sns.heatmap(mat, ax=ax, cmap='viridis', vmin=0, vmax=100,
                linewidths=0.5, linecolor='white',
                annot=annotate, fmt='.0f',
                cbar_kws={'label': '% detected'})
    ax.set_title('Species detection (% of assemblies detected)')
    ax.set_xlabel('Assay')
    ax.set_ylabel('Species group')
    plt.tight_layout()
    for ext in ('pdf', 'png'):
        fig.savefig(f"{figures_dir}/species_detection_heatmap.{ext}", dpi=150)
    plt.close(fig)


def _md5(path: str) -> str:
    h = hashlib.md5()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(65536), b''):
            h.update(chunk)
    return h.hexdigest()


def _tool_version(cmd: list) -> str:
    try:
        return subprocess.check_output(cmd, stderr=subprocess.STDOUT,
                                       text=True).strip().split('\n')[0]
    except Exception:
        return 'unknown'


def write_run_manifest(manifest_path: str, params: dict, assay_table: str) -> None:
    """Write run manifest with parameters, tool versions, and checksums."""
    lines = [
        f"# Run manifest — generated {datetime.now().isoformat()}",
        "",
        "## Parameters",
    ]
    for k, v in params.items():
        lines.append(f"  {k}: {v}")
    lines += [
        "",
        "## Tool versions",
        f"  BLAST: {_tool_version(['blastn', '-version'])}",
        f"  Python: {_tool_version(['python', '--version'])}",
        "",
        "## Checksums",
        f"  assay_table.csv  MD5: {_md5(assay_table)}",
    ]
    Path(manifest_path).write_text('\n'.join(lines) + '\n')


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--amplicons-dir', required=True)
    p.add_argument('--metadata', required=True)
    p.add_argument('--assay-table', required=True)
    p.add_argument('--reports-dir', required=True)
    p.add_argument('--max-primer-mismatches', type=int, default=2)
    p.add_argument('--prime3-exact-nt', type=int, default=3)
    p.add_argument('--max-probe-mismatches', type=int, default=1)
    p.add_argument('--max-amplicon-size', type=int, default=500)
    p.add_argument('--store-amplicon-sequences',
                   type=lambda x: x.lower() == 'true', default=True)
    p.add_argument('--keep-blast', type=lambda x: x.lower() == 'true', default=False)
    p.add_argument('--keep-logs', type=lambda x: x.lower() == 'true', default=False)
    args = p.parse_args()

    reports = Path(args.reports_dir)
    figures = reports / 'figures'
    per_assay = reports / 'per_assay'
    per_assay.mkdir(parents=True, exist_ok=True)
    figures.mkdir(parents=True, exist_ok=True)

    det = load_detection_results(args.amplicons_dir)
    meta = pd.read_csv(args.metadata, encoding='utf-8-sig')
    joined = join_metadata(det, meta)

    # Per-assay tables
    for assay, grp in joined.groupby('assay'):
        sorted_grp = grp.assign(
            _call_order=pd.Categorical(
                grp['detection_call'], categories=DETECTION_CALL_ORDER, ordered=True
            )
        ).sort_values('_call_order')[PER_ASSAY_COLS]
        sorted_grp.to_csv(per_assay / f"{assay}_results.csv", index=False)

    # Species summary — wide species × assay matrix of pct_detected
    species_summary = build_species_summary(joined)
    build_species_matrix(species_summary).to_csv(reports / 'species_summary.csv', index=False)

    # Assay summary (long) — one row per assay × species_group, all groups incl. misses
    assay_long = build_assay_summary_long(species_summary)
    assay_long.to_csv(reports / 'assay_summary_long.csv', index=False)
    with pd.ExcelWriter(reports / 'assay_summary.xlsx', engine='openpyxl') as writer:
        for assay, grp in assay_long.groupby('assay'):
            grp.to_excel(writer, sheet_name=str(assay)[:31], index=False)

    # Per-assembly detection joined with input metadata
    binary = build_detection_matrix(joined)
    build_detection_by_assembly(binary, meta).to_csv(reports / 'detection_by_assembly.csv', index=False)

    # Species-level heatmap
    write_species_heatmap(build_species_matrix(species_summary), str(figures))

    # Run manifest
    params = {
        'max_primer_mismatches': args.max_primer_mismatches,
        'prime3_exact_nt': args.prime3_exact_nt,
        'max_probe_mismatches': args.max_probe_mismatches,
        'max_amplicon_size': args.max_amplicon_size,
        'store_amplicon_sequences': args.store_amplicon_sequences,
        'keep_blast': args.keep_blast,
        'keep_logs': args.keep_logs,
    }
    write_run_manifest(str(reports / 'run_manifest.txt'), params, args.assay_table)

    print(f"Reports written to {args.reports_dir}")


if __name__ == '__main__':
    main()
