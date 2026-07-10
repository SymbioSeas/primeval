import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent / "workflow" / "scripts"))

import pytest
import pandas as pd
import numpy as np
from summarize import (
    load_detection_results,
    join_metadata,
    build_species_summary,
    build_species_matrix,
    build_assay_summary_long,
    build_detection_matrix,
    build_detection_by_assembly,
    write_species_heatmap,
    write_run_manifest,
)


# --- Fixtures ---

def make_detection_df(accessions, assays, call='Detected', n_amplicons=1):
    """Build a minimal detection DataFrame (like output of run_ispcr)."""
    rows = []
    for acc in accessions:
        for assay in assays:
            rows.append({
                'accession': acc, 'assay': assay, 'detection_call': call,
                'n_amplicons': n_amplicons, 'multi_amplicon_flag': n_amplicons > 1,
                'amplicon_sizes': '150', 'contig_ids': 'c1',
                'fwd_mismatches': 0, 'rev_mismatches': 0,
                'probe_mismatches': 0, 'probe_strand': '+',
            })
    return pd.DataFrame(rows)


def make_metadata_df(accessions, species='Vibrio harveyi',
                     ani_match_status='species_match', ani_taxonomy_check='OK'):
    rows = []
    for acc in accessions:
        rows.append({
            'accession': acc,
            'ani_best_match_organism': species,
            'ani_match_status': ani_match_status,
            'ani_taxonomy_check': ani_taxonomy_check,
            'infraspecific_strain': 'strain_X',
            'bs_strain': 'bs_X',
            'bs_host': 'shrimp',
            'bs_isolation_source': 'water',
            'assembly_level': 'Complete Genome',
        })
    return pd.DataFrame(rows)


# --- load_detection_results ---

def test_load_detection_results(tmp_path):
    """Load multiple per-assembly CSVs into one DataFrame."""
    assays = ['VhPath', 'Valg']
    for i, acc in enumerate(['GCF_001', 'GCF_002']):
        df = make_detection_df([acc], assays)
        df.to_csv(tmp_path / f"{acc}.csv", index=False)
    combined = load_detection_results(str(tmp_path))
    assert len(combined) == 4  # 2 assemblies × 2 assays
    assert set(combined['accession']) == {'GCF_001', 'GCF_002'}


def test_load_detection_results_empty_dir(tmp_path):
    """Empty directory returns empty DataFrame."""
    df = load_detection_results(str(tmp_path))
    assert df.empty


def test_load_detection_results_skips_empty_files(tmp_path):
    """Zero-byte CSV files are skipped without error."""
    (tmp_path / "empty.csv").write_text("")
    df = make_detection_df(['GCF_001'], ['VhPath'])
    df.to_csv(tmp_path / "GCF_001.csv", index=False)
    combined = load_detection_results(str(tmp_path))
    assert len(combined) == 1
    assert combined.iloc[0]['accession'] == 'GCF_001'


# --- join_metadata ---

def test_join_metadata_basic():
    det = make_detection_df(['GCF_001'], ['VhPath'])
    meta = make_metadata_df(['GCF_001'])
    result = join_metadata(det, meta)
    assert 'ani_best_match_organism' in result.columns
    assert result.iloc[0]['ani_best_match_organism'] == 'Vibrio harveyi'
    assert result.iloc[0]['bs_host'] == 'shrimp'
    assert result.iloc[0]['ani_confidence'] == 'High'


def test_join_metadata_unclassified():
    """Assemblies with no metadata match get 'Unclassified' species and 'Unknown' other fields."""
    det = make_detection_df(['GCF_999'], ['VhPath'])
    meta = make_metadata_df(['GCF_001'])  # GCF_999 not in metadata
    result = join_metadata(det, meta)
    assert result.iloc[0]['ani_best_match_organism'] == 'Unclassified'
    assert result.iloc[0]['bs_host'] == 'Unknown'
    assert result.iloc[0]['assembly_level'] == 'Unknown'
    assert result.iloc[0]['ani_confidence'] == 'Low'


# --- build_assay_summary_long / build_species_matrix ---

def test_build_assay_summary_long_includes_misses_and_sorts():
    det = pd.concat([
        make_detection_df(['GCF_001', 'GCF_002'], ['VhPath'], call='Detected'),
        make_detection_df(['GCF_003'], ['VhPath'], call='Not Detected', n_amplicons=0),
    ])
    meta = pd.concat([
        make_metadata_df(['GCF_001', 'GCF_002'], species='Vibrio harveyi'),
        make_metadata_df(['GCF_003'], species='Vibrio sp.'),
    ])
    joined = join_metadata(det, meta)
    summary = build_species_summary(joined)
    long = build_assay_summary_long(summary)
    # Miss (pct_detected=0) is INCLUDED
    assert (long['pct_detected'] == 0).any()
    assert list(long.columns) == [
        'assay', 'species_group', 'n_assemblies', 'n_detected', 'n_primer_only',
        'n_not_detected', 'pct_detected', 'pct_detected_or_primer',
        'n_multi_amplicon', 'max_amplicons',
    ]
    # Sorted by assay, then pct_detected desc (hit before miss within VhPath)
    vh = long[long['assay'] == 'VhPath']['pct_detected'].tolist()
    assert vh == sorted(vh, reverse=True)


def test_build_species_matrix_is_wide():
    det = pd.concat([
        make_detection_df(['GCF_001'], ['VhPath'], call='Detected'),
        make_detection_df(['GCF_001'], ['Valg'], call='Not Detected', n_amplicons=0),
    ])
    meta = make_metadata_df(['GCF_001'])
    joined = join_metadata(det, meta)
    matrix = build_species_matrix(build_species_summary(joined))
    assert 'species_group' in matrix.columns
    assert 'n_assemblies' in matrix.columns
    assert 'VhPath' in matrix.columns and 'Valg' in matrix.columns
    r = matrix[matrix['species_group'] == 'Vibrio harveyi'].iloc[0]
    assert r['VhPath'] == 100.0 and r['Valg'] == 0.0


# --- build_species_summary ---

def test_build_species_summary_counts():
    det = pd.concat([
        make_detection_df(['GCF_001', 'GCF_002'], ['VhPath'], call='Detected'),
        make_detection_df(['GCF_003'], ['VhPath'], call='Not Detected'),
        make_detection_df(['GCF_004'], ['VhPath'], call='Primer Only'),
    ])
    meta = make_metadata_df(['GCF_001', 'GCF_002', 'GCF_003', 'GCF_004'])
    joined = join_metadata(det, meta)
    summary = build_species_summary(joined)
    row = summary[(summary['species_group'] == 'Vibrio harveyi') &
                  (summary['assay'] == 'VhPath')].iloc[0]
    assert row['n_assemblies'] == 4
    assert row['n_detected'] == 2
    assert row['n_not_detected'] == 1
    assert row['n_primer_only'] == 1
    assert abs(row['pct_detected'] - 50.0) < 0.01
    assert abs(row['pct_detected_or_primer'] - 75.0) < 0.01
    assert row['n_multi_amplicon'] == 0
    assert row['max_amplicons'] == 1


# --- per-assay sort order ---

def test_per_assay_sort_order(tmp_path):
    """Per-assay CSVs are sorted Detected → Primer Only → Not Detected."""
    import sys
    from pathlib import Path as P
    sys.path.insert(0, str(P(__file__).parent.parent / "workflow" / "scripts"))
    from summarize import DETECTION_CALL_ORDER
    import pandas as pd

    # Build joined df with all three call types in scrambled order
    det = pd.concat([
        make_detection_df(['GCF_003'], ['VhPath'], call='Not Detected'),
        make_detection_df(['GCF_002'], ['VhPath'], call='Primer Only'),
        make_detection_df(['GCF_001'], ['VhPath'], call='Detected'),
    ])
    meta = make_metadata_df(['GCF_001', 'GCF_002', 'GCF_003'])
    joined = join_metadata(det, meta)

    per_assay_dir = tmp_path / "per_assay"
    per_assay_dir.mkdir()
    for assay, grp in joined.groupby('assay'):
        from summarize import PER_ASSAY_COLS
        sorted_grp = grp.assign(
            _call_order=pd.Categorical(
                grp['detection_call'], categories=DETECTION_CALL_ORDER, ordered=True
            )
        ).sort_values('_call_order')[PER_ASSAY_COLS]
        sorted_grp.to_csv(per_assay_dir / f"{assay}_results.csv", index=False)

    result = pd.read_csv(per_assay_dir / "VhPath_results.csv")
    calls = result['detection_call'].tolist()
    assert calls == ['Detected', 'Primer Only', 'Not Detected']


# --- build_detection_matrix ---

def test_build_detection_matrix_binary():
    det = pd.concat([
        make_detection_df(['GCF_001'], ['VhPath'], call='Detected'),
        make_detection_df(['GCF_001'], ['Valg'], call='Not Detected', n_amplicons=0),
    ])
    joined = join_metadata(det, make_metadata_df(['GCF_001']))
    binary = build_detection_matrix(joined)
    assert binary.loc['GCF_001', 'VhPath'] == 1
    assert binary.loc['GCF_001', 'Valg'] == 0


def test_build_detection_matrix_duplicate_raises():
    det = pd.concat([
        make_detection_df(['GCF_001'], ['VhPath'], call='Detected'),
        make_detection_df(['GCF_001'], ['VhPath'], call='Not Detected'),
    ])
    joined = join_metadata(det, make_metadata_df(['GCF_001']))
    with pytest.raises(ValueError, match="Duplicate"):
        build_detection_matrix(joined)


def test_build_detection_by_assembly_joins_metadata():
    det = pd.concat([
        make_detection_df(['GCF_001'], ['VhPath'], call='Detected'),
        make_detection_df(['GCF_001'], ['Valg'], call='Not Detected', n_amplicons=0),
    ])
    meta = make_metadata_df(['GCF_001'])
    joined = join_metadata(det, meta)
    binary = build_detection_matrix(joined)
    out = build_detection_by_assembly(binary, meta)
    assert len(out) == 1
    row = out.iloc[0]
    assert row['accession'] == 'GCF_001'
    assert row['bs_host'] == 'shrimp'         # metadata carried through
    assert row['VhPath'] == 1 and row['Valg'] == 0


def test_write_species_heatmap_creates_files(tmp_path):
    det = pd.concat([
        make_detection_df(['GCF_001'], ['VhPath'], call='Detected'),
        make_detection_df(['GCF_002'], ['VhPath'], call='Not Detected', n_amplicons=0),
        make_detection_df(['GCF_001', 'GCF_002'], ['Valg'], call='Detected'),
    ])
    meta = make_metadata_df(['GCF_001', 'GCF_002'])
    joined = join_metadata(det, meta)
    matrix = build_species_matrix(build_species_summary(joined))
    figures_dir = str(tmp_path / "figures")
    write_species_heatmap(matrix, figures_dir)
    assert (tmp_path / "figures" / "species_detection_heatmap.pdf").exists()
    assert (tmp_path / "figures" / "species_detection_heatmap.png").exists()


# --- write_run_manifest ---

def test_write_run_manifest(tmp_path):
    params = {
        'max_primer_mismatches': 2,
        'prime3_exact_nt': 3,
        'max_probe_mismatches': 1,
        'max_amplicon_size': 500,
        'store_amplicon_sequences': True,
    }
    assay_table_path = tmp_path / "assay_table.csv"
    assay_table_path.write_text("assay,probe,fwd,rev\n")
    write_run_manifest(str(tmp_path / "manifest.txt"), params, str(assay_table_path))
    content = (tmp_path / "manifest.txt").read_text()
    assert 'max_primer_mismatches: 2' in content
    assert 'prime3_exact_nt: 3' in content
    assert 'assay_table.csv' in content
    assert 'md5:' in content.lower() or 'MD5' in content
    assert 'blast' in content.lower() or 'BLAST' in content


# --- ANI confidence tiers ---

def test_ani_confidence_high():
    """species_match + OK → High."""
    det = make_detection_df(['GCF_001'], ['VhPath'])
    meta = make_metadata_df(['GCF_001'], ani_match_status='species_match', ani_taxonomy_check='OK')
    result = join_metadata(det, meta)
    assert result.iloc[0]['ani_confidence'] == 'High'


def test_ani_confidence_genus():
    """genus_match + OK → Genus."""
    det = make_detection_df(['GCF_001'], ['VhPath'])
    meta = make_metadata_df(['GCF_001'], ani_match_status='genus_match', ani_taxonomy_check='OK')
    result = join_metadata(det, meta)
    assert result.iloc[0]['ani_confidence'] == 'Genus'


def test_ani_confidence_low_mismatch():
    """mismatch + Inconclusive → Low."""
    det = make_detection_df(['GCF_001'], ['VhPath'])
    meta = make_metadata_df(['GCF_001'], ani_match_status='mismatch', ani_taxonomy_check='Inconclusive')
    result = join_metadata(det, meta)
    assert result.iloc[0]['ani_confidence'] == 'Low'


def test_species_summary_separates_genus_match():
    """Genus-confidence assemblies appear under 'species (genus match)', not the plain species name."""
    det = pd.concat([
        make_detection_df(['GCF_001'], ['VhPath'], call='Detected'),
        make_detection_df(['GCF_002'], ['VhPath'], call='Detected'),
    ])
    meta = pd.concat([
        make_metadata_df(['GCF_001'], ani_match_status='species_match', ani_taxonomy_check='OK'),
        make_metadata_df(['GCF_002'], ani_match_status='genus_match', ani_taxonomy_check='OK'),
    ])
    joined = join_metadata(det, meta)
    summary = build_species_summary(joined)
    groups = set(summary['species_group'])
    assert 'Vibrio harveyi' in groups
    assert 'Vibrio harveyi (genus match)' in groups
    # Each should count only its own assemblies
    high = summary[summary['species_group'] == 'Vibrio harveyi'].iloc[0]
    genus = summary[summary['species_group'] == 'Vibrio harveyi (genus match)'].iloc[0]
    assert high['n_assemblies'] == 1
    assert genus['n_assemblies'] == 1
