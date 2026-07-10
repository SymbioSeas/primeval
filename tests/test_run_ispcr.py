import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent / "workflow" / "scripts"))

import pytest
import pandas as pd
from run_ispcr import (
    check_3prime_exact,
    filter_primer_hits,
    find_valid_amplicons,
    check_probe_in_amplicons,
    call_detection,
    load_blast_results,
    iupac_match,
    count_iupac_mismatches,
    BLAST_COLS,
)


def make_hit(**kwargs) -> dict:
    """BLAST hit dict with sensible defaults (fwd primer, + strand, perfect match)."""
    defaults = {
        'qseqid': 'TestAssay_fwd', 'sseqid': 'contig1',
        'pident': 100.0, 'length': 17, 'mismatch': 0, 'gapopen': 0,
        'qstart': 1, 'qend': 17, 'sstart': 100, 'send': 116,
        'evalue': 0.001, 'bitscore': 32.0,
        'qseq': 'AGCCGAGCGTTACCAGC', 'sseq': 'AGCCGAGCGTTACCAGC',
    }
    defaults.update(kwargs)
    return defaults


# --- check_3prime_exact ---

def test_3prime_perfect_match():
    assert check_3prime_exact('AGCCGAGCGTTACCAGC', 'AGCCGAGCGTTACCAGC') is True


def test_3prime_last_nt_mismatch():
    assert check_3prime_exact('AGCCGAGCGTTACCAGC', 'AGCCGAGCGTTACCAGT') is False


def test_3prime_third_from_end_mismatch():
    assert check_3prime_exact('AGCCGAGCGTTACCAGC', 'AGCCGAGCGTTACCTGC') is False


def test_3prime_gap_fails():
    assert check_3prime_exact('AGCCGAGCGTTACCAGC', 'AGCCGAGCGTTACCA-C') is False


def test_3prime_mismatch_outside_window_passes():
    assert check_3prime_exact('AGCCGAGCGTTACCAGC', 'TGCCGAGCGTTACCAGC') is True


# --- filter_primer_hits ---

def test_filter_primer_perfect_forward():
    hits = pd.DataFrame([make_hit(qend=17, sstart=100, send=116, mismatch=0)])
    result = filter_primer_hits(hits, primer_seq='AGCCGAGCGTTACCAGC', max_mismatch=2, prime3_exact=3)
    assert len(result) == 1


def test_filter_primer_too_many_mismatches():
    # sseq has 3 mismatches vs primer (positions 0, 5, 10 changed)
    hits = pd.DataFrame([make_hit(qend=17, mismatch=3, sstart=100, send=116,
                                   sseq='TGCCGTGCGTTGCCAGC')])
    result = filter_primer_hits(hits, primer_seq='AGCCGAGCGTTACCAGC', max_mismatch=2, prime3_exact=3)
    assert len(result) == 0


def test_filter_primer_partial_alignment_rejected():
    """3' end absent (qend < primer_len) → rejected."""
    hits = pd.DataFrame([make_hit(qend=14, mismatch=0, qseq='AGCCGAGCGTTACC',
                                   sseq='AGCCGAGCGTTACC', sstart=100, send=113)])
    result = filter_primer_hits(hits, primer_seq='AGCCGAGCGTTACCAGC', max_mismatch=2, prime3_exact=3)
    assert len(result) == 0


def test_filter_primer_5prime_truncated_rejected():
    """3' end present but 5' end missing (qstart > 1) → rejected even with zero mismatches.

    This is the partial-alignment false-positive case: only the last 11 of 20 nt aligned,
    matching what NCBI BLAST reports as <100% query cover.
    """
    hits = pd.DataFrame([make_hit(qstart=7, qend=17, length=11, mismatch=0,
                                   qseq='GTTACCAGC', sseq='GTTACCAGC',
                                   sstart=106, send=116)])
    result = filter_primer_hits(hits, primer_seq='AGCCGAGCGTTACCAGC', max_mismatch=2, prime3_exact=3)
    assert len(result) == 0


def test_filter_primer_both_strands_accepted():
    """filter_primer_hits accepts both + and - strand hits; strand pairing is done in find_valid_amplicons."""
    plus_hit = make_hit(qend=17, mismatch=0, sstart=100, send=116)
    minus_hit = make_hit(qend=17, mismatch=0, sstart=116, send=100)
    hits = pd.DataFrame([plus_hit, minus_hit])
    result = filter_primer_hits(hits, primer_seq='AGCCGAGCGTTACCAGC', max_mismatch=2, prime3_exact=3)
    assert len(result) == 2


def test_filter_primer_minus_strand_accepted():
    hits = pd.DataFrame([make_hit(qend=22, mismatch=0, sstart=300, send=279,
                                   qseq='CGAACGCAATGATTCTCTGAGC',
                                   sseq='CGAACGCAATGATTCTCTGAGC')])
    result = filter_primer_hits(hits, primer_seq='CGAACGCAATGATTCTCTGAGC', max_mismatch=2, prime3_exact=3)
    assert len(result) == 1


def test_filter_primer_3prime_mismatch_rejected():
    hits = pd.DataFrame([make_hit(qend=17, mismatch=1,
                                   qseq='AGCCGAGCGTTACCAGT',
                                   sseq='AGCCGAGCGTTACCAGC',
                                   sstart=100, send=116)])
    result = filter_primer_hits(hits, primer_seq='AGCCGAGCGTTACCAGT', max_mismatch=2, prime3_exact=3)
    assert len(result) == 0


# --- find_valid_amplicons ---

def test_find_amplicons_valid():
    fwd = pd.DataFrame([make_hit(sseqid='c1', sstart=100, send=116, mismatch=0)])
    rev = pd.DataFrame([make_hit(qseqid='TestAssay_rev', sseqid='c1',
                                  sstart=300, send=279, mismatch=0)])
    amps = find_valid_amplicons(fwd, rev, max_amplicon_size=500)
    assert len(amps) == 1
    assert amps[0]['amplicon_size'] == 201  # 300 - 100 + 1


def test_find_amplicons_exceeds_size():
    fwd = pd.DataFrame([make_hit(sseqid='c1', sstart=100, send=116, mismatch=0)])
    rev = pd.DataFrame([make_hit(qseqid='TestAssay_rev', sseqid='c1',
                                  sstart=700, send=679, mismatch=0)])
    amps = find_valid_amplicons(fwd, rev, max_amplicon_size=500)
    assert len(amps) == 0


def test_find_amplicons_different_contigs():
    fwd = pd.DataFrame([make_hit(sseqid='c1', sstart=100, send=116, mismatch=0)])
    rev = pd.DataFrame([make_hit(qseqid='TestAssay_rev', sseqid='c2',
                                  sstart=300, send=279, mismatch=0)])
    amps = find_valid_amplicons(fwd, rev, max_amplicon_size=500)
    assert len(amps) == 0


def test_find_amplicons_inverted_rejected():
    fwd = pd.DataFrame([make_hit(sseqid='c1', sstart=300, send=316, mismatch=0)])
    rev = pd.DataFrame([make_hit(qseqid='TestAssay_rev', sseqid='c1',
                                  sstart=150, send=129, mismatch=0)])
    amps = find_valid_amplicons(fwd, rev, max_amplicon_size=500)
    assert len(amps) == 0


def test_find_amplicons_multi():
    fwd = pd.DataFrame([make_hit(sseqid='c1', sstart=100, send=116, mismatch=0)])
    rev = pd.DataFrame([
        make_hit(qseqid='TestAssay_rev', sseqid='c1', sstart=300, send=279, mismatch=0),
        make_hit(qseqid='TestAssay_rev', sseqid='c1', sstart=450, send=429, mismatch=1,
                 qseq='CGAACGCAATGATTCTCTGAGC', sseq='CGAACGCAATGATTCTCTGAGC'),
    ])
    amps = find_valid_amplicons(fwd, rev, max_amplicon_size=500)
    assert len(amps) == 2


def test_find_amplicons_reverse_orientation():
    """Target gene on - strand: fwd hits - strand, rev hits + strand → valid amplicon."""
    # fwd on - strand: sstart=2000, send=1981 (sstart > send)
    fwd = pd.DataFrame([make_hit(sseqid='c1', sstart=2000, send=1981, mismatch=0,
                                  qseq='AGCCGAGCGTTACCAGC', sseq='AGCCGAGCGTTACCAGC')])
    # rev on + strand: sstart=1720, send=1741 (sstart < send)
    rev = pd.DataFrame([make_hit(qseqid='TestAssay_rev', sseqid='c1',
                                  sstart=1720, send=1741, mismatch=0, qend=22,
                                  qseq='CGAACGCAATGATTCTCTGAGC',
                                  sseq='CGAACGCAATGATTCTCTGAGC')])
    amps = find_valid_amplicons(fwd, rev, max_amplicon_size=500)
    assert len(amps) == 1
    assert amps[0]['amplicon_start'] == 1720  # min(2000, 1720)
    assert amps[0]['amplicon_end'] == 2000    # max(2000, 1720)
    assert amps[0]['amplicon_size'] == 281    # 2000 - 1720 + 1


def test_find_amplicons_same_strand_rejected():
    """Both primers on the same strand cannot form an amplicon."""
    fwd = pd.DataFrame([make_hit(sseqid='c1', sstart=100, send=116, mismatch=0)])
    rev = pd.DataFrame([make_hit(qseqid='TestAssay_rev', sseqid='c1',
                                  sstart=300, send=316, mismatch=0)])  # both + strand
    amps = find_valid_amplicons(fwd, rev, max_amplicon_size=500)
    assert len(amps) == 0


# --- check_probe_in_amplicons ---

def _make_amplicon(contig='c1', start=100, end=300, fwd_mm=0, rev_mm=0):
    return {
        'contig_id': contig, 'amplicon_start': start, 'amplicon_end': end,
        'amplicon_size': end - start + 1, 'fwd_mismatches': fwd_mm,
        'rev_mismatches': rev_mm, 'probe_found': False,
        'probe_mismatches': None, 'probe_strand': None, 'amplicon_sequence': '',
    }


_PROBE_SEQ = 'ACGGGACAAAAAGGATGGCGAGTAC'  # 25 nt


def test_probe_within_amplicon():
    amps = [_make_amplicon()]
    probe_hits = pd.DataFrame([make_hit(qseqid='TestAssay_probe', sseqid='c1',
                                         qend=25, sstart=150, send=174, mismatch=0,
                                         qseq=_PROBE_SEQ, sseq=_PROBE_SEQ)])
    result = check_probe_in_amplicons(amps, probe_hits, probe_seq=_PROBE_SEQ, max_probe_mismatches=1)
    assert result[0]['probe_found'] is True
    assert result[0]['probe_mismatches'] == 0
    assert result[0]['probe_strand'] == '+'


def test_probe_outside_amplicon():
    amps = [_make_amplicon()]
    probe_hits = pd.DataFrame([make_hit(qseqid='TestAssay_probe', sseqid='c1',
                                         qend=25, sstart=500, send=524, mismatch=0,
                                         qseq=_PROBE_SEQ, sseq=_PROBE_SEQ)])
    result = check_probe_in_amplicons(amps, probe_hits, probe_seq=_PROBE_SEQ, max_probe_mismatches=1)
    assert result[0]['probe_found'] is False


def test_probe_too_many_mismatches():
    amps = [_make_amplicon()]
    # 2 true mismatches in sseq (positions 23,24 changed) so IUPAC count == 2
    sseq_mm2 = _PROBE_SEQ[:-2] + 'TT'
    probe_hits = pd.DataFrame([make_hit(qseqid='TestAssay_probe', sseqid='c1',
                                         qend=25, sstart=150, send=174, mismatch=2,
                                         qseq=_PROBE_SEQ, sseq=sseq_mm2)])
    result = check_probe_in_amplicons(amps, probe_hits, probe_seq=_PROBE_SEQ, max_probe_mismatches=1)
    assert result[0]['probe_found'] is False


def test_probe_minus_strand_detected():
    amps = [_make_amplicon()]
    probe_hits = pd.DataFrame([make_hit(qseqid='TestAssay_probe', sseqid='c1',
                                         qend=25, sstart=174, send=150, mismatch=0,
                                         qseq=_PROBE_SEQ, sseq=_PROBE_SEQ)])
    result = check_probe_in_amplicons(amps, probe_hits, probe_seq=_PROBE_SEQ, max_probe_mismatches=1)
    assert result[0]['probe_found'] is True
    assert result[0]['probe_strand'] == '-'


# --- call_detection ---

def test_call_detected():
    amps = [{'probe_found': True, 'amplicon_size': 142, 'contig_id': 'c1'}]
    call, n, multi, sizes, contigs = call_detection(amps)
    assert call == 'Detected'
    assert n == 1
    assert multi is False
    assert '142' in sizes


def test_call_primer_only():
    amps = [{'probe_found': False, 'amplicon_size': 142, 'contig_id': 'c1'}]
    call, n, multi, sizes, contigs = call_detection(amps)
    assert call == 'Primer Only'
    assert n == 1


def test_call_not_detected():
    call, n, multi, sizes, contigs = call_detection([])
    assert call == 'Not Detected'
    assert n == 0
    assert sizes == ''


def test_call_multi_amplicon():
    amps = [
        {'probe_found': True, 'amplicon_size': 142, 'contig_id': 'c1'},
        {'probe_found': True, 'amplicon_size': 198, 'contig_id': 'c1'},
    ]
    call, n, multi, sizes, contigs = call_detection(amps)
    assert call == 'Detected'
    assert n == 2
    assert multi is True
    assert '142' in sizes
    assert '198' in sizes


# --- load_blast_results ---

def test_load_blast_empty(tmp_path):
    f = tmp_path / "empty.tsv"
    f.write_text('')
    df = load_blast_results(str(f))
    assert df.empty
    assert list(df.columns) == BLAST_COLS


def test_load_blast_nonexistent_file(tmp_path):
    df = load_blast_results(str(tmp_path / "no_such_file.tsv"))
    assert df.empty
    assert list(df.columns) == BLAST_COLS


# --- run_ispcr integration ---

def test_run_ispcr_integration(tmp_path):
    """End-to-end test: one assay, one assembly, full detection."""
    from run_ispcr import run_ispcr

    # Minimal assay table — VhPath only, no IDT modifications in these seqs
    assay_csv = tmp_path / "assay_table.csv"
    assay_csv.write_text(
        "assay,probe,fwd,rev\n"
        "VhPath,ACGGGACAAAAAGGATGGCGAGTAC,AGCCGAGCGTTACCAGC,CGAACGCAATGATTCTCTGAGC\n"
    )

    # Fake .fna — single contig, enough sequence to contain the amplicon
    fna_path = tmp_path / "GCF_000001.fna"
    # 400 nt sequence; amplicon will be at positions 100-300
    seq = "A" * 400
    fna_path.write_text(f">contig1\n{seq}\n")

    # BLAST TSV: fwd hit (+strand), rev hit (-strand), probe hit within amplicon
    # fwd: VhPath_fwd, contig1, 100% id, len=17, 0 mm, sstart=100, send=116
    # rev: VhPath_rev, contig1, 100% id, len=22, 0 mm, sstart=300, send=279  (rev strand: sstart>send)
    # probe: VhPath_probe, contig1, within amplicon, sstart=150, send=174
    fwd_seq = "AGCCGAGCGTTACCAGC"
    rev_seq = "CGAACGCAATGATTCTCTGAGC"
    probe_seq = "ACGGGACAAAAAGGATGGCGAGTAC"
    blast_tsv = tmp_path / "blast.tsv"
    blast_tsv.write_text(
        "\t".join(["VhPath_fwd", "contig1", "100.0", "17", "0", "0",
                   "1", "17", "100", "116", "0.001", "32.0", fwd_seq, fwd_seq]) + "\n" +
        "\t".join(["VhPath_rev", "contig1", "100.0", "22", "0", "0",
                   "1", "22", "300", "279", "0.001", "44.0", rev_seq, rev_seq]) + "\n" +
        "\t".join(["VhPath_probe", "contig1", "100.0", "25", "0", "0",
                   "1", "25", "150", "174", "0.001", "50.0", probe_seq, probe_seq]) + "\n"
    )

    det_df = run_ispcr(
        blast_tsv=str(blast_tsv),
        assay_table=str(assay_csv),
        fna_path=str(fna_path),
        max_primer_mismatches=2,
        prime3_exact_nt=3,
        max_probe_mismatches=1,
        max_amplicon_size=500,
        store_amplicon_sequences=False,
    )

    assert len(det_df) == 1
    row = det_df.iloc[0]
    assert row['assay'] == 'VhPath'
    assert row['detection_call'] == 'Detected'
    assert row['n_amplicons'] == 1
    assert not row['multi_amplicon_flag']
    assert '201' in str(row['amplicon_sizes'])  # 300 - 100 + 1
    assert '100' in str(row['amplicon_starts']) and '300' in str(row['amplicon_ends'])


def test_run_ispcr_not_detected(tmp_path):
    """No BLAST hits → Not Detected for all assays."""
    from run_ispcr import run_ispcr

    assay_csv = tmp_path / "assay_table.csv"
    assay_csv.write_text(
        "assay,probe,fwd,rev\n"
        "VhPath,ACGGGACAAAAAGGATGGCGAGTAC,AGCCGAGCGTTACCAGC,CGAACGCAATGATTCTCTGAGC\n"
    )
    fna_path = tmp_path / "GCF_000001.fna"
    fna_path.write_text(">contig1\n" + "A" * 400 + "\n")
    blast_tsv = tmp_path / "blast.tsv"
    blast_tsv.write_text("")  # empty file

    det_df = run_ispcr(
        blast_tsv=str(blast_tsv),
        assay_table=str(assay_csv),
        fna_path=str(fna_path),
        max_primer_mismatches=2,
        prime3_exact_nt=3,
        max_probe_mismatches=1,
        max_amplicon_size=500,
        store_amplicon_sequences=False,
    )

    assert len(det_df) == 1
    assert det_df.iloc[0]['detection_call'] == 'Not Detected'
    # merged single frame carries the joined amplicon columns (empty when no amplicons)
    for col in ('amplicon_starts', 'amplicon_ends', 'amplicon_sequences'):
        assert col in det_df.columns


# --- call_detection mixed-probe multi-amplicon ---

def test_call_mixed_probe_multi_amplicon():
    """Two amplicons, only one probe-positive → Detected, multi_flag=True."""
    amps = [
        {'probe_found': True, 'amplicon_size': 142, 'contig_id': 'c1'},
        {'probe_found': False, 'amplicon_size': 198, 'contig_id': 'c1'},
    ]
    call, n, multi, sizes, contigs = call_detection(amps)
    assert call == 'Detected'
    assert n == 2
    assert multi is True


# --- IUPAC primitives ---

def test_iupac_match_exact_bases():
    """Standard bases match themselves."""
    for b in 'ACGT':
        assert iupac_match(b, b) is True

def test_iupac_match_degenerate_valid():
    """Degenerate codes match all bases in their set."""
    assert iupac_match('R', 'A') is True
    assert iupac_match('R', 'G') is True
    assert iupac_match('Y', 'C') is True
    assert iupac_match('Y', 'T') is True
    assert iupac_match('N', 'A') is True
    assert iupac_match('N', 'T') is True

def test_iupac_match_degenerate_invalid():
    """Degenerate codes don't match bases outside their set."""
    assert iupac_match('R', 'C') is False
    assert iupac_match('R', 'T') is False
    assert iupac_match('Y', 'A') is False
    assert iupac_match('Y', 'G') is False

def test_iupac_match_gap_never_matches():
    """Gap character never matches any base."""
    assert iupac_match('A', '-') is False
    assert iupac_match('N', '-') is False

def test_count_iupac_mismatches_perfect():
    assert count_iupac_mismatches('ACGT', 'ACGT') == 0

def test_count_iupac_mismatches_degenerate_no_mismatch():
    """R opposite A or G = 0 mismatches."""
    assert count_iupac_mismatches('ACRG', 'ACAG') == 0
    assert count_iupac_mismatches('ACRG', 'ACGG') == 0

def test_count_iupac_mismatches_degenerate_mismatch():
    """R opposite C or T = 1 mismatch."""
    assert count_iupac_mismatches('ACRG', 'ACCG') == 1

def test_count_iupac_mismatches_gap_counts():
    """Gap in sseq counts as a mismatch."""
    assert count_iupac_mismatches('ACGT', 'AC-T') == 1

def test_check_3prime_exact_iupac_pass():
    """Degenerate base at 3' position matches compatible subject base."""
    assert check_3prime_exact('ACGR', 'ACGA') is True
    assert check_3prime_exact('ACGR', 'ACGG') is True

def test_check_3prime_exact_iupac_fail():
    """Degenerate base at 3' position fails on incompatible subject base."""
    assert check_3prime_exact('ACGR', 'ACGC') is False
    assert check_3prime_exact('ACGR', 'ACGT') is False

# --- filter_primer_hits with IUPAC and gapopen ---

def test_filter_primer_hits_iupac_match_accepted():
    """Degenerate primer base matching subject within threshold is accepted."""
    # Primer: AGCCGAGCGTTACCAGR (17 nt, R at end)
    # Subject: AGCCGAGCGTTACCAGA — R matches A → 0 IUPAC mismatches
    hits = pd.DataFrame([make_hit(
        qend=17, mismatch=1,  # BLAST overcounts: reports 1 mismatch for R vs A
        qseq='AGCCGAGCGTTACCAGR',
        sseq='AGCCGAGCGTTACCAGA',
        sstart=100, send=116,
    )])
    result = filter_primer_hits(hits, primer_seq='AGCCGAGCGTTACCAGR',
                                 max_mismatch=0, prime3_exact=3)
    assert len(result) == 1
    assert result.iloc[0]['mismatch'] == 0  # IUPAC-corrected count


def test_filter_primer_hits_iupac_true_mismatch_rejected():
    """Degenerate primer base NOT matching subject is still a mismatch."""
    # Primer: AGCCGAGCGTTACCAGR (R at end)
    # Subject: AGCCGAGCGTTACCAGC — R does NOT match C → 1 IUPAC mismatch
    hits = pd.DataFrame([make_hit(
        qend=17, mismatch=1,
        qseq='AGCCGAGCGTTACCAGR',
        sseq='AGCCGAGCGTTACCAGC',
        sstart=100, send=116,
    )])
    result = filter_primer_hits(hits, primer_seq='AGCCGAGCGTTACCAGR',
                                 max_mismatch=0, prime3_exact=3)
    assert len(result) == 0


def test_filter_primer_hits_gapopen_rejected():
    """Hits with gapopen > 0 are rejected regardless of mismatch count."""
    hits = pd.DataFrame([{**make_hit(qend=17, mismatch=0,
                                      qseq='AGCCGAGCGTTACCAGC',
                                      sseq='AGCCGAGCGTTACCAGC',
                                      sstart=100, send=116),
                           'gapopen': 1}])
    result = filter_primer_hits(hits, primer_seq='AGCCGAGCGTTACCAGC',
                                 max_mismatch=2, prime3_exact=3)
    assert len(result) == 0


# --- call_detection: probe-free assays ---

def test_call_detection_no_probe_with_amplicon():
    """has_probe=False: any valid amplicon → Detected (not Primer Only)."""
    amps = [{'contig_id': 'c1', 'amplicon_size': 200, 'probe_found': False,
             'amplicon_start': 100, 'amplicon_end': 299}]
    call, n, multi, sizes, contigs = call_detection(amps, has_probe=False)
    assert call == 'Detected'
    assert n == 1
    assert not multi


def test_call_detection_no_probe_no_amplicon():
    """has_probe=False + no amplicons → Not Detected."""
    call, n, multi, sizes, contigs = call_detection([], has_probe=False)
    assert call == 'Not Detected'
    assert n == 0


def test_call_detection_no_probe_multi_amplicon():
    """has_probe=False + 2 amplicons → Detected, multi_flag True."""
    amps = [
        {'contig_id': 'c1', 'amplicon_size': 200, 'probe_found': False,
         'amplicon_start': 100, 'amplicon_end': 299},
        {'contig_id': 'c1', 'amplicon_size': 210, 'probe_found': False,
         'amplicon_start': 500, 'amplicon_end': 709},
    ]
    call, n, multi, sizes, contigs = call_detection(amps, has_probe=False)
    assert call == 'Detected'
    assert n == 2
    assert multi


def test_check_probe_empty_seq_passthrough():
    """Empty probe_seq returns amplicons unchanged (probe_found stays False)."""
    amps = [{'contig_id': 'c1', 'amplicon_start': 100, 'amplicon_end': 300,
             'amplicon_size': 201, 'fwd_mismatches': 0, 'rev_mismatches': 0,
             'probe_found': False, 'probe_mismatches': None, 'probe_strand': None,
             'amplicon_sequence': ''}]
    probe_hits = pd.DataFrame(columns=['qseqid', 'sseqid', 'pident', 'length',
                                        'mismatch', 'gapopen', 'qstart', 'qend',
                                        'sstart', 'send', 'evalue', 'bitscore',
                                        'qseq', 'sseq'])
    result = check_probe_in_amplicons(amps, probe_hits, probe_seq='', max_probe_mismatches=1)
    assert result[0]['probe_found'] is False
    assert result[0]['probe_mismatches'] is None
