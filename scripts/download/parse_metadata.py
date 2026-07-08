#!/usr/bin/env python3
"""
parse_metadata.py

Build a metadata CSV from the NCBI `datasets summary` JSON-Lines output
produced by download_assemblies.sh. Called automatically by that script;
can also be run standalone.

Usage
-----
    python3 parse_metadata.py \
        --jsonl  assemblies/.tmp/summary_raw.jsonl \
        --out    assemblies/metadata.csv
"""

import argparse
import csv
import json
import sys
from pathlib import Path


def g(d, *keys, default=""):
    """Safely walk a chain of nested dict keys, returning `default` if missing."""
    cur = d
    for k in keys:
        if not isinstance(cur, dict):
            return default
        cur = cur.get(k)
        if cur is None:
            return default
    return cur


def biosample_attrs(rec):
    """Return {attribute_name: value} from the BioSample attributes list."""
    attrs = g(rec, "assembly_info", "biosample", "attributes", default=[])
    out = {}
    if isinstance(attrs, list):
        for a in attrs:
            name = a.get("name")
            if name:
                out[name] = a.get("value", "")
    return out


# BioSample attribute names we surface as dedicated columns (harmonized names).
# Multiple aliases are checked in order; first non-empty wins.
BIOSAMPLE_FIELDS = {
    "bs_strain":           ["strain"],
    "bs_isolate":          ["isolate"],
    "bs_host":             ["host", "specific_host"],
    "bs_isolation_source": ["isolation_source", "isolation source"],
    "bs_collection_date":  ["collection_date", "collection date"],
    "bs_geo_loc_name":     ["geo_loc_name", "geographic location (country and/or sea)"],
    "bs_lat_lon":          ["lat_lon"],
    "bs_sample_type":      ["sample_type"],
    "bs_env_broad_scale":  ["env_broad_scale"],
    "bs_env_local_scale":  ["env_local_scale"],
    "bs_env_medium":       ["env_medium"],
    "bs_serovar":          ["serovar"],
    "bs_source_type":      ["source_type"],
}


def build_row(rec):
    attrs = biosample_attrs(rec)

    def attr(aliases):
        for a in aliases:
            if attrs.get(a):
                return attrs[a]
        return ""

    row = {
        # --- identifiers -----------------------------------------------------
        "accession":            rec.get("accession", ""),
        "paired_accession":     rec.get("paired_accession", ""),
        "source_database":      rec.get("source_database", ""),
        # --- taxonomy / organism --------------------------------------------
        "organism_name":        g(rec, "organism", "organism_name"),
        "tax_id":               g(rec, "organism", "tax_id"),
        "infraspecific_strain": g(rec, "organism", "infraspecific_names", "strain"),
        # --- assembly --------------------------------------------------------
        "assembly_name":        g(rec, "assembly_info", "assembly_name"),
        "assembly_level":       g(rec, "assembly_info", "assembly_level"),
        "assembly_status":      g(rec, "assembly_info", "assembly_status"),
        "assembly_type":        g(rec, "assembly_info", "assembly_type"),
        "assembly_method":      g(rec, "assembly_info", "assembly_method"),
        "release_date":         g(rec, "annotation_info", "release_date"),
        "submitter":            g(rec, "assembly_info", "biosample", "owner", "name"),
        "bioproject_accession": g(rec, "assembly_info", "bioproject_accession"),
        "biosample_accession":  g(rec, "assembly_info", "biosample", "accession"),
        "wgs_project":          g(rec, "wgs_info", "wgs_project_accession"),
        # --- annotation ------------------------------------------------------
        "annotation_pipeline":  g(rec, "annotation_info", "pipeline"),
        "gene_total":           g(rec, "annotation_info", "stats", "gene_counts", "total"),
        "gene_protein_coding":  g(rec, "annotation_info", "stats", "gene_counts", "protein_coding"),
        "gene_pseudogene":      g(rec, "annotation_info", "stats", "gene_counts", "pseudogene"),
        # --- assembly stats --------------------------------------------------
        "total_sequence_length": g(rec, "assembly_stats", "total_sequence_length"),
        "gc_percent":            g(rec, "assembly_stats", "gc_percent"),
        "number_of_contigs":     g(rec, "assembly_stats", "number_of_contigs"),
        "contig_n50":            g(rec, "assembly_stats", "contig_n50"),
        "number_of_scaffolds":   g(rec, "assembly_stats", "number_of_scaffolds"),
        "scaffold_n50":          g(rec, "assembly_stats", "scaffold_n50"),
        "total_chromosomes":     g(rec, "assembly_stats", "total_number_of_chromosomes"),
        "genome_coverage":       g(rec, "assembly_stats", "genome_coverage"),
        # --- quality: CheckM -------------------------------------------------
        "checkm_completeness":   g(rec, "checkm_info", "completeness"),
        "checkm_contamination":  g(rec, "checkm_info", "contamination"),
        "checkm_marker_set":     g(rec, "checkm_info", "checkm_marker_set"),
        # --- quality: ANI ----------------------------------------------------
        "ani_match_status":      g(rec, "average_nucleotide_identity", "match_status"),
        "ani_taxonomy_check":    g(rec, "average_nucleotide_identity", "taxonomy_check_status"),
        "ani_best_match_organism": g(rec, "average_nucleotide_identity", "best_ani_match", "organism_name"),
        "ani_best_match_value":  g(rec, "average_nucleotide_identity", "best_ani_match", "ani"),
    }

    # BioSample harmonized attributes (curated subset)
    for col, aliases in BIOSAMPLE_FIELDS.items():
        row[col] = attr(aliases)

    return row


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--jsonl", required=True,
                    help="datasets summary JSON-Lines file (produced by download_assemblies.sh)")
    ap.add_argument("--out", required=True,
                    help="output metadata CSV path")
    args = ap.parse_args()

    jsonl_path = Path(args.jsonl)
    if not jsonl_path.is_file():
        sys.exit(f"ERROR: JSONL not found: {jsonl_path}")

    rows = []
    bad = 0
    with open(jsonl_path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            try:
                rows.append(build_row(json.loads(line)))
            except json.JSONDecodeError:
                bad += 1

    if not rows:
        sys.exit("ERROR: no records parsed -- is the JSONL empty?")

    # Stable column order: identifiers/taxonomy/assembly first, then the rest.
    fieldnames = list(rows[0].keys())
    extra = sorted({k for r in rows for k in r} - set(fieldnames))
    fieldnames += extra

    out_path = Path(args.out)
    with open(out_path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)

    n_with_org = sum(1 for r in rows if r["organism_name"])
    print(f"Wrote {len(rows)} rows x {len(fieldnames)} columns -> {out_path}")
    print(f"  rows with organism name: {n_with_org}/{len(rows)}")
    if bad:
        print(f"  WARNING: {bad} malformed JSON line(s) skipped")


if __name__ == "__main__":
    main()
