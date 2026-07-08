# assay-design — clade-specific ortholog discovery

Companion tool to [primeval](../README.md). `assay-design` parses a
[Panaroo](https://gtonkinhill.github.io/panaroo/) pangenome to find orthologs
that are **conserved within** a clade and **specific to** it (absent from the
rest of the dataset), then extracts representative protein and nucleotide
sequences. These clade-specific genes are the candidate targets from which the
Vpop dPCR assays were designed; the resulting amplicons are then validated with
primeval.

## Workflow

```
Panaroo output ──▶ Stage 1: identify conserved & specific orthologs ──▶ Stage 2: extract representative sequences
                   (parse_group_specific_orthologs.py)                   (extract_proteins.py)
```

Both stages are wrapped by `run_pipeline.py`.

## Environment

Uses the repository's single top-level environment — no separate install:

```bash
conda env create -f ../environment.yaml   # if not already created
conda activate primeval
```

Only `python`, `pandas`, and `numpy` are needed for this tool (all included in
that environment).

## Inputs

| Input | Format |
|-------|--------|
| **Presence/absence matrix** (`--matrix`) | Panaroo `gene_presence_absence` (`.csv` or `.tsv`). Metadata columns (`Gene`, `Annotation`, …) are detected automatically; all other columns are treated as genomes. Each cell holds the locus ID for that genome (empty = gene absent). |
| **Inclusion groups** (`--isolates-dir`) | A directory of `.txt` files, one per group. Each file lists the genome names (one per line) that make up that inclusion group. Names must match matrix column headers. The filename stem names the group in all outputs. |
| **Gene data** (`--gene-data`) | Panaroo `gene_data.csv` with columns `gff_file`, `annotation_id`, `prot_sequence`, `dna_sequence`, `gene_name`, `description`. Only needed for Stage 2. |
| **Representatives** (`--representatives`) | TSV with columns `group_stem` and `representative_assembly`. Chooses which genome's sequences are extracted per group. `group_stem` must match an inclusion-group filename stem; `representative_assembly` must be a matrix column. Only needed for Stage 2. |

For a given group, the **exclusion set** is every genome column in the matrix
that is not in that group.

## Usage

Run both stages on one dataset:

```bash
python run_pipeline.py \
    --matrix          gene_presence_absence.csv \
    --isolates-dir    isolate_groups \
    --gene-data       gene_data.csv \
    --representatives representatives.tsv \
    --output-dir      results
```

Useful flags:

- `--threshold 1.0` — strict all-or-nothing (present in *all* inclusion, absent from *all* exclusion)
- `--skip-stage2` — ortholog identification only (no `--gene-data`/`--representatives` needed)
- `--min-match-rate 0.5` — Stage 2 guard (see below)

Or run the stages individually:

```bash
python parse_group_specific_orthologs.py --matrix gene_presence_absence.csv --isolates-dir isolate_groups
python extract_proteins.py --orthologs-dir <dataset>_orthologs --gene-data gene_data.csv --representatives representatives.tsv
```

## Worked example

A tiny, ready-to-run example (subsampled from the FMS2026 dataset) lives in
[`example/`](example/):

```bash
cd example
bash run_example.sh
```

It runs in seconds and produces conserved/specific ortholog tables and
representative FASTAs for three *V. mediterranei* clades. See
[`example/README.md`](example/README.md).

## Outputs

Per inclusion group:

| File | Contents |
|------|----------|
| `<group>_conserved_orthologs.tsv` | Orthogroups conserved in the group, with inclusion/exclusion coverage metrics (`n_/frac_inclusion_present`, `n_/frac_exclusion_present`, `exclusion_assemblies_present`). |
| `<group>_specific_orthologs.tsv` | The conserved subset that is also specific to the group (clade-specific candidates). |
| `<group>_conserved_proteins.faa` | Representative protein sequences (Stage 2). |
| `<group>_conserved_genes.fna` | Representative nucleotide sequences (Stage 2). |

## Thresholds and rationale

**`--threshold` (default 0.9)** governs both calls:

- **Conserved** — present in **≥ 90%** of the group's inclusion genomes. Requiring
  near-universal (rather than strictly universal) presence tolerates draft-genome
  incompleteness and occasional annotation/clustering dropout, while still
  demanding the gene be a stable feature of the clade — a prerequisite for a
  sensitive assay.
- **Specific** — conserved **and** present in **≤ 10%** of the exclusion genomes
  (i.e. `1 − threshold`). Allowing a small exclusion fraction absorbs rare
  horizontal gene transfer and annotation noise without letting genuinely shared
  genes through — a prerequisite for a specific assay.

Set `--threshold 1.0` for the strict definition (present in every inclusion
genome, absent from every exclusion genome). Values below ~0.8 are not
recommended: they admit orthologs too patchy to make reliable assay targets.

**`--min-match-rate` (default 0.5)** is a Stage 2 safety guard. Sequences are
extracted by matching each conserved ortholog's representative locus ID against
`gene_data` `annotation_id`s. If a genome's annotation-ID format in `gene_data`
differs from the presence/absence matrix, fewer than this fraction of IDs will
match and the group is **skipped with a warning** rather than producing a
misleadingly truncated FASTA. Raise it to be stricter about extraction
completeness.

## Note on functional annotation

Downstream functional annotation of the extracted proteins (e.g. with
[EggNOG-mapper](http://eggnog-mapper.embl.de/)) is outside the scope of this
tool. The `.faa` files produced by Stage 2 are directly usable as input to such
tools if desired.
