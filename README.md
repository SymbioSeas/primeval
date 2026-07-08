# primeval: in silico PCR assay validation

**primeval** evaluates the sensitivity and specificity of PCR assays (including probe-based assays designed for ddPCR/qPCR) against a user-provided set of genome assemblies. For each assay, it identifies valid amplicons using BLAST-based primer alignment and reports detection calls, mismatch counts, and amplicon sizes per assembly.

## Two tools in this repository

| Tool                         | Location                                    | Purpose                                                                                                                                                                                                                                               |
| ---------------------------- | ------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **primeval** (pipeline)      | repository root (`workflow/`, `config/`, …) | Validate PCR assays in silico against a set of genome assemblies.                                                                                                                                                                                     |
| **assay-design** (companion) | [`assay-design/`](assay-design/)            | Identify clade-conserved / clade-specific orthologs from a Panaroo pangenome and extract representative sequences. This is the candidate-gene discovery step used to develop the Vpop assays. See [`assay-design/README.md`](assay-design/README.md). |

Both tools share the single conda environment defined in `environment.yaml`.

## Features

- Evaluates primer and probe binding across a user-provided set of genome assemblies
- Supports probe-based assays (hydrolysis probes, e.g., ddPCR/qPCR) and probe-free assays (SYBR/dsDNA)
- IUPAC degenerate base support in all primer and probe sequences
- Handles primer binding on either strand of an assembly
- Configurable mismatch tolerances and 3′-exact match requirements
- Outputs per-assay detection summaries, amplicon details, and detection heatmaps
- Runs locally or on SLURM HPC clusters via Snakemake profiles

## Requirements

- [Conda](https://docs.conda.io/en/latest/) or [Mamba](https://mamba.readthedocs.io/)
- [NCBI datasets CLI](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/) (for downloading assemblies)
- To install NCBI datasets CLI with conda:
```
conda install -c conda-forge ncbi-datasets-cli
```

## Installation

```bash
git clone https://github.com/SymbioSeas/primeval.git
cd primeval
conda env create -f environment.yaml
conda activate primeval
```

The pipeline runs inside this activated environment; the Snakemake profiles set
`use-conda: false` so no per-rule environments are built.

## Quick start

### 1. Download assemblies

Use the provided helper script to download RefSeq assemblies for your taxon of interest:

```bash
bash scripts/download/download_assemblies.sh -t "Vibrionaceae" -o assemblies/
```

This downloads all RefSeq assemblies (complete through contig level) for the specified taxon and writes a `metadata.csv` to the output directory. See [Downloading assemblies](#downloading-assemblies) for options and HPC usage.

### 2. Configure the pipeline

Edit `config/config.yaml` to point at your assemblies and set detection thresholds:

```yaml
assembly_dir: "assemblies"          # directory containing .fna files
metadata: "assemblies/metadata.csv" # metadata CSV from download step
results_dir: "results"              # all outputs written here
assay_table: "assay_table.csv"      # your assay definitions

max_primer_mismatches: 2            # mismatches allowed per primer
prime3_exact_nt: 3                  # 3′-terminal bases that must match exactly
max_probe_mismatches: 1             # mismatches allowed in probe
max_amplicon_size: 500              # maximum amplicon size (bp)
```

### 3. Prepare your assay table

Create a CSV file named `assay_table.csv` with one row per assay (see [Assay table format](#assay-table-format)).

### 4. Run

```bash
# Local (8 cores)
snakemake --profile workflow/profiles/local

# SLURM cluster
snakemake --profile workflow/profiles/slurm
```

---

## Assay table format

The assay table is a CSV file with the following columns:

| Column  | Required | Description                                                                               |
| ------- | -------- | ----------------------------------------------------------------------------------------- |
| `assay` | Yes      | Unique assay name (used in all output files)                                              |
| `fwd`   | Yes      | Forward primer sequence (5′→3′)                                                           |
| `rev`   | Yes      | Reverse primer sequence (5′→3′, same orientation as fwd - primeval handles RC internally) |
| `probe` | No       | Probe sequence (5′→3′). Leave empty for probe-free (SYBR) assays                          |

**Sequence notation:**
- Standard IUPAC ambiguity codes are supported (R, Y, S, W, K, M, B, D, H, V, N)
- Modifications can be noted inline using `/ModName/` or `[ModName]` notation; these are stripped before alignment (e.g., `/56-FAM/ACGT[BHQ1]` → `ACGT`)

Example:

```
assay,probe,fwd,rev
Assay1,ACGGGACAAAAAGGATGGCGAGTAC,AGCCGAGCGTTACCAGC,CGAACGCAATGATTCTCTGAGC
Assay2,,GCTACGCCCTCCATCATCC,GCGCGTGATTATCTGATAGC
```

---

## Detection thresholds

primeval reports a detection call per assay per assembly using thresholds set in
`config/config.yaml`. The defaults reflect PCR biochemistry:

| Parameter               | Default | Rationale                                                                                                                                                                                                                         |
| ----------------------- | ------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `max_primer_mismatches` | 2       | A primer with one or two **internal** mismatches still typically primes efficiently; this tolerates strain-level SNPs while excluding poor binders. Counted IUPAC-aware (a degenerate base matches any of its represented bases). |
| `prime3_exact_nt`       | 3       | Mismatches at the 3′ terminus inhibit polymerase extension, so the last several bases must match exactly regardless of `max_primer_mismatches`.                                                                                   |
| `max_probe_mismatches`  | 1       | Used for probe-based assays only. Hydrolysis probes tolerate less mismatch than primers, so the default is stricter.                                                                                                              |
| `max_amplicon_size`     | 500     | Typical qPCR/ddPCR amplicons are ~70–200 bp; 500 bp captures valid products while rejecting spurious long-range primer pairings.                                                                                                  |

**Detection calls:**
- `Detected` — valid amplicon found with probe contained within it (probe assays), or valid amplicon found (probe-free assays)
- `Primer Only` — valid amplicon found but probe not detected within it
- `Not Detected` — no valid amplicon found

The BLAST search itself is run with deliberately permissive settings
(`evalue=1000`, `perc_identity=70`, `word_size=7`, tuned for short oligo
queries) so that no candidate binding site is missed; stringency is enforced
downstream by the mismatch, 3′-exact, and amplicon-size filters above.

---

## Outputs

All outputs are written to `results_dir` (default: `results/`):

```
results/
├── blast/                        # raw BLAST output per assembly (intermediate)
├── amplicons/
│   ├── {accession}.csv           # per-assembly detection calls (one row per assay)
│   └── {accession}_amplicons.csv # per-assembly amplicon details
└── reports/
    ├── species_summary.csv       # detection rates by species and assay
    ├── assay_summary.csv         # detection rates by assay across species groups
    ├── assay_summary.xlsx        # same, Excel format
    ├── detection_matrix.xlsx     # detection call matrix (assemblies × assays)
    ├── run_manifest.txt          # parameter log and tool versions
    └── figures/
        ├── heatmap_binary.pdf    # binary detection heatmap
        └── heatmap_verbose.pdf   # heatmap with mismatch details
```

---

## Downloading assemblies

### Basic usage

```bash
bash scripts/download/download_assemblies.sh -t "Taxon name" -o assemblies/
```

Options:

| Flag | Description | Default |
|------|-------------|---------|
| `-t TAXON` | Taxon name or NCBI tax ID (required; **repeatable** — see below) | — |
| `-o OUTDIR` | Output directory | `assemblies` |
| `-l LEVELS` | Assembly levels (comma-separated) | `complete,chromosome,scaffold,contig` |
| `-e EMAIL` | NCBI e-mail (or set `NCBI_EMAIL` env var) | — |
| `-k API_KEY` | NCBI API key for higher rate limits (or set `NCBI_API_KEY` env var) | — |

### Multiple taxa

Pass `-t` more than once to download the **de-duplicated union** of several taxa
into a single output directory:

```bash
bash scripts/download/download_assemblies.sh \
    -t "Vibrio jasicida" \
    -t "Vibrio owensii" \
    -o assemblies/
```

Each taxon is queried in turn; assemblies shared between taxa (e.g. when one
query nests inside another) are downloaded and written to `metadata.csv` only
once. Mixing names and tax IDs is fine, e.g. `-t "Vibrio owensii" -t 661487`.

### Setting your NCBI API key once

An [NCBI API key](https://www.ncbi.nlm.nih.gov/account/) raises your download
rate limit from 3 to 10 requests/sec. This is worth setting if you're downloading large sets of assemblies (i.e., hundreds or thousands of assemblies). Rather than passing `-k` every time, save it once:

```bash
cp config/ncbi_credentials.example.sh config/ncbi_credentials.sh
# edit config/ncbi_credentials.sh and paste your key into NCBI_API_KEY
```

`download_assemblies.sh` sources this file automatically on every run.
To keep the key elsewhere, point `PRIMEVAL_CREDENTIALS` at your own file.

The key is resolved as: **`-k` flag → `NCBI_API_KEY` environment variable →
credentials file** (first one set wins).

### HPC / SLURM

Local available storage requirements for primeval are directly scaled by the assembly dataset provided (i.e., you need space to store the downloaded assemblies you provide primeval!). If needed, primeval runs can easily be submitted in a SLURM environment using the wrapper below.

Wrap the script in an sbatch job for large downloads:

```bash
sbatch --time=24:00:00 --mem=8G \
  --wrap="bash scripts/download/download_assemblies.sh -t Vibrionaceae -o assemblies/"
```

The script is resume-aware: if interrupted, re-running it will skip assemblies already successfully downloaded.

---

## Assay specificity validation

primeval tests assay sensitivity and specificity against **your input assembly dataset**. The scope of specificity testing is therefore determined by which assemblies you provide.

**Recommended workflow for specificity screening:**

1. **Single-primer BLAST screen** (NCBI web interface): Individually BLAST each primer and probe sequence against the NCBI `nt` database, *excluding* your target taxon. This identifies any off-target binding sites outside your group of interest. If no hits are returned for any oligo, off-target amplification outside the taxon is extremely unlikely (a primer must bind for any amplicon to form).

2. **Expand the input dataset if needed**: If step 1 returns hits in a non-target taxon, download assemblies from that taxon and add them to your `assembly_dir`. primeval will then determine whether those single-primer hits form complete, detectable amplicons.

This two-stage approach is computationally efficient, such that you only download and evaluate assemblies in taxa where off-target primer binding is possible.

> **Planned feature:** A future release will support BLASTing directly against NCBI pre-built reference databases (e.g., `ref_prok_rep_genomes`) as a single-step broader specificity check, without requiring manual assembly downloads.

---

## Test dataset

A small set of 5 assemblies is provided for validating your installation:

```bash
bash test_data/download_test_data.sh
```

Then update `config/config.yaml`:

```yaml
assembly_dir: "test_data/assemblies"
metadata: "test_data/assemblies/metadata.csv"
```

And run:

```bash
snakemake --profile workflow/profiles/local
```

The test dataset covers all detection scenarios: `Detected` (including via minus-strand primer binding), `Primer Only`, and `Not Detected`.

---

## Citation

If you use primeval in your research, please cite:

> Smith S, et al. (2026) *[manuscript title]*. *[journal]*. doi:[doi]

---

## License

MIT — see [LICENSE](LICENSE).
