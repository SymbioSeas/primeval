#!/usr/bin/env bash
# Run the assay-design pipeline on the bundled worked example.
# From the example/ directory:  bash run_example.sh
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ASSAY_DESIGN_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

python3 "${ASSAY_DESIGN_DIR}/run_pipeline.py" \
    --matrix          "${SCRIPT_DIR}/example_gene_presence_absence.csv" \
    --isolates-dir    "${SCRIPT_DIR}/isolate_groups" \
    --gene-data       "${SCRIPT_DIR}/gene_data.csv" \
    --representatives "${SCRIPT_DIR}/representatives.tsv" \
    --output-dir      "${SCRIPT_DIR}/output"

echo ""
echo "Outputs written to ${SCRIPT_DIR}/output/"
