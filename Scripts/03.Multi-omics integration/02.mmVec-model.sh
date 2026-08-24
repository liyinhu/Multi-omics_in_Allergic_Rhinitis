#!/usr/bin/env bash

# Usage:
#   bash 02.mmVec-model.sh [INPUT_DIR] [METADATA_DIR] [OUTPUT_DIR]
#   bash 02.mmVec-model.sh -i INPUT_DIR -m METADATA_DIR -o OUTPUT_DIR
#
# Defaults:
#   INPUT_DIR=current directory
#   METADATA_DIR=metadata
#   OUTPUT_DIR=current directory
#
# Notes:
#   Run this script in a qiime2 environment with `biom` and `qiime` available.
#   The mmVec/qiime2 command logic and model parameters are kept unchanged.

set -euo pipefail

usage() {
    sed -n '2,13p' "$0"
}

INPUT_DIR="."
METADATA_DIR="metadata"
OUTPUT_DIR="."

POSITIONAL_ARGS=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        -i|--input-dir)
            INPUT_DIR="$2"
            shift 2
            ;;
        -m|--metadata-dir)
            METADATA_DIR="$2"
            shift 2
            ;;
        -o|--output-dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            POSITIONAL_ARGS+=("$1")
            shift
            ;;
    esac
done

if [[ ${#POSITIONAL_ARGS[@]} -gt 3 ]]; then
    usage
    exit 1
fi

if [[ ${#POSITIONAL_ARGS[@]} -ge 1 ]]; then
    INPUT_DIR="${POSITIONAL_ARGS[0]}"
fi
if [[ ${#POSITIONAL_ARGS[@]} -ge 2 ]]; then
    METADATA_DIR="${POSITIONAL_ARGS[1]}"
fi
if [[ ${#POSITIONAL_ARGS[@]} -ge 3 ]]; then
    OUTPUT_DIR="${POSITIONAL_ARGS[2]}"
fi

if ! command -v biom >/dev/null 2>&1; then
    echo "Error: biom command not found. Please activate the required environment first." >&2
    exit 1
fi

if ! command -v qiime >/dev/null 2>&1; then
    echo "Error: qiime command not found. Please activate the required environment first." >&2
    exit 1
fi

required_inputs=(
    "Microbiome.sample_aligned.txt"
    "Proteome.sample_aligned.txt"
    "Metabolome.sample_aligned.txt"
    "Lipidome.sample_aligned.txt"
)

for file in "${required_inputs[@]}"; do
    if [[ ! -f "${INPUT_DIR}/${file}" ]]; then
        echo "Error: required input file not found: ${INPUT_DIR}/${file}" >&2
        exit 1
    fi
done

required_metadata=(
    "Proteome_metadata.txt"
    "Metabolome_metadata.txt"
    "Lipidome_metadata.txt"
    "Microbiome_metadata.txt"
)

for file in "${required_metadata[@]}"; do
    if [[ ! -f "${METADATA_DIR}/${file}" ]]; then
        echo "Error: required metadata file not found: ${METADATA_DIR}/${file}" >&2
        exit 1
    fi
done

mkdir -p "${OUTPUT_DIR}"

# Under qiime2-2020 environment

# Import microbiome, proteome, metabolome and lipidome
biom convert -i "${INPUT_DIR}/Microbiome.sample_aligned.txt" -o "${OUTPUT_DIR}/Microbiome.biom" --table-type="Taxon table" --to-json
biom convert -i "${INPUT_DIR}/Proteome.sample_aligned.txt" -o "${OUTPUT_DIR}/Proteome.biom" --table-type="Metabolite table" --to-json
biom convert -i "${INPUT_DIR}/Metabolome.sample_aligned.txt" -o "${OUTPUT_DIR}/Metabolome.biom" --table-type="Metabolite table" --to-json
biom convert -i "${INPUT_DIR}/Lipidome.sample_aligned.txt" -o "${OUTPUT_DIR}/Lipidome.biom" --table-type="Metabolite table" --to-json

# File transformation
qiime tools import \
        --input-path "${OUTPUT_DIR}/Microbiome.biom" \
        --input-format BIOMV100Format \
        --output-path "${OUTPUT_DIR}/Microbiome.qza" \
        --type FeatureTable[Frequency]

qiime tools import \
        --input-path "${OUTPUT_DIR}/Proteome.biom" \
        --input-format BIOMV100Format \
        --output-path "${OUTPUT_DIR}/Proteome.qza" \
        --type FeatureTable[Frequency]

qiime tools import \
        --input-path "${OUTPUT_DIR}/Metabolome.biom" \
        --input-format BIOMV100Format \
        --output-path "${OUTPUT_DIR}/Metabolome.qza" \
        --type FeatureTable[Frequency]

qiime tools import \
        --input-path "${OUTPUT_DIR}/Lipidome.biom" \
        --input-format BIOMV100Format \
        --output-path "${OUTPUT_DIR}/Lipidome.qza" \
        --type FeatureTable[Frequency]

## mmVec model training
# Microbiome-Proteome
qiime mmvec paired-omics \
        --i-microbes "${OUTPUT_DIR}/Microbiome.qza" \
        --i-metabolites "${OUTPUT_DIR}/Proteome.qza" \
        --p-epochs 100 \
        --p-latent-dim 3 \
        --p-learning-rate 0.001 \
        --p-summary-interval 1 \
        --p-no-arm-the-gpu \
        --output-dir "${OUTPUT_DIR}/Microbiome_Proteome_Summary"

# Microbiome-Metabolome
qiime mmvec paired-omics \
        --i-microbes "${OUTPUT_DIR}/Microbiome.qza" \
        --i-metabolites "${OUTPUT_DIR}/Metabolome.qza" \
        --p-epochs 100 \
        --p-latent-dim 3 \
        --p-learning-rate 0.001 \
        --p-summary-interval 1 \
        --p-no-arm-the-gpu \
        --output-dir "${OUTPUT_DIR}/Microbiome_Metabolome_Summary"

# Microbiome-Lipidome
qiime mmvec paired-omics \
        --i-microbes "${OUTPUT_DIR}/Microbiome.qza" \
        --i-metabolites "${OUTPUT_DIR}/Lipidome.qza" \
        --p-epochs 100 \
        --p-latent-dim 3 \
        --p-learning-rate 0.001 \
        --p-summary-interval 1 \
        --p-no-arm-the-gpu \
        --output-dir "${OUTPUT_DIR}/Microbiome_Lipidome_Summary"

## Export feature relations
# Microbiome-Proteome
qiime metadata tabulate \
        --m-input-file "${OUTPUT_DIR}/Microbiome_Proteome_Summary/conditionals.qza" \
        --o-visualization "${OUTPUT_DIR}/Microbiome_Proteome_conditionals_viz.qzv"

qiime emperor biplot \
        --i-biplot "${OUTPUT_DIR}/Microbiome_Proteome_Summary/conditional_biplot.qza" \
        --m-sample-metadata-file "${METADATA_DIR}/Proteome_metadata.txt" \
        --m-feature-metadata-file "${METADATA_DIR}/Microbiome_metadata.txt" \
        --p-number-of-features 15 \
        --p-ignore-missing-samples \
        --o-visualization "${OUTPUT_DIR}/Microbiome_Proteome_emperor.qzv"

qiime mmvec summarize-single \
        --i-model-stats "${OUTPUT_DIR}/Microbiome_Proteome_Summary/model_stats.qza" \
        --o-visualization "${OUTPUT_DIR}/Microbiome_Proteome_model_summary.qzv"

# Microbiome-Metabolome
qiime metadata tabulate \
        --m-input-file "${OUTPUT_DIR}/Microbiome_Metabolome_Summary/conditionals.qza" \
        --o-visualization "${OUTPUT_DIR}/Microbiome_Metabolome_conditionals_viz.qzv"

qiime emperor biplot \
        --i-biplot "${OUTPUT_DIR}/Microbiome_Metabolome_Summary/conditional_biplot.qza" \
        --m-sample-metadata-file "${METADATA_DIR}/Metabolome_metadata.txt" \
        --m-feature-metadata-file "${METADATA_DIR}/Microbiome_metadata.txt" \
        --p-number-of-features 15 \
        --p-ignore-missing-samples \
        --o-visualization "${OUTPUT_DIR}/Microbiome_Metabolome_emperor.qzv"

qiime mmvec summarize-single \
        --i-model-stats "${OUTPUT_DIR}/Microbiome_Metabolome_Summary/model_stats.qza" \
        --o-visualization "${OUTPUT_DIR}/Microbiome_Metabolome_model_summary.qzv"

# Microbiome-Lipidome
qiime metadata tabulate \
        --m-input-file "${OUTPUT_DIR}/Microbiome_Lipidome_Summary/conditionals.qza" \
        --o-visualization "${OUTPUT_DIR}/Microbiome_Lipidome_conditionals_viz.qzv"

qiime emperor biplot \
        --i-biplot "${OUTPUT_DIR}/Microbiome_Lipidome_Summary/conditional_biplot.qza" \
        --m-sample-metadata-file "${METADATA_DIR}/Lipidome_metadata.txt" \
        --m-feature-metadata-file "${METADATA_DIR}/Microbiome_metadata.txt" \
        --p-number-of-features 15 \
        --p-ignore-missing-samples \
        --o-visualization "${OUTPUT_DIR}/Microbiome_Lipidome_emperor.qzv"

qiime mmvec summarize-single \
        --i-model-stats "${OUTPUT_DIR}/Microbiome_Lipidome_Summary/model_stats.qza" \
        --o-visualization "${OUTPUT_DIR}/Microbiome_Lipidome_model_summary.qzv"

## Q2 calculation
# Microbiome-Proteome
qiime mmvec paired-omics \
        --i-microbes "${OUTPUT_DIR}/Microbiome.qza" \
        --i-metabolites "${OUTPUT_DIR}/Proteome.qza" \
        --p-latent-dim 0 \
        --p-summary-interval 1 \
        --output-dir "${OUTPUT_DIR}/Proteome_Summary"

qiime mmvec summarize-paired \
        --i-model-stats "${OUTPUT_DIR}/Microbiome_Proteome_Summary/model_stats.qza" \
        --i-baseline-stats "${OUTPUT_DIR}/Proteome_Summary/model_stats.qza" \
        --o-visualization "${OUTPUT_DIR}/Proteome_paired-summary.qzv"

# Microbiome-Metabolome
qiime mmvec paired-omics \
        --i-microbes "${OUTPUT_DIR}/Microbiome.qza" \
        --i-metabolites "${OUTPUT_DIR}/Metabolome.qza" \
        --p-latent-dim 0 \
        --p-summary-interval 1 \
        --output-dir "${OUTPUT_DIR}/Metabolome_Summary"

qiime mmvec summarize-paired \
        --i-model-stats "${OUTPUT_DIR}/Microbiome_Metabolome_Summary/model_stats.qza" \
        --i-baseline-stats "${OUTPUT_DIR}/Metabolome_Summary/model_stats.qza" \
        --o-visualization "${OUTPUT_DIR}/Metabolome_paired-summary.qzv"

# Microbiome-Lipidome
qiime mmvec paired-omics \
        --i-microbes "${OUTPUT_DIR}/Microbiome.qza" \
        --i-metabolites "${OUTPUT_DIR}/Lipidome.qza" \
        --p-latent-dim 0 \
        --p-summary-interval 1 \
        --output-dir "${OUTPUT_DIR}/Lipidome_Summary"

qiime mmvec summarize-paired \
        --i-model-stats "${OUTPUT_DIR}/Microbiome_Lipidome_Summary/model_stats.qza" \
        --i-baseline-stats "${OUTPUT_DIR}/Lipidome_Summary/model_stats.qza" \
        --o-visualization "${OUTPUT_DIR}/Lipidome_paired-summary.qzv"
