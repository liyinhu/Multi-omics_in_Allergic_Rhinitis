# Environment Notes

This repository has **not yet been fully validated in a single clean environment**. Individual scripts have been updated to improve command-line execution where applicable, but full end-to-end smoke tests with real or example data are still pending.

The project currently depends on a mix of:

- R scripts
- Python scripts
- external command-line tools for the QIIME2 / mmVec workflow

Exact package versions were **not fully pinned in the original repository**. Users should record their local `sessionInfo()` for R and the installed Python package versions used for each run. A formal `requirements.txt`, environment file, or R package list may be added later.

## R Dependencies by Module

### Proteome

- `OlinkAnalyze`
- `robustbase`
- `vegan`
- `ggplot2`
- `ggrepel`
- `clusterProfiler`
- `org.Hs.eg.db`
- `enrichplot`
- `dplyr`
- `stringr`
- `plotly`
- `grid`
- `forcats`

### Metabolome and Lipidome

- `robustbase`
- `vegan`
- `ggplot2`
- `ggrepel`
- `plotly`
- `grid`

### Multi-omics Integration

- `mixOmics`
- `MetaNet`
- `WGCNA`
- `igraph`
- `mediation`
- `dplyr`
- `tidyr`
- `tibble`
- `purrr`

### Organ Analysis

- `glmnet`
- `caret`
- `lavaan`
- `semPlot`
- `psych`
- `qgraph`
- `bootnet`
- `corrplot`
- `RColorBrewer`
- `ggplot2`
- `dplyr`
- `tidyr`
- `stringr`

### AR Risk Model and Stratification

- `glmnet`
- `caret`
- `FSA`
- `RColorBrewer`
- `umap`
- `cluster`
- `factoextra`
- `ggplot2`
- `dplyr`
- `tidyr`
- `vegan`
- `Metrics`
- `ModelMetrics`

### Additional R Note

- Some environments may also require `NMF` if related exploratory work or extensions are added locally, but it is not currently documented as a direct dependency in the updated script interfaces.

## Python Dependencies by Module

### Organ Analysis

- `pandas`
- `numpy`
- `scikit-learn`
- `joblib`
- `tqdm`

### AR Risk Model and Stratification

- `pandas`
- `numpy`
- `scikit-learn`
- `matplotlib`
- `shap`

## External Tools

The multi-omics integration shell workflow depends on:

- `QIIME2`
- `mmvec`
- `biom`

These tools must be available in the local environment before running `Scripts/04.Multi-omics integration/02.mmVec-model.sh`.

## Current Execution Status

- Scripts now support command-line execution where applicable.
- Original `file.choose()` behavior is preserved in relevant R scripts as the interactive fallback.
- Example data and full smoke tests will be added separately.
- Users should expect to validate local paths, package availability, and input formats before running module-level analyses.

## Multiple-testing Correction

Documentation across the repository refers to multiple-testing adjustment as **Benjamini-Hochberg (BH) correction**.

In R, `p.adjust(..., method = "fdr")` corresponds to **BH-adjusted FDR**. Standardizing the label to BH does not change the underlying calculation.
