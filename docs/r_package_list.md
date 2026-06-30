# R Package List

This file lists the main R and Bioconductor packages referenced by the current repository scripts.

These dependencies have **not** been fully validated in a single clean environment. After running analyses locally, users should record their R `sessionInfo()` and the Python package versions used for that run.

## Proteome

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

## Metabolome and Lipidome

- `robustbase`
- `vegan`
- `ggplot2`
- `ggrepel`
- `plotly`
- `grid`

## Multi-omics Integration

- `mixOmics`
- `MetaNet`
- `WGCNA`
- `igraph`
- `mediation`
- `dplyr`
- `tidyr`
- `tibble`
- `purrr`

## Organ Analysis

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

## AR Risk Model and Stratification

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

## Additional Note

- `NMF` is included here as a package to track for environment planning, even though it is not currently documented as a direct dependency in the updated script interfaces.
