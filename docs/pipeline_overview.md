# Pipeline Overview

This repository provides analysis scripts for an allergic rhinitis multi-omics study. It is not a fully automated end-to-end workflow: scripts can now be run from the command line where applicable, original `file.choose()` behavior is preserved for interactive use, output file names are preserved for backward compatibility, and statistical models plus core parameters were not intentionally changed.

Example input files will be added separately after script interfaces are finalized. Scripts still require local R/Python environment validation with real or example data.

## Main Modules

### 1. Proteome

Goal: normalize Olink/proteome data, test AR vs HC differences, and run KEGG GSEA.

Recommended order:
1. `01.Olink bridging normalization.R`
2. `02.Proteomic data normalization.R`
3. `03.Differential protein analysis.R`
4. `04.GSEA enrichment.R`

Key manual handoff:
- `01.Olink bridging normalization.R` writes a bridge-normalized Olink table, but the downstream normalization step still expects a curated `sample x protein` matrix. The Olink long-format export usually needs manual reshaping first.
- `02.Proteomic data normalization.R` writes `Agesex_adjusted_protein_levels.csv`, but downstream differential analysis usually also needs `Group` merged back into the adjusted table.

### 2. Metabolome and Lipidome

Goal: normalize metabolite/lipid data and test AR vs HC differences.

Recommended order:
- Metabolome
  1. `01.Metabolic data normalization.R`
  2. `02.Differential metabolite analysis.R`
- Lipidome
  1. `03.Lipidomic data normalization.R`
  2. `04.Differential lipids analysis.R`

Key manual handoff:
- `Agesex_adjusted_metabolites.csv` and `Agesex_adjusted_lipidites.csv` usually need `Group` merged back in before differential analysis.
- The differential scripts still assume the first two columns are metadata and features start from column 3.

### 3. Multi-omics Integration

Goal: build cross-omics associations, detect modules, extract signaling triplets, and test mediation.

Recommended order:
1. `01.DIABOLO-model.R`
2. `02.mmVec-model.sh`
3. `03.Module assignment and  eigenvalue calculation.R`
4. `04.Cross-omics subpath isolation.R`
5. `05.Mediation effect analysis.R`

Key manual handoff:
- `01.DIABOLO-model.R` exports `DIABLO_network_edges.csv`, but `02.mmVec-model.sh` does not directly export `edges_from_mmVec.csv`; mmVec outputs still need manual conversion into an edge list.
- `03.Module assignment and  eigenvalue calculation.R` expects prepared `edges_from_DIABOLO.csv`, `edges_from_mmVec.csv`, `multi_omics_signatures.csv`, and `feature_type.csv`.
- `Node_module.csv` is not yet the final node table for subpath extraction. It still needs feature-type / group annotation to become `multi-omics network.nodes`.
- `04.Cross-omics subpath isolation.R` requires a manually consolidated `multi-omics network.edge` and `multi-omics network.nodes`.
- `05.Mediation effect analysis.R` expects a `sample x feature` signature table, which may require manual reformatting if an upstream step used `feature x sample`.

### 4. Organ Analysis

Goal: estimate organ-level risk indices, select organ-linked proteins, compare organ correlation networks, and fit a structural equation model.

Recommended order:
1. `01.Organ risk index.py`
2. `02.LASSO analysis.R`
3. `03.GLASSO analysis.R`
4. `04.Structural Equations Model.R`

Key manual handoff:
- `01.Organ risk index.py` produces organ-level outputs used downstream, but these still rely on curated protein matrices, protein-to-organ mappings, and metadata.
- `02.LASSO analysis.R` writes `Top_Lasso_Features_by_Organ.csv`, but no standardized downstream matrix is created automatically.
- The organ branch still uses the `1HC` / `2AR` label convention, which differs from the `HC` / `AR` convention used by the risk-model scripts.

### 5. AR Risk Model and Stratification

Goal: train omics-specific and multi-omics classifiers, rank features with SHAP, perform feature selection, and stratify patients.

Recommended order:
1. `01.AR risk model-Micorbiome.py`
2. `02.AR risk model-Metabolome.py`
3. `03.AR risk model-Lipidome.py`
4. `04.AR risk model-Proteome.py`
5. `05.Feature selection.R`
6. `06.AR risk model-Multi-omics.py`
7. `07.AR stratification.R`

Key manual handoff:
- The AR risk model scripts expect curated target tables and do not construct them from upstream normalization outputs.
- `05.Feature selection.R` writes `Lasso_feature.csv`, but this is a feature-selection output, not the final input matrix for downstream classifiers or stratification.
- `07.AR stratification.R` expects a curated numeric matrix and does not read `Lasso_feature.csv` directly.

## P-value Correction

Across the repository, multiple-testing adjustment is documented as **Benjamini-Hochberg (BH) correction**.

In R, `p.adjust(..., method = "fdr")` corresponds to **BH-adjusted FDR**. Where scripts were standardized, the explicit `method = "BH"` form is used without changing the calculation result.

## Current Interface Status

- Scripts now support command-line execution where applicable.
- Original `file.choose()` behavior is preserved as the interactive fallback in relevant R scripts.
- Output file names are preserved for backward compatibility.
- Statistical models, thresholds, and core parameters were not intentionally changed during interface cleanup.

## Important Practical Note

This repository should currently be treated as a set of analysis scripts plus documented interfaces, not as a one-click reproducible workflow. The most common points requiring manual intervention are:

- reshaping Olink long-format output into a `sample x protein` matrix;
- merging `Group` back into `Agesex_adjusted_*` tables before differential analysis;
- converting DIABOLO and mmVec outputs into compatible edge lists;
- augmenting `Node_module.csv` into `multi-omics network.nodes`;
- constructing curated target matrices for risk modeling and stratification.
