# Multi-omics Analysis in Allergic Rhinitis

<p align="center">
  <img src="Figures/Overview.png" width="950">
</p>

<p align="center">
Overview of the allergic rhinitis multi-omics study.
</p>

## Overview

This repository contains **allergic rhinitis (AR) multi-omics analysis scripts** and related documentation. The scripts cover proteomic, metabolomic, lipidomic, microbiome-related, organ-level, and integrative analyses used to characterize molecular alterations, cross-omics associations, organ-level effects, risk prediction, and molecular stratification in AR.


## Modules

The repository includes five main analysis modules:

* **Proteome analysis**
  Olink bridging normalization, proteomic normalization, differential protein analysis, and KEGG GSEA.

* **Metabolome and lipidome analysis**
  Age/sex-adjusted metabolite and lipid profiles, differential analysis, and volcano plots.

* **Multi-omics integration**
  DIABLO, mmVec-related integration, network module analysis, cross-omics subpath extraction, and mediation analysis.

* **Organ analysis**
  Organ risk index construction, LASSO feature selection, GLASSO network analysis, and structural equation modeling.

* **AR risk model and stratification**
  Single-omics and multi-omics random forest models, SHAP interpretation, LASSO feature selection, and unsupervised stratification.

## Repository Structure

```text
Figures/                                       Main overview figure and related images
Scripts/01.Proteome                            Proteome analysis scripts
Scripts/02.Metabolome and lipidome             Metabolome and lipidome analysis scripts
Scripts/03.Multi-omics integration             Multi-omics integration scripts
Scripts/04.Organ analysis                      Organ-level analysis scripts
Scripts/05.AR risk model and stratification    Risk model and stratification scripts
docs/                                          Pipeline, script manifest, and environment notes
examples/                                      Synthetic minimal example input files
requirements.txt                               Main Python dependencies
```

## Documentation

Detailed documentation is available in:

* [Pipeline overview](docs/pipeline_overview.md)
  Overview of analysis modules, recommended script order, and manually curated intermediate steps.

* [Scripts manifest](docs/scripts_manifest.csv)
  Script-level summary of inputs, outputs, command-line arguments, retained interactive behavior, key parameters, and potential issues.

* [Environment notes](docs/environment_notes.md)
  Notes on R/Python dependencies, external tools, and environment reproducibility.

* [R package list](docs/r_package_list.md)
  Module-level list of major R and Bioconductor dependencies.

* [Example data README](examples/README.md)
  Description of the synthetic minimal example input files.

## Example Data

Synthetic minimal example input files are provided in:

```text
examples/
```

These files are intended only to illustrate expected input formats, including sample-by-feature matrices, metadata tables, edge lists, node annotation files, and model input tables. The example files do not contain real participant information, real sample IDs, or real measured omics values.

## Usage

Most scripts support command-line execution where applicable, while the original `file.choose()` interactive behavior is preserved in R scripts as a fallback.

For example, differential protein analysis can be run as:

```bash
Rscript "Scripts/01.Proteome/03.Differential protein analysis.R" \
  --input examples/proteome/adjusted_proteome.csv \
  --outdir results/proteome_diff \
  --group-column Group
```

Some intermediate files still require manual formatting, reshaping, or annotation before downstream steps can be run. Examples include converting Olink long-format outputs to sample-by-protein matrices, merging group labels back into adjusted omics matrices, and preparing network edge/node files for multi-omics integration.

Please refer to the [pipeline overview](docs/pipeline_overview.md) and [scripts manifest](docs/scripts_manifest.csv) for detailed input and output requirements.

## Dependencies

The repository uses both R and Python scripts.

Python dependencies are listed in:

```text
requirements.txt
```

R and Bioconductor dependencies are summarized in:

```text
docs/r_package_list.md
```

Some integration scripts also require external tools, including:

* QIIME2
* mmvec
* biom

Please refer to [Environment notes](docs/environment_notes.md) for details.

## License

This project is released under the MIT License.

## Contact

For questions, collaborations, or bug reports, please open an issue or contact:

[xjy005351@siat.ac.cn](mailto:xjy005351@siat.ac.cn)
