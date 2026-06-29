# Multi-omics Analysis in Allergic Rhinitis

<p align="center">
  <img src="Figures/Overview.png" width="950">
</p>

<p align="center">
Overview of the allergic rhinitis multi-omics study.
</p>

## Overview

This repository contains **allergic rhinitis multi-omics analysis scripts** and related documentation.

The current repository covers these main analysis modules:

- Proteome
- Metabolome and lipidome
- Multi-omics integration
- Organ analysis
- AR risk model and stratification

## Documentation

Project details are documented in:

- [Pipeline overview](docs/pipeline_overview.md)
- [Scripts manifest](docs/scripts_manifest.csv)
- [Environment notes](docs/environment_notes.md)

## Usage Overview

- Scripts now support command-line execution where applicable.
- Original `file.choose()` behavior is preserved as the interactive fallback in relevant R scripts.
- Output file names are preserved for backward compatibility.
- Some intermediate files still require manual formatting, reshaping, or annotation before downstream steps can be run.

## Reproducibility Notes

- Statistical models and core parameters were not intentionally changed during the current interface cleanup.
- Benjamini-Hochberg (BH) correction is used for P-value adjustment where applicable.
- Exact package versions were not fully pinned in the original repository.
- Scripts still require local R/Python validation with real or example data.
- Example input files and smoke tests will be added separately.

## Contact

Email: xjy005351@siat.ac.cn. For questions, collaborations, or bug reports, please open an issue or contact via email.

## License

This project is released under the MIT License.
