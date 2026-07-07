# IL-4/IL-13 Prostate Fibroblast ceRNA Analysis

## Purpose

This repository contains code, compact processed inputs, and a reproducible workflow for the revised IL-4/IL-13 prostate fibroblast ceRNA analysis.

## Repository contents

- `config/`: analysis settings, thresholds, contrasts, and input/output paths.
- `data/metadata/`: sample metadata for the mRNA, miRNA, and circRNA analyses.
- `data/processed/`: compact processed inputs required to rerun the public workflow.
- `R/`: shared R helper functions.
- `scripts/`: assay-specific analysis, enrichment, miRanda, network, and summary scripts.
- `results/`: empty output folders preserved with `.gitkeep`; populated by running the workflow.
- `run_all.R`: one-command workflow driver.
- `renv.lock`: recorded R package environment.
- `CODE_AND_DATA_AVAILABILITY.md`: data-source and workflow-scope summary.

## Requirements

- R
- `renv`
- GNU Make, optional

## Reproducibility

Restore the R environment:

```sh
Rscript -e 'install.packages("renv"); renv::restore()'
```

Run the workflow:

```sh
Rscript run_all.R
```

Optional Make target:

```sh
make run
```

## Generated outputs

The workflow generates:

- manuscript-numbered figure files in `results/figures/`
- supplementary/source tables in `results/tables/`
- logs in `results/logs/`
- session information in `results/session/`

## Analysis summary

- mRNA: edgeR with TMM normalization, GLM fitting, and likelihood-ratio tests.
- miRNA: `log2(CPM + 1)` values with Welch tests.
- circRNA: recovered Arraystar differential-expression source tables.
- ceRNA network: miRanda-predicted interactions from the current revised integrated network, including prioritized axes.

## Data availability

- miRNA SRA BioProject: `PRJNA1425456`
- mRNA SRA BioProject: `PRJNA1425902`
- circRNA GEO accession: `GSE324030`

## License and citation

License terms are provided in `LICENSE`. Citation metadata are provided in `CITATION.cff`.
