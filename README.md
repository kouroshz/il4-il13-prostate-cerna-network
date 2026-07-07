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

From a fresh clone:

```sh
git clone https://github.com/kouroshz/il4-il13-prostate-cerna-network.git
cd il4-il13-prostate-cerna-network
```

Restore the R environment with a noninteractive CRAN mirror:

```sh
Rscript -e 'options(repos = c(CRAN = "https://cloud.r-project.org")); if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv"); renv::restore(prompt = FALSE)'
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

- supplementary/source tables in `results/tables/`
- manuscript-numbered figure files in `results/figures/`:
  - `Figure_1A_miRNA_overlap`
  - `Figure_1B_miRNA_IL4_volcano`
  - `Figure_1C_miRNA_IL13_volcano`
  - `Figure_2A_circRNA_overlap`
  - `Figure_2B_circRNA_IL4_volcano`
  - `Figure_2C_circRNA_IL13_volcano`
  - `Figure_3A_mRNA_overlap`
  - `Figure_3B_mRNA_IL4_volcano`
  - `Figure_3C_mRNA_IL13_volcano`
  - `Figure_3D_mRNA_IL4_enrichment`
  - `Figure_3E_mRNA_IL13_enrichment`
  - `Figure_4A_IL4_ceRNA_network`
  - `Figure_4B_IL13_ceRNA_network`
  - `Figure_4C_shared_ceRNA_network`
  - `Figure_4D_focused_ceRNA_network`
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
