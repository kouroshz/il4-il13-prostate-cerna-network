# IL-4/IL-13 Prostate Fibroblast ceRNA Analysis

This repository contains the public, reproducible computational workflow for the revised IL-4/IL-13 prostate fibroblast ceRNA analysis. It includes compact processed inputs, sample metadata, R code, configuration files, and scripts that regenerate the manuscript computational figures and supplementary/source tables.

## Quick Start

From a fresh clone:

```sh
git clone https://github.com/kouroshz/il4-il13-prostate-cerna-network.git
cd il4-il13-prostate-cerna-network
```

Restore the recorded R package environment:

```sh
Rscript -e 'options(repos = c(CRAN = "https://cloud.r-project.org")); if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv"); renv::restore(prompt = FALSE)'
```

Run the full workflow:

```sh
Rscript run_all.R
```

The workflow starts from the compact processed inputs in `data/processed/` and writes regenerated outputs to `results/`.

## Optional Makefile

GNU Make is optional. The Makefile is a convenience wrapper for users who prefer Make-based commands:

```sh
make run
```

`make run` runs the same workflow as `Rscript run_all.R`. The Makefile also includes lightweight helper targets for validation and cleaning generated outputs; it does not change analysis settings or thresholds.

## Expected Outputs

After a successful run:

- `results/figures/` contains manuscript-numbered PNG and PDF figures.
- `results/tables/` contains supplementary/source-data TSV tables.
- `results/logs/` contains workflow and miRanda input-check logs.
- `results/session/` contains the R session information for reproducibility.

The generated computational figures are:

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

The generated supplementary/source tables are:

- `Table_S1_miRNA_differential_expression.tsv`
- `Table_S2_circRNA_differential_expression.tsv`
- `Table_S3_mRNA_differential_expression.tsv`
- `Table_S4_enrichment_results.tsv`
- `Table_S5_integrated_ceRNA_network.tsv`
- `Table_S6_prioritized_axes_miRanda_scores.tsv`
- `Table_S7_summary_counts.tsv`

## Repository Contents

- `config/`: analysis settings, thresholds, contrasts, and paths.
- `data/metadata/`: sample metadata used by the workflow.
- `data/processed/`: compact processed matrices, recovered circRNA source tables, miRanda inputs/outputs, and curated network inputs.
- `R/`: shared helper functions.
- `scripts/`: numbered workflow scripts for differential expression, enrichment, miRanda input checks, network output, and summary counts.
- `results/`: empty output folders in the repository; populated by `Rscript run_all.R`.
- `run_all.R`: one-command workflow driver.
- `renv.lock`: recorded R package environment.
- `CODE_AND_DATA_AVAILABILITY.md`: data-source and workflow-scope summary.

## Analysis Scope

- mRNA: edgeR with TMM normalization, GLM fitting, and likelihood-ratio tests.
- miRNA: `log2(CPM + 1)` values with Welch tests.
- circRNA: recovered Arraystar differential-expression source tables used as the manuscript Figure 2 source.
- Enrichment: mRNA enrichment panels for Figure 3D/E.
- ceRNA network: current revised integrated network using miRanda-predicted interactions and prioritized axes.

The recovered miRanda outputs were generated with miRanda v3.3a using mature miRNA query FASTA files with RNA `U` bases, reference `all_3putr_circ.fa.txt`, gap open penalty -9, gap extend penalty -4, minimum miRanda alignment score ≥ 140, energy threshold 1.0 kcal/mol, and scaling parameter 4. Source tables retain miRanda score, energy, and site-count information; target-position-like fields from recovered miRanda output are not treated as genomic coordinates.

The one-command workflow does not perform a de novo transcriptome-wide miRanda rerun. It uses the recovered miRanda outputs and curated integrated network source tables included in the repository.

## Data Availability

Raw sequencing and array data are deposited separately under the manuscript accessions:

- miRNA SRA BioProject: `PRJNA1425456`
- mRNA SRA BioProject: `PRJNA1425902`
- circRNA GEO accession: `GSE324030`

This repository starts from compact processed inputs rather than raw FASTQ, featureCounts, or full raw microarray export files.

## Citation and License

Citation metadata are provided in `CITATION.cff`. License terms are provided in `LICENSE`.
