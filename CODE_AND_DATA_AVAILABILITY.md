# Code and Data Availability

This repository contains compact processed inputs, sample metadata, R scripts, configuration files, and the output-folder scaffold for the IL-4/IL-13 prostate fibroblast ceRNA analysis.

The public workflow starts from compact processed matrices and recovered analysis tables, not raw FASTQ or featureCounts files. Run the workflow from the repository root with:

```sh
Rscript run_all.R
```

The repository includes `renv.lock` for package restoration. The workflow writes session information to `results/session/sessionInfo.txt`; the included root-level `sessionInfo.txt` records the validation environment for this release.

Raw sequencing and microarray data are deposited separately at the manuscript accessions:

- miRNA SRA BioProject `PRJNA1425456`
- mRNA SRA BioProject `PRJNA1425902`
- circRNA GEO accession `GSE324030`

Reviewer tokens and accession release status should be confirmed before final manuscript submission.

The repository includes the recovered miRanda prediction outputs and curated integrated network source tables used for the manuscript; de novo transcriptome-wide miRanda rerun is outside the scope of this release. The recovered miRanda outputs used miRanda v3.3a with minimum alignment score ≥ 140, energy threshold 1.0 kcal/mol, gap open penalty -9, gap extend penalty -4, and scaling parameter 4.
