# Code And Data Availability

This repository is self-contained from compact processed counts and normalized expression matrices onward. It includes:

- compact mRNA count matrix
- verified miRNA CPM expression matrix
- Arraystar normalized circRNA matrix and separate provider workbook
- explicit sample metadata for each assay
- final regenerated differential-expression, enrichment, miRanda, and restricted network outputs
- figure source-data files
- R scripts, configuration, `renv.lock`, `sessionInfo.txt`, and a one-command workflow

The public statistical workflow does not use precomputed release-candidate differential-expression tables as inputs. Differential expression is recomputed from the compact matrices by `scripts/01_mrna_differential_expression.R`, `scripts/02_mirna_differential_expression.R`, and `scripts/03_circrna_differential_expression.R`.

The individual featureCounts files are not included in the public release. If available locally, they can be converted into the compact mRNA count matrix with:

```sh
Rscript scripts/00_build_mrna_count_matrix_from_featurecounts.R
```

Historical browser-based g:Profiler exports and the final Cytoscape `.cys` session were not recovered locally. Recovered restricted miRanda outputs and small historical-provenance records are included. Complete miRanda rerun remains pending complete target resources and a compatible miRanda v3.3a binary.

The manuscript records these deposition targets: miRNA SRA BioProject `PRJNA1425456`, mRNA SRA BioProject `PRJNA1425902`, and circRNA GEO accession `GSE324030`. Confirm final public release status and reviewer tokens before manuscript submission.
