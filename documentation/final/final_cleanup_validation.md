# Final Cleanup Validation

Date: 2026-07-07

## Public Release Size

`github_release/` contains 181 files before rerunning generated outputs.

- Size from `du -sk`: 48,376 KiB
- Byte-sum size: 46.832 MB
- Files > 5 MB: 4
- Files > 25 MB: 0

Files flagged > 5 MB:

- `data/processed/final/circrna/circRNA_IL13_vs_Vehicle_complete.tsv` (6.404 MB)
- `data/processed/final/circrna/circRNA_IL4_vs_Vehicle_complete.tsv` (6.378 MB)
- `data/processed/final/mrna/mRNA_IL13_vs_Vehicle_complete.tsv` (6.192 MB)
- `data/processed/final/mrna/mRNA_IL4_vs_Vehicle_complete.tsv` (6.175 MB)

These are canonical complete final differential-expression tables retained once under `data/processed/final/`. No public release file exceeds 25 MB.

## Cleanup Completed

- Removed the full internal `documentation/final/recomputed_vs_accepted_v2.tsv` from `github_release/`; retained only `documentation/final/recomputed_vs_accepted_v2_summary.tsv`.
- Removed duplicate top-level `results/source_data/`, `results/networks/`, and `results/tables/` public copies; canonical release outputs now live under `results/final/`.
- Kept complete DE tables once under `data/processed/final/` and kept minimal figure source-data tables under `results/final/source_data/`.
- Removed duplicate mature miRNA FASTA copies. The only mature miRNA FASTA in the public release is `data/processed/miranda_inputs/mature_hsa.fa`.
- Confirmed the mature miRNA FASTA sequence lines contain 0 `T` bases and 14,707 `U` bases.
- Confirmed no FASTQ, FQ, BAM, SAM, `.RData`, `.rds`, `.zip`, `.DS_Store`, or generated PDF files remain in the initial public `github_release/` tree.

## Clean-Copy Validation

Validation was run from a temporary clean copy outside the project tree.

Command:

```sh
Rscript run_all.R
```

Result: completed successfully.

Post-run checks in the temporary copy:

- Size after regenerating outputs: 72,512 KiB
- File count after regenerating outputs: 290
- Files > 25 MB: 0
- Files > 5 MB: 6
- Raw/cache file types (`FASTQ`, `FQ`, `BAM`, `SAM`, `.RData`, `.rds`, `.zip`): 0
- Absolute local path hits in workflow code/config checked: 0, excluding the project name string `IL4_IL13_scirep`
- Mature miRNA FASTA files in initial release: `data/processed/miranda_inputs/mature_hsa.fa` only

## Reproduced Counts

Clean-copy `run_all.R` reproduced the revised manuscript computational counts:

- miRNA IL4 vs Vehicle: 45 nominal, 20 up, 25 down, 0 FDR significant
- miRNA IL13 vs Vehicle: 28 nominal, 15 up, 13 down, 0 FDR significant
- circRNA IL4 vs Vehicle: 321 nominal, 318 up, 3 down, 0 FDR significant
- circRNA IL13 vs Vehicle: 42 nominal, 42 up, 0 down, 0 FDR significant
- mRNA IL4 vs Vehicle: 769 nominal, 385 up, 384 down, 536 FDR significant
- mRNA IL13 vs Vehicle: 697 nominal, 359 up, 338 down, 445 FDR significant

## Final Strategy Checks

- miRNA volcano source data include explicit `biologically_prioritized`, `experimentally_evaluated`, and `candidate_highlight` columns for miR-140-3p, miR-135b-5p, and miR-625-3p.
- No hidden manual threshold change is used for candidate miRNA display; candidate highlighting is an explicit category.
- The integrated network table is included as `results/final/source_data/Supplementary_Table_integrated_ceRNA_network.tsv` and `results/final/networks/integrated_ceRNA_network_curated.tsv`.
- `results/final/source_data/curated_integrated_candidate_axis_edges.tsv` contains the five prioritized candidate edges with miRanda scores, predicted binding energies, and coordinate-status fields.
- Original coordinate-like workbook fields are preserved but not publicly interpreted as binding locations unless Quentin confirms their meaning.
- The final README describes one revised manuscript workflow and does not describe a conflicting competing-analysis framework.
