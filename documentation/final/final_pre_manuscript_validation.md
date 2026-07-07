# Final Pre-Manuscript Validation

Date: 2026-07-07

## Completed This Pass

- Archived the pre-cleanup state under `archive/pre_quick_final_cleanup_2026_07_07/` with SHA-256 checksums.
- Parsed the Jill/Quentin integrated ceRNA workbook as the authoritative Figure 4D/Supplementary Table X network source.
- Created `documentation/final/figure_1_3_reproducibility_audit.tsv` for main computational Figures 1-3.
- Created `documentation/final/integrated_network_coordinate_audit.tsv` and suppressed public binding-location claims because the workbook coordinate fields are not independently interpreted.
- Created `documentation/final/candidate_axis_summary.tsv` for miR-140-3p/CEMIP/circPAPPA, miR-135b-5p/TNC/circLIFR, and miR-625-3p/LIFR.
- Created `documentation/final/quentin_circrna_enrichment_handoff.md`, `documentation/final/quentin_todo.md`, and `handoff_to_quentin/`.
- Updated manuscript-facing methods, results, figure legend, reviewer-response, README, and Quarto docs to describe one revised manuscript-aligned workflow.
- Refreshed `github_release/`, removed duplicate public result copies, and reran clean-copy validation successfully.

## Current Release Status

- `github_release/` files: 181
- Initial release size: 48,376 KiB by `du -sk` (46.832 MB byte sum)
- Files > 25 MB: 0
- Files > 5 MB: 4
- Clean-copy `Rscript run_all.R`: passed
- Mature miRNA FASTA check: 0 `T` bases and 14,707 `U` bases in sequence lines

## Micro-Cleanup Note

- Assay-specific p-value wording was corrected; no analytical outputs changed.

## Manuscript-Ready Claims

- miRanda scores and predicted binding energies from the integrated network can be reported for the five prioritized candidate edges.
- Coordinate-like workbook fields should be described only as preserved source fields until Quentin confirms their interpretation.
- CircRNA host-gene enrichment should not be described as adjusted-significant unless Quentin reproduces significance using the intended original settings.
- miRNA candidate labels should be described as biologically prioritized and experimentally evaluated, not as hidden threshold-relaxed discoveries.

## Remaining Author/Quentin Checks

- Quentin: verify circRNA host-gene enrichment inputs/settings and decide whether enrichment panels should be revised, moved, or labelled as exploratory top terms.
- Quentin: confirm whether the original coordinate columns can be interpreted as target-coordinate positions.
- Jill/Quentin: confirm final Figure 4D and Figure 5A labels for circPAPPA/hsa-circ_0002052, CEMIP, TNC, LIFR, and the tested miRNAs.
- Jill: provide final wet-lab replicate counts, statistical tests, error-bar definitions, and outcome labels for the functional panels.
