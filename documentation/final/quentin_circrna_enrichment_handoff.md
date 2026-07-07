# Quentin CircRNA Enrichment Handoff

## Purpose

This handoff gives Quentin the exact reproducible circRNA host-gene enrichment inputs, backgrounds, and current gprofiler2 results so he can verify or reproduce the manuscript Figure 2 enrichment panels. The current GitHub workflow uses recovered original circRNA DEG tables for Figure 2 counts and does not depend on any external handoff materials.

## Final CircRNA Differential-Expression Inputs

- IL-4 recovered original DEG input: `data/processed/circrna/circRNA_original_IL4_vs_Vehicle_DEG.tsv`
- IL-13 recovered original DEG input: `data/processed/circrna/circRNA_original_IL13_vs_Vehicle_DEG.tsv`
- Final circRNA complete tables:
  - `data/processed/final/circrna/circRNA_IL4_vs_Vehicle_complete.tsv`
  - `data/processed/final/circrna/circRNA_IL13_vs_Vehicle_complete.tsv`
- Final count summary: `data/processed/final/circrna/circRNA_summary_counts.tsv`

Final counts:

- IL-4: 321 nominal circRNAs total, 318 up, 3 down, 0 FDR significant.
- IL-13: 42 nominal circRNAs total, 42 up, 0 down, 0 FDR significant.
- Shared same-direction nominal circRNAs: 25, all upregulated.

## Exact Enrichment Input Lists

- IL-4 upregulated nominal circRNA host-gene input list:
  - `data/processed/final/enrichment_inputs/circrna/circRNA_IL4_vs_Vehicle_up_nominal_host_genes_input_symbols.txt`
  - current unique host-gene count: 232
- IL-13 upregulated nominal circRNA host-gene input list:
  - `data/processed/final/enrichment_inputs/circrna/circRNA_IL13_vs_Vehicle_up_nominal_host_genes_input_symbols.txt`
  - current unique host-gene count: 32

## Exact Background Lists

- IL-4 circRNA host-gene background:
  - `data/processed/final/enrichment_inputs/circrna/circRNA_IL4_vs_Vehicle_up_nominal_host_genes_background_symbols.txt`
  - current background count: 6228
- IL-13 circRNA host-gene background:
  - `data/processed/final/enrichment_inputs/circrna/circRNA_IL13_vs_Vehicle_up_nominal_host_genes_background_symbols.txt`
  - current background count: 6229

The background is all unique annotated host genes represented on the analyzed Arraystar circRNA platform for the relevant contrast.

## Current Internal gprofiler2 Result

Summary file:

`results/final/enrichment/circrna/circRNA_enrichment_summary.tsv`

Current reproducible gprofiler2 result:

| contrast | analysis_id | input genes | background genes | adjusted-significant terms |
|---|---|---:|---:|---:|
| IL-4 vs Vehicle | `circRNA_IL4_vs_Vehicle_up_nominal_host_genes` | 232 | 6228 | 0 |
| IL-13 vs Vehicle | `circRNA_IL13_vs_Vehicle_up_nominal_host_genes` | 32 | 6229 | 0 |

Per-contrast output folders:

- `results/final/enrichment/circrna/circRNA_IL4_vs_Vehicle_up_nominal_host_genes/`
- `results/final/enrichment/circrna/circRNA_IL13_vs_Vehicle_up_nominal_host_genes/`

Each folder contains:

- `*_parameters.tsv`
- `*_input_list.txt`
- `*_background_list.txt`
- `*_complete_results.tsv`
- `*_significant_results.tsv`
- `*_figure_source_data.tsv`
- `*_plotted_terms.tsv`

## Current Analysis Parameters

From the `*_parameters.tsv` files:

- organism: `hsapiens`
- identifier type: HGNC gene symbol
- sources: GO:MF, GO:BP, GO:CC, KEGG
- correction method: g:Profiler `g_SCS`
- ordered query: FALSE
- domain scope: custom

## Manuscript Figure 2 Enrichment Issue

The current reproducible gprofiler2 rerun found no adjusted-significant circRNA host-gene terms. If manuscript Figure 2 currently shows significant circRNA enrichment panels from browser g:Profiler or a different historical workflow, Quentin should verify:

1. The exact uploaded IL-4 and IL-13 host-gene lists.
2. The exact background/domain setting.
3. The identifier type.
4. The g:Profiler correction method.
5. The database version/date if available.
6. Whether panels should be retained, revised, moved to supplement, or described as exploratory/non-significant top-term displays.

## What Would Require Manuscript Change

The manuscript should change if Quentin cannot reproduce adjusted-significant circRNA host-gene enrichment using the intended original inputs/settings. In that case:

- Do not describe circRNA host-gene enrichment terms as statistically significant.
- State that the reproducible gprofiler2 rerun did not identify adjusted-significant circRNA host-gene terms.
- If top terms are shown for transparency, label them explicitly as non-significant top adjusted terms.

## Files To Provide Quentin

The `handoff_to_quentin/` folder contains a compact copy of the exact input/background lists, current output summaries, curated integrated network TSV, candidate-axis summary, coordinate audit, and exact verification questions.
