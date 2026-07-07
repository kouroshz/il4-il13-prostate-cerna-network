# Final Core Results Summary

## Analysis Definitions

The final public workflow recomputes or regenerates manuscript source data from compact processed inputs:

- mRNA-seq: `data/processed/mrna/all_counts.xlsx`
- miRNA-seq: `data/processed/mirna/Expression_Browser_CPM.csv`
- circRNA microarray: `data/processed/circrna/circRNA_original_IL4_vs_Vehicle_DEG.tsv`, `data/processed/circrna/circRNA_original_IL13_vs_Vehicle_DEG.tsv`, and `data/processed/circrna/CircRNA_normalized_Data_KZ.xlsx`
- Integrated network: `data/processed/final/integrated_ceRNA_network_curated.tsv`

The independently cultured and treated biological replicate is the statistical unit.

Final significance/candidate classes:

- mRNA `nominal_discovery`: nominal P < 0.05 and absolute fold change >= 1.5.
- circRNA `nominal_discovery`: nominal P < 0.05 and absolute fold change >= 1.5.
- miRNA `nominal_discovery`: nominal P < 0.05.
- `FDR_significant`: BH FDR < 0.05, with the 1.5-fold rule retained for mRNA/circRNA but not forced for miRNA.
- `biologically_prioritized`: selected using direction, predicted interaction, fibrosis relevance, biological rationale, and suitability for experimental testing.

FDR values are retained in complete tables for transparency but are not the primary candidate-generation gate.

## Final Feature Counts

| assay | contrast | nominal discovery | nominal up | nominal down | FDR significant | FDR up | FDR down |
|---|---|---:|---:|---:|---:|---:|---:|
| mRNA-seq | IL4 vs Vehicle | 769 | 385 | 384 | 536 | 317 | 219 |
| mRNA-seq | IL13 vs Vehicle | 697 | 359 | 338 | 445 | 281 | 164 |
| miRNA-seq | IL4 vs Vehicle | 45 | 20 | 25 | 0 | 0 | 0 |
| miRNA-seq | IL13 vs Vehicle | 28 | 15 | 13 | 0 | 0 | 0 |
| circRNA original workflow | IL4 vs Vehicle | 321 | 318 | 3 | 0 | 0 | 0 |
| circRNA original workflow | IL13 vs Vehicle | 42 | 42 | 0 | 0 | 0 | 0 |

## Shared-Feature Definition

Shared same-direction features must pass the stated threshold in both IL4-vs-Vehicle and IL13-vs-Vehicle and have the same Treatment-minus-Vehicle direction in both contrasts.

| assay | threshold | shared same-direction features | shared up | shared down |
|---|---|---:|---:|---:|
| mRNA-seq | nominal_discovery | 546 | 303 | 243 |
| mRNA-seq | FDR_significant | 379 | 258 | 121 |
| miRNA-seq | nominal_discovery | 10 | 3 | 7 |
| miRNA-seq | FDR_significant | 0 | 0 | 0 |
| circRNA original workflow | nominal_discovery | 25 | 25 | 0 |
| circRNA original workflow | FDR_significant | 0 | 0 | 0 |

## Key mRNA Values

CEMIP, TNC, and LIFR have positive Treatment-minus-Vehicle mRNA fold changes after both IL-4 and IL-13 and are FDR significant in both contrasts.

| gene | contrast | log2FC | nominal P | BH FDR | nominal discovery | FDR significant |
|---|---|---:|---:|---:|---|---|
| CEMIP | IL4 vs Vehicle | 2.672904 | 2.638673e-08 | 1.634226e-06 | yes | yes |
| CEMIP | IL13 vs Vehicle | 2.642033 | 1.671424e-07 | 1.115317e-05 | yes | yes |
| TNC | IL4 vs Vehicle | 3.801322 | 5.738259e-43 | 3.824687e-40 | yes | yes |
| TNC | IL13 vs Vehicle | 3.490508 | 2.467849e-37 | 1.729098e-34 | yes | yes |
| LIFR | IL4 vs Vehicle | 1.590842 | 4.323631e-27 | 1.476045e-24 | yes | yes |
| LIFR | IL13 vs Vehicle | 1.524361 | 2.342462e-47 | 2.984084e-44 | yes | yes |

## Candidate miRNAs

miR-140-3p, miR-135b-5p, and miR-625-3p are highlighted in volcano source data using an explicit `biologically_prioritized` and `candidate_highlight` category. This category is separate from nominal-P coloring and is not a hidden relaxed threshold.

| miRNA | IL4 log2FC | IL4 P | IL4 nominal | IL13 log2FC | IL13 P | IL13 nominal | interpretation |
|---|---:|---:|---|---:|---:|---|---|
| miR-140-3p | -0.292238 | 0.022247 | yes | -0.180649 | 0.091811 | no | primary biologically prioritized/tested CEMIP candidate |
| miR-135b-5p | -1.468787 | 0.020109 | yes | -1.766442 | 0.004561 | yes | additional evaluated TNC/circLIFR candidate |
| miR-625-3p | -0.616516 | 0.038658 | yes | -0.377775 | 0.049382 | yes | additional evaluated LIFR candidate |

## Integrated Network

The curated integrated IL-4/IL-13 ceRNA network table contains 270 miRanda-predicted interactions and is the Figure 4D/Supplementary Table X source. The Figure 4D network is exploratory, prediction-based, and biologically filtered around fibrosis-associated genes CEMIP, TNC, and LIFR.

| candidate edge | miRanda score | predicted binding energy | coordinate status | source |
|---|---:|---:|---|---|
| miR-140-3p / CEMIP | 152 | -16.17 | original workbook coordinate fields preserved; interpretation pending Quentin confirmation | `data/processed/final/integrated_ceRNA_network_curated.tsv` |
| miR-140-3p / circPAPPA | 147 | -18.27 | original workbook coordinate fields preserved; interpretation pending Quentin confirmation | `data/processed/final/integrated_ceRNA_network_curated.tsv` |
| miR-135b-5p / TNC | 152 | -11.55 | original workbook coordinate fields preserved; interpretation pending Quentin confirmation | `data/processed/final/integrated_ceRNA_network_curated.tsv` |
| miR-135b-5p / circLIFR | 143 | -9.43 | original workbook coordinate fields preserved; interpretation pending Quentin confirmation | `data/processed/final/integrated_ceRNA_network_curated.tsv` |
| miR-625-3p / LIFR | 157 | -14.26 | original workbook coordinate fields preserved; interpretation pending Quentin confirmation | `data/processed/final/integrated_ceRNA_network_curated.tsv` |

Edges should be interpreted as miRanda-predicted interactions, not as experimentally validated physical competition for every edge.
