# Computational Figure Caption Drafts

## Figure 1

Corrected miRNA-seq analysis of N1 fibroblasts treated with Vehicle, IL-4, or IL-13 using three independently cultured and treated biological replicates per condition. Points in volcano plots represent individual miRNAs tested by two-sided Welch tests on log2(CPM + 1) values. Nominal-discovery miRNAs meet nominal P < 0.05; no miRNAs were FDR significant under the final BH FDR < 0.05 criterion. Shared lists require the same miRNA to meet the threshold in both cytokine comparisons with the same direction of change. Labeled candidate miRNAs were prioritized for experimental evaluation based on integrated ceRNA network position, predicted interactions with fibrosis-relevant targets, and downstream validation experiments, rather than fold-change magnitude alone.

## Figure 2

Recovered original circRNA microarray analysis from three biological replicates per condition. Points in volcano plots represent individual circRNAs from the recovered original DEG tables with Treatment-minus-Vehicle fold changes. Heatmaps display Arraystar normalized log2 values for selected nominal-discovery circRNAs. Enrichment panels show gprofiler2 analyses of unique host genes for upregulated nominal-discovery circRNAs using all unique host genes represented on the analyzed Arraystar platform as background; no circRNA host-gene terms passed adjusted P < 0.05 in the final reproducible analysis.

## Figure 3

Corrected mRNA-seq analysis of three biological replicates per condition using edgeR. Points in volcano plots represent tested Ensembl genes after CPM filtering and TMM normalization. Nominal-discovery genes meet nominal P < 0.05 and absolute fold change >= 1.5, and FDR-significant genes meet BH FDR < 0.05 and absolute fold change >= 1.5. Enrichment panels show gprofiler2 analyses of upregulated nominal-discovery mRNAs using all retained tested genes as background.

## Figure 4

Curated integrated IL-4/IL-13 ceRNA candidate network. Network edges indicate miRanda-predicted miRNA-target interactions with minimum alignment score 140. The biologically focused pro-fibrotic subnetwork is centered on CEMIP, LIFR, and TNC and highlights prioritized axes selected for experimental evaluation. Edge weights are not shown; miRanda scores and predicted binding energies are provided in Supplementary Table X and the GitHub source-data table. The network is exploratory and hypothesis-generating, and predicted edges should not be described as experimentally proven physical competition unless directly validated.

## Figure 5

Experimental evaluation of selected biologically prioritized candidate interactions. Add panel-specific biological n, technical-replicate handling, statistical unit, error-bar definition, and statistical test from Jill's records. Figure 5A should label miR-140-3p, CEMIP mRNA, and circPAPPA / hsa-circ_0002052 if shown, and should distinguish predicted components from experimentally supported components.
