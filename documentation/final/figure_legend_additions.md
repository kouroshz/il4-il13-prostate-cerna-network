# Figure-Legend Additions

## General Computational Figure Language

For mRNA-seq, miRNA-seq, and circRNA microarray panels, add: data were generated from three independently cultured and treated biological replicates per condition. Each library or microarray represents one biological replicate, and the biological replicate is the statistical unit. Fold changes are Treatment minus Vehicle. mRNA nominal-discovery features are defined as edgeR likelihood-ratio-test P < 0.05 and absolute fold change >= 1.5; circRNA nominal-discovery features are defined as recovered original workflow P < 0.05 and absolute fold change >= 1.5; miRNA nominal-discovery features are defined as Welch-test P < 0.05. Benjamini-Hochberg FDR values are reported in complete supplementary tables.

## Figure 1

Volcano plot points represent individual miRNAs tested in the corrected log2(CPM + 1) two-sided Welch-test analysis. Colored points should distinguish FDR-significant features, nominal-discovery features, and features not meeting threshold. State that no miRNAs were FDR significant under the final threshold. Shared-list panels should state that overlap requires the same miRNA to meet the threshold in both cytokine comparisons with the same direction of change.

## Figure 2

Volcano plot points represent individual circRNAs from the recovered original circRNA differential-expression tables, with fold changes expressed as Treatment minus Vehicle. Heatmaps show selected nominal-discovery circRNAs and Arraystar normalized log2 expression values. GO and KEGG panels should state that host genes of upregulated nominal-discovery circRNAs were analyzed with gprofiler2 using all unique host genes represented on the analyzed platform as background. If non-significant top terms are displayed, state that no plotted circRNA host-gene terms passed adjusted P < 0.05.

## Figure 3

Volcano plot points represent individual Ensembl genes tested by edgeR. PCA or sample-correlation panels, where shown, use log2(CPM + 1) values derived from retained mRNA features. GO and KEGG panels should state that upregulated nominal-discovery mRNAs were analyzed with gprofiler2 using all genes retained for differential-expression testing as background.

## Figure 4

Network nodes represent cytokine-regulated miRNAs, mRNAs, and circRNAs included in the curated integrated IL-4/IL-13 candidate network. Edges indicate miRanda-predicted miRNA-target interactions. Edge weights are not shown; corresponding miRanda scores and predicted binding energies are provided in Supplementary Table X and the GitHub repository. Panel D represents a biologically focused pro-fibrotic subnetwork centered on CEMIP, LIFR, and TNC. The network is exploratory and hypothesis-generating, and predicted edges should not be described as experimentally proven physical competition unless directly supported by validation data.

## Figure 5

Add panel-specific wet-lab details from Jill: biological n, technical-replicate handling, statistical unit, error-bar definition, and statistical test for each miRNA knockdown, miRNA overexpression, western blot, qRT-PCR, RNA pulldown, and circPAPPA panel. For Figure 5A, explicitly label miRNA: miR-140-3p; mRNA: CEMIP mRNA; circRNA: circPAPPA / hsa-circ_0002052. State that circPAPPA was the predicted circRNA partner for the miR-140-3p/CEMIP axis but was not functionally supported if the circPAPPA knockdown remained negative. Do not state that every predicted network edge was experimentally validated.
