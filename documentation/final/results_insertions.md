# Results Insertions

Use these as minimal additions or replacements within computational Results subsections.

## Opening Replicate Statement

RNA profiling was performed using three independently cultured and treated biological replicates per condition, with predefined IL-4-versus-Vehicle and IL-13-versus-Vehicle comparisons. The independently cultured and treated sample was the statistical unit for each omics assay.

## miRNA Results

Using the two-sided Welch-test analysis of log2(CPM + 1) miRNA expression values, 45 miRNAs met the exploratory nominal P < 0.05 screen after IL-4 treatment and 28 met this screen after IL-13 treatment. No miRNAs were FDR significant at BH FDR < 0.05. Ten miRNAs met the nominal-P screen in both cytokine comparisons with the same direction of change. miR-140-3p showed a modest reduction following both IL-4 and IL-13 treatment, reached nominal P < 0.05 after IL-4 but not after IL-13, and was prioritized because of its predicted interaction with strongly induced CEMIP, its relevance to fibrosis, and subsequent functional support.

## circRNA Results

The final reproducible circRNA analysis used the recovered original circRNA differential-expression tables. IL-4 treatment produced 321 nominal-discovery circRNAs, including 318 upregulated and 3 downregulated circRNAs. IL-13 treatment produced 42 nominal-discovery circRNAs, all upregulated. Twenty-five circRNAs were shared upregulated nominal-discovery features. No circRNAs were FDR significant after BH correction of the recovered P values. Arraystar provider subset tables and a limma analysis were retained as provenance/internal check outputs and were not used as the primary Figure 2 count source.

## mRNA Results

The recovered edgeR analysis identified 769 IL-4 nominal-discovery mRNAs and 697 IL-13 nominal-discovery mRNAs. Among these, 536 IL-4 and 445 IL-13 mRNAs were FDR significant. Same-direction shared mRNA responses included 546 nominal-discovery genes and 379 FDR-significant genes. CEMIP, TNC, and LIFR were positively induced by both IL-4 and IL-13 and were FDR significant in both contrasts.

## Enrichment Results

For reproducibility, enrichment was rerun programmatically with gprofiler2 using exact final input lists and assay-specific backgrounds. Upregulated nominal-discovery mRNAs produced 68 adjusted-significant enriched terms for IL-4 and 44 for IL-13. Upregulated nominal-discovery circRNA host genes did not produce adjusted-significant terms in the final gprofiler2 run; top adjusted terms are shown only for transparency and should not be described as significant pathway enrichment.

If the original circRNA enrichment panels are retained, Quentin should verify the exact input lists, background lists, g:Profiler settings, and panel interpretation against the current reproducible gprofiler2 handoff. If that verification is not available, state that the reproducible gprofiler2 rerun did not identify adjusted-significant circRNA host-gene terms.

## Network Results

The curated integrated IL-4/IL-13 ceRNA network table contains miRanda-predicted miRNA-target interactions with scores and predicted binding energies. The candidate network was filtered around the pro-fibrotic genes CEMIP, TNC, and LIFR. This identified miR-140-3p/CEMIP/circPAPPA, miR-135b-5p/TNC/circLIFR, and miR-625-3p/LIFR as prioritized axes for experimental evaluation. This directly supports Reviewer 2's request to explain why miR-140-3p was selected. The network is exploratory and hypothesis-generating: edges represent miRanda-predicted interactions, not experimentally proven physical competition for every predicted edge.
