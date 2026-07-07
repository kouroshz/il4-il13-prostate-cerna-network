# Reviewer-Response Bioinformatics Notes

## Reviewer 1

We clarified the computational design by documenting sample metadata for all omics assays, the independently cultured biological replicate as the statistical unit, and the one-command reproducible workflow in the public repository. The revised workflow uses compact processed inputs, assay-specific scripts, exported source-data tables, and session information.

For mRNA-seq, the revised workflow uses the recovered edgeR GLM analysis from the compact count matrix. For miRNA-seq, the workflow uses log2(CPM + 1) values with two-sided Welch tests and nominal P < 0.05 for exploratory screening. For circRNA, the workflow uses the recovered original circRNA DEG tables as the Figure 2 source, with Arraystar provider subset tables and limma retained as provenance/internal check outputs. BH FDR values are exported in complete tables for transparency.

## Reviewer 2

The rationale for miR-140-3p selection is now stated explicitly. miR-140-3p is not presented as the strongest fold-change miRNA. It is highlighted as a biologically prioritized/tested candidate because it is directionally reduced, connects in the curated integrated ceRNA network to CEMIP and circPAPPA, CEMIP is strongly induced and FDR significant in both cytokine contrasts, and downstream functional experiments support the CEMIP relationship.

Additional evaluated candidate axes are documented: miR-135b-5p/TNC/circLIFR and miR-625-3p/LIFR. These are described as evaluated candidates with negative or less supportive results, not as validated mechanisms.

Figure 4D should be described as an exploratory, prediction-based, biologically filtered network centered on CEMIP, TNC, and LIFR. Edges indicate miRanda-predicted interactions; edge weights are not shown. miRanda scores and predicted binding energies are provided in Supplementary Table X. Original coordinate-like workbook fields are preserved in the source table but should not be described as coordinate fields unless Quentin verifies their interpretation.

Figure 5A should label miRNA: miR-140-3p, mRNA: CEMIP mRNA, and circRNA: circPAPPA / hsa-circ_0002052. If the circPAPPA knockdown remains negative, the text should state that circPAPPA is the predicted circRNA partner but was not functionally supported.
