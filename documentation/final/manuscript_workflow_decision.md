# Final Manuscript Workflow Decision

## Decision

The public GitHub workflow is now organized as one revised manuscript-aligned workflow, not as a historical analysis plus a competing improved analysis.

Where the recovered original workflow is valid, it is used as the final manuscript source. Genuine errors or ambiguities are corrected in code, source data, and manuscript-ready text drafts, and the final figure/text counts are expected to match the corrected outputs.

## Assay-Specific Rules

- mRNA-seq: the recovered edgeR GLM workflow is used for final Figure 3 values. Exploratory nominal-discovery genes are nominal P < 0.05 with absolute fold change >= 1.5. BH FDR is calculated and exported for every tested gene.
- miRNA-seq: the verified CPM matrix is analyzed with log2(CPM + 1) and two-sided Welch tests. The exploratory screen is nominal P < 0.05. The 1.5-fold cutoff is not forced for miRNA candidate highlighting because the manuscript candidate choices were network- and validation-driven. BH FDR is retained in complete tables.
- circRNA microarray: the recovered original circRNA DEG workbooks are used as the final Figure 2 source. Exploratory nominal-discovery circRNAs are nominal P < 0.05 with absolute fold change >= 1.5. BH FDR is computed from the recovered P values for transparency.
- Network: `IL4 IL13 Integrated ceRNA Network.xlsx`, parsed to `data/processed/final/integrated_ceRNA_network_curated.tsv`, is the curated Figure 4D and Supplementary Table X source.

## Candidate Selection

Candidate miRNAs are highlighted explicitly as biologically prioritized/tested candidates, not by hidden relaxed volcano thresholds. The highlighted candidates are miR-140-3p, miR-135b-5p, and miR-625-3p.

Figure captions should state that labeled candidate miRNAs were prioritized for experimental evaluation based on integrated ceRNA network position, predicted interactions with fibrosis-relevant targets, and downstream validation experiments, rather than fold-change magnitude alone.

## FDR Reporting

FDR values are calculated, retained, and exported for transparency. FDR is not the primary candidate-generation gate for this exploratory manuscript workflow.

## Provenance/Internal Checks

Arraystar provider subset tables and the limma circRNA reanalysis remain available as provenance/internal check outputs. They are not the main Figure 2 count source and should not be presented as a competing public analysis.
