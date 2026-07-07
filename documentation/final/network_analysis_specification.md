# Network Analysis Specification

## Scope

The canonical network source for the revised manuscript package is the curated integrated IL-4/IL-13 ceRNA interaction table:

`data/processed/final/integrated_ceRNA_network_curated.tsv`

This table was parsed from Jill/Quentin's integrated Excel workbook and preserves the original miRNA, target, miRanda score, predicted hybridization energy, coordinate-like workbook fields, and target-type fields while adding standardized columns for downstream documentation. The coordinate-like fields are preserved for Quentin review but should not be described as verified coordinate fields in public wording. The table should be treated as the curated interaction table underlying the final pro-fibrotic shared Cytoscape network, especially Figure 4D.

The raw IL-4 and IL-13 miRanda DOCX files are supporting provenance and are incomplete relative to the curated integrated Excel. They can be used to cross-check selected interactions, but they are not the sole authoritative network source.

## Figure 4D Interpretation

The full predicted network was biologically filtered to interactions connected to the pro-fibrotic genes CEMIP, TNC, and LIFR. This focused subnetwork nominated three experimentally evaluated candidate axes:

- miR-140-3p / CEMIP / circPAPPA
- miR-135b-5p / TNC / circLIFR
- miR-625-3p / LIFR

Network diagrams were generated in Cytoscape from curated miRanda-predicted interactions. The network is exploratory and hypothesis-generating. Edges represent predicted miRNA-target interactions, not experimentally proven physical competition for every edge. miRanda scores and predicted binding energies are provided in the integrated network source table and supplementary/source-data table; they are not shown as edge weights in the Cytoscape figure.

## Canonical Files

`data/processed/final/integrated_ceRNA_network_curated.tsv`: curated integrated IL-4/IL-13 ceRNA interaction table with standardized fields.

`documentation/final/candidate_axis_summary.tsv`: focused summary of the three prioritized axes and their candidate edges, including scores, energies, coordinate-status audit result, source file, experimental status, and interpretation.

`results/final/source_data/Supplementary_Table_integrated_ceRNA_network.tsv`: manuscript/source-data copy of the curated integrated network table.

`github_release/results/final/networks/integrated_ceRNA_network_curated.tsv`: compact public GitHub release network-output copy.

`github_release/results/final/source_data/Supplementary_Table_integrated_ceRNA_network.tsv`: public supplementary/source-data copy.

## Inclusion Rules

mRNA and circRNA nominal discovery is nominal P < 0.05 plus absolute fold change >= 1.5. miRNA nominal discovery is nominal P < 0.05, with candidate prioritization additionally using directionality, predicted interaction, fibrosis relevance, and experimental suitability. FDR values are retained for transparency. Nominal-discovery features may enter exploratory candidate networks, but they must not be described as genome-wide significant unless they also satisfy the FDR criterion.

Predicted edges are hypothesis-generating. Experimental support for a candidate does not retroactively make the discovery-assay miRNA or target statistically significant.
