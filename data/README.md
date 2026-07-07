# Data Directory

This directory contains compact processed inputs and regenerated final outputs.

## Included Inputs

- `metadata/mrna_samples.tsv`: mRNA sample metadata.
- `metadata/mirna_samples.tsv`: miRNA sample metadata.
- `metadata/circrna_samples.tsv`: circRNA sample metadata.
- `processed/mrna/all_counts.xlsx`: compact mRNA count matrix.
- `processed/mirna/Expression_Browser_CPM.csv`: verified miRNA CPM expression matrix.
- `processed/circrna/CircRNA_normalized_Data_KZ.xlsx`: Arraystar normalized circRNA log2 expression matrix.
- `processed/circrna/circRNA_original_IL4_vs_Vehicle_DEG.tsv`: recovered original IL-4-versus-Vehicle circRNA DEG table used as the Figure 2 source.
- `processed/circrna/circRNA_original_IL13_vs_Vehicle_DEG.tsv`: recovered original IL-13-versus-Vehicle circRNA DEG table used as the Figure 2 source.
- `processed/circrna/Differentially_Expressed_CircRNAs_Arraystar_provider.xlsx`: Arraystar provider differential-expression subset retained as separate provenance.
- `processed/miranda_inputs/`: recovered restricted miRanda input FASTA files.
- `processed/miranda_outputs/`: recovered restricted miRanda output files.

## Generated Final Outputs

- `processed/final/mrna/`: final edgeR mRNA differential-expression tables and summaries.
- `processed/final/mirna/`: final miRNA Welch-test tables, overlaps, and miR-140-3p values.
- `processed/final/circrna/`: final recovered-original circRNA tables plus provider/limma provenance tables.
- `processed/final/enrichment_inputs/`: exact gprofiler2 input and background lists.
- `processed/final/miranda/`: copied miRanda inputs, focused pending input, and parsed predictions.
- `processed/final/networks/`: canonical restricted network node and edge tables.

## Excluded Raw Files

The individual featureCounts files are excluded from the public release. Use `scripts/00_build_mrna_count_matrix_from_featurecounts.R` to reconstruct the compact count matrix if the raw count files are available locally.
