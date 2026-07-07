# Data Directory

This directory contains compact processed inputs required by the public workflow.

## Included Inputs

- `metadata/mrna_samples.tsv`: mRNA sample metadata.
- `metadata/mirna_samples.tsv`: miRNA sample metadata.
- `metadata/circrna_samples.tsv`: circRNA sample metadata.
- `processed/mrna/all_counts.xlsx`: compact mRNA count matrix.
- `processed/mirna/Expression_Browser_CPM.csv`: verified miRNA CPM expression matrix.
- `processed/circrna/CircRNA_normalized_Data_KZ.xlsx`: Arraystar normalized circRNA log2 expression matrix.
- `processed/circrna/circRNA_original_IL4_vs_Vehicle_DEG.tsv`: recovered original IL-4-versus-Vehicle circRNA DEG table used as the Figure 2 source.
- `processed/circrna/circRNA_original_IL13_vs_Vehicle_DEG.tsv`: recovered original IL-13-versus-Vehicle circRNA DEG table used as the Figure 2 source.
- `processed/circrna/Differentially_Expressed_CircRNAs_Arraystar_provider.xlsx`: Arraystar provider differential-expression subset retained as a supporting check input.
- `processed/miranda_inputs/`: recovered miRanda input FASTA files.
- `processed/miranda_outputs/`: recovered miRanda output files.
- `processed/network/integrated_ceRNA_network_curated.tsv`: compact curated integrated ceRNA network input.

## Generated Outputs

Generated workflow outputs are written under `results/tables/`, `results/figures/`, `results/logs/`, and `results/session/`. Those folders are kept empty in the repository except for `.gitkeep` placeholders.

## Excluded Raw Files

The individual featureCounts files are excluded from the public release. Use `scripts/00_build_mrna_count_matrix_from_featurecounts.R` to reconstruct the compact count matrix if the raw count files are available locally.
