#!/usr/bin/env Rscript

# Sync the curated integrated ceRNA network TSV into manuscript source-data
# locations. The original Excel workbook is archived outside the public release;
# the compact TSV is the self-contained GitHub input/output for Figure 4D and
# Supplementary Table X.

source(file.path(getwd(), "R", "final_helpers.R"))
load_pkgs(c("dplyr"))

network_tsv <- p("data", "processed", "final", "integrated_ceRNA_network_curated.tsv")
if (!file.exists(network_tsv)) {
  stop("Missing curated integrated network TSV: ", network_tsv)
}

network <- read_tsv(network_tsv)
required <- c(
  "miRNA", "target", "score", "hyb_energy", "target_type", "gene_symbol",
  "coordinate_status", "coordinate_interpretation"
)
missing <- setdiff(required, names(network))
if (length(missing)) stop("Curated integrated network is missing required columns: ", paste(missing, collapse = ", "))

write_tsv(network, p("results", "final", "source_data", "Supplementary_Table_integrated_ceRNA_network.tsv"))
write_tsv(network, p("results", "final", "networks", "integrated_ceRNA_network_curated.tsv"))

candidate_edges <- network |>
  dplyr::filter(!is.na(candidate_axis), nzchar(candidate_axis)) |>
  dplyr::transmute(
    candidate_axis,
    miRNA,
    target,
    target_type,
    gene_symbol,
    score,
    hybridization_energy = hyb_energy,
    original_coordinate_columns,
    coordinate_1,
    coordinate_2,
    coordinate_3,
    coordinate_4,
    coordinate_status,
    source_file = "data/processed/final/integrated_ceRNA_network_curated.tsv",
    experimentally_evaluated,
    experimental_result,
    interpretation = dplyr::case_when(
      candidate_axis == "miR-140-3p / CEMIP / circPAPPA" & target_type == "mRNA" ~ "primary biologically prioritized and experimentally supported candidate edge",
      candidate_axis == "miR-140-3p / CEMIP / circPAPPA" & target_type == "circRNA" ~ "predicted circRNA partner; functional support not established if knockdown remains negative",
      TRUE ~ "additional experimentally evaluated candidate edge with negative or less supportive result"
    )
  )
write_tsv(candidate_edges, p("data", "processed", "final", "networks", "curated_integrated_candidate_axis_edges.tsv"))
write_tsv(candidate_edges, p("results", "final", "source_data", "curated_integrated_candidate_axis_edges.tsv"))

message("Integrated network source-data copies refreshed.")
