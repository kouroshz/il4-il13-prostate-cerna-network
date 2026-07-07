#!/usr/bin/env Rscript

# Write publication-facing integrated ceRNA network tables.

source(file.path(getwd(), "R", "final_helpers.R"))
load_pkgs(c("dplyr"))

table_dir <- ensure_dir(p("results", "tables"))
network_tsv <- p("data", "processed", "network", "integrated_ceRNA_network_curated.tsv")
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

write_tsv(
  network |>
    dplyr::arrange(miRNA, target_type, target),
  file.path(table_dir, "Table_S5_integrated_ceRNA_network.tsv")
)

priority_axes <- network |>
  dplyr::filter(!is.na(prioritized_axis), nzchar(prioritized_axis)) |>
  dplyr::transmute(
    prioritized_axis,
    miRNA,
    target,
    target_type,
    gene_symbol,
    miRanda_score = score,
    hybridization_energy = hyb_energy,
    coordinate_status,
    experimentally_evaluated,
    experimental_result,
    interpretation = dplyr::case_when(
      prioritized_axis == "miR-140-3p / CEMIP / circPAPPA" & target_type == "mRNA" ~ "primary biologically prioritized and experimentally supported edge",
      prioritized_axis == "miR-140-3p / CEMIP / circPAPPA" & target_type == "circRNA" ~ "predicted circRNA partner; functional support not established if knockdown remains negative",
      TRUE ~ "additional experimentally evaluated edge with negative or less supportive result"
    )
  ) |>
  dplyr::arrange(prioritized_axis, target_type, target)
write_tsv(priority_axes, file.path(table_dir, "Table_S6_prioritized_axes_miRanda_scores.tsv"))

message("Table S5 and Table S6 outputs written.")
