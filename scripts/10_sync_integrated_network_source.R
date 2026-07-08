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
  "miRNA", "target", "score", "hyb_energy", "gene.stp", "target_type", "gene_symbol",
  "prioritized_axis", "figure4D_focus", "experimentally_evaluated", "experimental_result"
)
missing <- setdiff(required, names(network))
if (length(missing)) stop("Curated integrated network is missing required columns: ", paste(missing, collapse = ", "))

clean_text <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x
}

site_count_from_target_positions <- function(x) {
  x <- clean_text(x)
  vapply(strsplit(trimws(x), "\\s+"), function(parts) {
    parts <- parts[nzchar(parts)]
    if (!length(parts)) return(NA_integer_)
    length(parts)
  }, integer(1))
}

public_experimental_result <- function(axis, target, target_type, evaluated) {
  axis <- clean_text(axis)
  target <- clean_text(target)
  evaluated <- clean_text(evaluated)
  dplyr::case_when(
    axis == "miR-140-3p / CEMIP / circPAPPA" & target == "CEMIP" ~
      "strongest functional support; pursued further",
    axis == "miR-140-3p / CEMIP / circPAPPA" & target == "circPAPPA" ~
      "predicted circRNA partner; circPAPPA knockdown did not significantly alter CEMIP under tested conditions",
    axis == "miR-135b-5p / TNC / circLIFR" & target == "circLIFR" ~
      "predicted circRNA partner in the prioritized miR-135b-5p/TNC axis",
    axis == "miR-135b-5p / TNC / circLIFR" & target == "TNC" ~
      "did not show the same supportive response as miR-140-3p/CEMIP under tested conditions",
    axis == "miR-625-3p / LIFR" & target == "LIFR" ~
      "did not show the same supportive response as miR-140-3p/CEMIP under tested conditions",
    TRUE ~ ""
  )
}

network_public <- network |>
  dplyr::transmute(
    miRNA,
    target,
    target_type,
    gene_symbol,
    miRanda_score = score,
    miRanda_energy = hyb_energy,
    site_count = site_count_from_target_positions(`gene.stp`),
    prioritized_axis = clean_text(prioritized_axis),
    figure4D_focus = clean_text(figure4D_focus),
    experimentally_evaluated = clean_text(experimentally_evaluated),
    experimental_result = public_experimental_result(prioritized_axis, target, target_type, experimentally_evaluated)
  )

write_tsv(
  network_public |>
    dplyr::arrange(miRNA, target_type, target),
  file.path(table_dir, "Table_S5_integrated_ceRNA_network.tsv")
)

priority_axes <- network_public |>
  dplyr::filter(nzchar(prioritized_axis)) |>
  dplyr::transmute(
    prioritized_axis,
    miRNA,
    target,
    target_type,
    gene_symbol,
    miRanda_score,
    miRanda_energy,
    site_count,
    experimentally_evaluated,
    experimental_result,
    interpretation = experimental_result
  ) |>
  dplyr::arrange(prioritized_axis, target_type, target)
write_tsv(priority_axes, file.path(table_dir, "Table_S6_prioritized_axes_miRanda_scores.tsv"))

message("Table S5 and Table S6 outputs written.")
