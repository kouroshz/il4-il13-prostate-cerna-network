#!/usr/bin/env Rscript

# Publication-facing summary counts table.

source(file.path(getwd(), "R", "final_helpers.R"))
load_pkgs(c("dplyr"))

table_dir <- ensure_dir(p("results", "tables"))

mirna <- read_tsv(file.path(table_dir, "Table_S1_miRNA_differential_expression.tsv"))
circrna <- read_tsv(file.path(table_dir, "Table_S2_circRNA_differential_expression.tsv"))
mrna <- read_tsv(file.path(table_dir, "Table_S3_mRNA_differential_expression.tsv"))

count_contrast <- function(tab, assay, contrast) {
  x <- tab |> dplyr::filter(contrast == !!contrast)
  dplyr::tibble(
    assay = assay,
    contrast = contrast,
    metric = c("total_features", "nominal_discovery", "nominal_up", "nominal_down", "FDR_significant"),
    value = c(
      nrow(x),
      sum(x$nominal_discovery, na.rm = TRUE),
      sum(x$nominal_discovery & x$direction == "up", na.rm = TRUE),
      sum(x$nominal_discovery & x$direction == "down", na.rm = TRUE),
      sum(x$FDR_significant, na.rm = TRUE)
    )
  )
}

shared_same_direction <- function(tab, assay) {
  x <- dplyr::inner_join(
    tab |> dplyr::filter(contrast == "IL4_vs_Vehicle", nominal_discovery) |> dplyr::select(feature_id, direction_IL4 = direction),
    tab |> dplyr::filter(contrast == "IL13_vs_Vehicle", nominal_discovery) |> dplyr::select(feature_id, direction_IL13 = direction),
    by = "feature_id"
  ) |>
    dplyr::filter(direction_IL4 == direction_IL13)
  dplyr::tibble(
    assay = assay,
    contrast = "shared_IL4_IL13",
    metric = "shared_same_direction_nominal_discovery",
    value = nrow(x)
  )
}

summary_counts <- dplyr::bind_rows(
  count_contrast(mirna, "miRNA-seq", "IL4_vs_Vehicle"),
  count_contrast(mirna, "miRNA-seq", "IL13_vs_Vehicle"),
  shared_same_direction(mirna, "miRNA-seq"),
  count_contrast(circrna, "circRNA microarray", "IL4_vs_Vehicle"),
  count_contrast(circrna, "circRNA microarray", "IL13_vs_Vehicle"),
  shared_same_direction(circrna, "circRNA microarray"),
  count_contrast(mrna, "mRNA-seq", "IL4_vs_Vehicle"),
  count_contrast(mrna, "mRNA-seq", "IL13_vs_Vehicle"),
  shared_same_direction(mrna, "mRNA-seq")
)

expected <- dplyr::tibble(
  assay = c("miRNA-seq", "miRNA-seq", "miRNA-seq", "circRNA microarray", "circRNA microarray", "circRNA microarray", "circRNA microarray", "circRNA microarray", "mRNA-seq", "mRNA-seq", "mRNA-seq"),
  contrast = c("IL4_vs_Vehicle", "IL13_vs_Vehicle", "shared_IL4_IL13", "IL4_vs_Vehicle", "IL4_vs_Vehicle", "IL13_vs_Vehicle", "IL13_vs_Vehicle", "shared_IL4_IL13", "IL4_vs_Vehicle", "IL13_vs_Vehicle", "shared_IL4_IL13"),
  metric = c("nominal_discovery", "nominal_discovery", "shared_same_direction_nominal_discovery", "nominal_up", "nominal_down", "nominal_up", "nominal_down", "shared_same_direction_nominal_discovery", "nominal_discovery", "nominal_discovery", "shared_same_direction_nominal_discovery"),
  expected_value = c(45, 28, 10, 318, 3, 42, 0, 25, 769, 697, 546)
)
check <- expected |>
  dplyr::left_join(summary_counts, by = c("assay", "contrast", "metric")) |>
  dplyr::mutate(matches_expected = value == expected_value)
if (any(!check$matches_expected, na.rm = TRUE) || any(is.na(check$value))) {
  stop("Summary count mismatch: ", paste(check$assay, check$contrast, check$metric, check$value, check$expected_value, sep = "=", collapse = "; "))
}

write_tsv(
  summary_counts |>
    dplyr::left_join(expected, by = c("assay", "contrast", "metric")) |>
    dplyr::mutate(matches_expected = dplyr::if_else(is.na(expected_value), NA, value == expected_value)),
  file.path(table_dir, "Table_S7_summary_counts.tsv")
)

message("Table S7 output written.")
