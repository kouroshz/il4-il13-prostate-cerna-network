#!/usr/bin/env Rscript

# Publication-facing circRNA differential-expression workflow.
# Uses recovered original circRNA DEG tables as the manuscript Figure 2 source.

source(file.path(getwd(), "R", "final_helpers.R"))
load_pkgs(c("dplyr", "ggplot2", "ggrepel", "openxlsx"))

table_dir <- ensure_dir(p("results", "tables"))
figure_dir <- ensure_dir(p("results", "figures"))

metadata <- read_tsv(p("data", "metadata", "circrna_samples.tsv"))
required_meta <- c("sample_name", "source_column", "condition", "biological_replicate")
missing_meta <- setdiff(required_meta, names(metadata))
if (length(missing_meta)) stop("Missing circRNA metadata columns: ", paste(missing_meta, collapse = ", "))

norm <- openxlsx::read.xlsx(p("data", "processed", "circrna", "CircRNA_normalized_Data_KZ.xlsx"))
if (!"circRNA" %in% names(norm)) stop("circRNA normalized matrix must contain circRNA")
if (!all(metadata$source_column %in% names(norm))) {
  stop("circRNA metadata references missing expression columns: ", paste(setdiff(metadata$source_column, names(norm)), collapse = ", "))
}

read_original_deg <- function(treatment) {
  contrast_name <- paste0(treatment, "_vs_Vehicle")
  tsv_rel <- file.path("data", "processed", "circrna", paste0("circRNA_original_", contrast_name, "_DEG.tsv"))
  tsv_path <- p(tsv_rel)
  if (!file.exists(tsv_path)) stop("Missing recovered original circRNA DEG input: ", tsv_path)
  x <- read_tsv(tsv_path)
  required <- c("circRNA", "circRNA_type", "GeneSymbol", "P", "log10FC", "reg", "sig")
  missing <- setdiff(required, names(x))
  if (length(missing)) stop("Missing original circRNA DEG columns: ", paste(missing, collapse = ", "))
  x |>
    dplyr::mutate(
      raw_p_value = as.numeric(P),
      original_log10_fold_change = as.numeric(log10FC),
      log2_fold_change = original_log10_fold_change / log10(2),
      original_direction = tolower(as.character(reg)),
      original_sig_flag = tolower(as.character(sig)) == "yes",
      source_input_file = tsv_rel
    ) |>
    dplyr::arrange(circRNA, raw_p_value) |>
    dplyr::distinct(circRNA, .keep_all = TRUE)
}

circ_label <- function(alias, feature_id, gene_symbol) {
  relevant <- gene_symbol %in% c("CEMIP", "LIFR", "TNC", "PAPPA")
  dplyr::case_when(
    relevant ~ paste0("circ", gene_symbol),
    TRUE ~ ""
  )
}

run_contrast <- function(treatment) {
  contrast_name <- paste0(treatment, "_vs_Vehicle")
  meta <- metadata |>
    dplyr::filter(condition %in% c("Vehicle", treatment)) |>
    dplyr::mutate(condition = factor(condition, levels = c("Vehicle", treatment))) |>
    dplyr::arrange(condition, biological_replicate)
  vehicle_cols <- meta$source_column[meta$condition == "Vehicle"]
  treatment_cols <- meta$source_column[meta$condition == treatment]
  ordered_cols <- c(vehicle_cols, treatment_cols)

  original <- read_original_deg(treatment)
  expr_ann <- norm |>
    dplyr::select(
      feature_id = circRNA,
      probeID,
      alias = Alias,
      source,
      chrom,
      strand,
      txStart,
      txEnd,
      best_transcript,
      circRNA_length,
      dplyr::all_of(ordered_cols)
    ) |>
    dplyr::distinct(feature_id, .keep_all = TRUE)

  out <- original |>
    dplyr::transmute(
      assay = "circRNA microarray",
      analysis_method = "recovered_original_circRNA_DEG_Treatment_minus_Vehicle",
      contrast = contrast_name,
      feature_id = circRNA,
      gene_symbol = GeneSymbol,
      circRNA_type,
      original_log10_fold_change,
      log2_fold_change,
      raw_p_value,
      BH_FDR = stats::p.adjust(raw_p_value, method = "BH"),
      original_direction,
      original_sig_flag,
      source_input_file
    ) |>
    dplyr::left_join(expr_ann, by = "feature_id") |>
    dplyr::mutate(
      direction = dplyr::case_when(
        original_direction %in% c("up", "down") ~ original_direction,
        TRUE ~ direction_from_logfc(log2_fold_change)
      ),
      nominal_discovery = !is.na(raw_p_value) & raw_p_value < 0.05 &
        !is.na(log2_fold_change) & abs(log2_fold_change) >= threshold_log2fc,
      FDR_significant = !is.na(BH_FDR) & BH_FDR < 0.05 &
        !is.na(log2_fold_change) & abs(log2_fold_change) >= threshold_log2fc,
      label_for_plot = circ_label(alias, feature_id, gene_symbol)
    )
  out$mean_log2_expr_vehicle <- rowMeans(as.matrix(out[, vehicle_cols, drop = FALSE]), na.rm = TRUE)
  out$mean_log2_expr_treatment <- rowMeans(as.matrix(out[, treatment_cols, drop = FALSE]), na.rm = TRUE)
  out <- out |>
    dplyr::select(
      assay, analysis_method, contrast, feature_id, probeID, alias, gene_symbol,
      circRNA_type, original_log10_fold_change, log2_fold_change, raw_p_value, BH_FDR,
      direction, nominal_discovery, FDR_significant, original_direction, original_sig_flag,
      mean_log2_expr_vehicle, mean_log2_expr_treatment,
      dplyr::all_of(ordered_cols),
      source, chrom, strand, txStart, txEnd, best_transcript, circRNA_length,
      label_for_plot, source_input_file
  )

  prefix <- if (treatment == "IL4") "Figure_2B_circRNA_IL4_volcano" else "Figure_2C_circRNA_IL13_volcano"
  plot_volcano(
    out,
    paste("circRNA", treatment, "vs Vehicle"),
    file.path(figure_dir, paste0(prefix, ".pdf")),
    label_features = c("circPAPPA", "circLIFR", "circTNC", "circCEMIP")
  )
  out
}

all_tabs <- list(
  IL4_vs_Vehicle = run_contrast("IL4"),
  IL13_vs_Vehicle = run_contrast("IL13")
)

summary_counts <- dplyr::bind_rows(lapply(names(all_tabs), function(contrast) {
  tab <- all_tabs[[contrast]]
  dplyr::tibble(
    assay = "circRNA microarray",
    contrast = contrast,
    total_features = nrow(tab),
    nominal_discovery = sum(tab$nominal_discovery, na.rm = TRUE),
    nominal_up = sum(tab$nominal_discovery & tab$direction == "up", na.rm = TRUE),
    nominal_down = sum(tab$nominal_discovery & tab$direction == "down", na.rm = TRUE),
    FDR_significant = sum(tab$FDR_significant, na.rm = TRUE)
  )
}))

shared_nominal <- dplyr::inner_join(
  dplyr::filter(all_tabs$IL4_vs_Vehicle, nominal_discovery) |>
    dplyr::select(feature_id, direction_IL4 = direction),
  dplyr::filter(all_tabs$IL13_vs_Vehicle, nominal_discovery) |>
    dplyr::select(feature_id, direction_IL13 = direction),
  by = "feature_id"
) |>
  dplyr::filter(direction_IL4 == direction_IL13)

expected <- c(IL4_up = 318, IL4_down = 3, IL13_up = 42, IL13_down = 0, shared_up = 25)
observed <- c(
  IL4_up = summary_counts$nominal_up[summary_counts$contrast == "IL4_vs_Vehicle"],
  IL4_down = summary_counts$nominal_down[summary_counts$contrast == "IL4_vs_Vehicle"],
  IL13_up = summary_counts$nominal_up[summary_counts$contrast == "IL13_vs_Vehicle"],
  IL13_down = summary_counts$nominal_down[summary_counts$contrast == "IL13_vs_Vehicle"],
  shared_up = sum(shared_nominal$direction_IL4 == "up", na.rm = TRUE)
)
if (!identical(as.integer(observed), as.integer(expected))) {
  stop("circRNA count mismatch: observed ", paste(names(observed), observed, sep = "=", collapse = "; "))
}

write_tsv(
  dplyr::bind_rows(all_tabs) |>
    dplyr::arrange(contrast, raw_p_value, feature_id),
  file.path(table_dir, "Table_S2_circRNA_differential_expression.tsv")
)

save_figure(
  plot_two_set_venn(
    title = "circRNA differential-expression summary",
    left_label = "IL-4 upregulated",
    right_label = "IL-13 upregulated",
    left_only = as.integer(observed["IL4_up"] - observed["shared_up"]),
    shared = as.integer(observed["shared_up"]),
    right_only = as.integer(observed["IL13_up"] - observed["shared_up"]),
    shared_note = "Shared upregulated",
    bottom_note = paste0("Downregulated circRNAs: IL-4 = ", observed["IL4_down"], "; IL-13 = ", observed["IL13_down"])
  ),
  "Figure_2A_circRNA_overlap",
  width = 7.2,
  height = 4.6
)

message("Table S2 and Figure 2A-C outputs written.")
