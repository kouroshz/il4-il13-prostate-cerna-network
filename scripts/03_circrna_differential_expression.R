#!/usr/bin/env Rscript

# Final circRNA differential-expression workflow.
# Uses the recovered original circRNA DEG tables as the Figure 2/manuscript
# source. Arraystar provider subset tables and a limma reanalysis are retained
# as provenance/internal check outputs, not as the primary figure-count
# source.

source(file.path(getwd(), "R", "final_helpers.R"))
load_pkgs(c("dplyr", "ggplot2", "ggrepel", "pheatmap", "openxlsx", "limma"))

out_dir <- ensure_dir(p("data", "processed", "final", "circrna"))
fig_dir <- ensure_dir(p("results", "final", "figures", "circrna"))
src_dir <- ensure_dir(p("results", "final", "source_data", "circrna"))

summarise_counts <- function(tab, contrast, source_label) {
  dplyr::tibble(
    assay = "circRNA microarray",
    analysis_source = source_label,
    contrast = contrast,
    total_features = nrow(tab),
    nominal_discovery = sum(tab$nominal_discovery, na.rm = TRUE),
    nominal_up = sum(tab$nominal_discovery & tab$direction == "up", na.rm = TRUE),
    nominal_down = sum(tab$nominal_discovery & tab$direction == "down", na.rm = TRUE),
    FDR_significant = sum(tab$FDR_significant, na.rm = TRUE),
    FDR_up = sum(tab$FDR_significant & tab$direction == "up", na.rm = TRUE),
    FDR_down = sum(tab$FDR_significant & tab$direction == "down", na.rm = TRUE)
  )
}

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
  if (file.exists(tsv_path)) {
    x <- read_tsv(tsv_path)
    source_file <- tsv_rel
  } else {
    excel_name <- if (treatment == "IL4") "IL4_vs_V_DEGs.xlsx" else "IL13_vs_V_DEGs.xlsx"
    source_rel <- file.path("data", "raw", "circrna_microarray", excel_name)
    source_path <- p(source_rel)
    if (!file.exists(source_path)) {
      stop("Missing recovered original circRNA DEG input: ", tsv_path, " or ", source_path)
    }
    x <- openxlsx::read.xlsx(source_path)
    source_file <- source_rel
  }
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
      source_input_file = source_file
    ) |>
    dplyr::arrange(circRNA, raw_p_value) |>
    dplyr::distinct(circRNA, .keep_all = TRUE)
}

run_original_contrast <- function(treatment) {
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
      Sequence,
      dplyr::all_of(ordered_cols)
    ) |>
    dplyr::distinct(feature_id, .keep_all = TRUE)

  out <- original |>
    dplyr::transmute(
      assay = "circRNA microarray",
      analysis_source = "recovered_original_circRNA_DEG_workbook_Treatment_minus_Vehicle",
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
      label_for_plot = ifelse(is.na(gene_symbol) | gene_symbol == "", feature_id, gene_symbol)
    )
  out$mean_log2_expr_vehicle <- rowMeans(as.matrix(out[, vehicle_cols, drop = FALSE]), na.rm = TRUE)
  out$mean_log2_expr_treatment <- rowMeans(as.matrix(out[, treatment_cols, drop = FALSE]), na.rm = TRUE)
  out <- out |>
    dplyr::select(
      assay, analysis_source, contrast, feature_id, probeID, alias, gene_symbol,
      circRNA_type, original_log10_fold_change, log2_fold_change, raw_p_value, BH_FDR,
      direction, nominal_discovery, FDR_significant, original_direction, original_sig_flag,
      mean_log2_expr_vehicle, mean_log2_expr_treatment,
      dplyr::all_of(ordered_cols),
      source, chrom, strand, txStart, txEnd, best_transcript, circRNA_length,
      label_for_plot, source_input_file
    )

  stem <- paste0("circRNA_", contrast_name)
  write_split_de(out, out_dir, stem)
  write_tsv(
    out |>
      dplyr::transmute(
        feature_id,
        gene_symbol,
        display_label = label_for_plot,
        log2_fold_change,
        original_log10_fold_change,
        raw_p_value,
        BH_FDR,
        nominal_discovery,
        FDR_significant,
        direction
      ),
    file.path(src_dir, paste0(stem, "_volcano_source_data.tsv"))
  )
  plot_volcano(out, paste("circRNA original", contrast_name), file.path(fig_dir, paste0(stem, "_volcano.pdf")))
  out
}

run_limma_check <- function(treatment) {
  contrast_name <- paste0(treatment, "_vs_Vehicle")
  meta <- metadata |>
    dplyr::filter(condition %in% c("Vehicle", treatment)) |>
    dplyr::mutate(condition = factor(condition, levels = c("Vehicle", treatment))) |>
    dplyr::arrange(condition, biological_replicate)
  vehicle_cols <- meta$source_column[meta$condition == "Vehicle"]
  treatment_cols <- meta$source_column[meta$condition == treatment]
  ordered_cols <- c(vehicle_cols, treatment_cols)

  expr <- as.matrix(norm[, ordered_cols, drop = FALSE])
  mode(expr) <- "numeric"
  row_ids <- paste0("row", seq_len(nrow(expr)))
  rownames(expr) <- row_ids

  group <- factor(ifelse(meta$condition == "Vehicle", "Vehicle", "Treatment"), levels = c("Vehicle", "Treatment"))
  design <- stats::model.matrix(~0 + group)
  colnames(design) <- levels(group)
  contrast <- limma::makeContrasts(Treatment_minus_Vehicle = Treatment - Vehicle, levels = design)
  fit <- limma::lmFit(expr, design)
  fit2 <- limma::eBayes(limma::contrasts.fit(fit, contrast))

  stats <- limma::topTable(fit2, number = Inf, adjust.method = "BH", sort.by = "none") |>
    tibble::rownames_to_column("row_id") |>
    dplyr::rename(log2_fold_change = logFC, raw_p_value = P.Value, BH_FDR = adj.P.Val)

  ann <- norm |>
    dplyr::mutate(row_id = row_ids) |>
    dplyr::transmute(
      row_id,
      feature_id = circRNA,
      probeID,
      alias = Alias,
      source,
      chrom,
      strand,
      txStart,
      txEnd,
      circRNA_type,
      best_transcript,
      gene_symbol = GeneSymbol,
      sequence = Sequence,
      circRNA_length,
      dplyr::across(dplyr::all_of(ordered_cols))
    )

  out <- stats |>
    dplyr::left_join(ann, by = "row_id")
  out$mean_log2_expr_vehicle <- rowMeans(as.matrix(out[, vehicle_cols, drop = FALSE]), na.rm = TRUE)
  out$mean_log2_expr_treatment <- rowMeans(as.matrix(out[, treatment_cols, drop = FALSE]), na.rm = TRUE)
  out |>
    dplyr::mutate(
      assay = "circRNA microarray",
      analysis_source = "internal_check_limma_on_Arraystar_normalized_log2_matrix",
      contrast = contrast_name,
      manual_log2FC_from_group_means = mean_log2_expr_treatment - mean_log2_expr_vehicle,
      label_for_plot = ifelse(is.na(gene_symbol) | gene_symbol == "", feature_id, gene_symbol)
    ) |>
    classify_de() |>
    dplyr::select(
      assay, analysis_source, contrast, feature_id, probeID, alias, gene_symbol,
      log2_fold_change, manual_log2FC_from_group_means, AveExpr, t, B,
      raw_p_value, BH_FDR, direction, nominal_discovery, FDR_significant,
      mean_log2_expr_vehicle, mean_log2_expr_treatment,
      dplyr::all_of(ordered_cols),
      source, chrom, strand, txStart, txEnd, circRNA_type, best_transcript,
      circRNA_length, label_for_plot
    )
}

run_arraystar_provider <- function(sheet, contrast_name) {
  provider_file <- p("data", "processed", "circrna", "Differentially_Expressed_CircRNAs_Arraystar_provider.xlsx")
  if (!file.exists(provider_file)) {
    message("Arraystar provider workbook not found; skipping provider provenance tables.")
    return(data.frame())
  }
  x <- openxlsx::read.xlsx(provider_file, sheet = sheet)
  vehicle_mean_col <- "group-V(normalized)"
  treatment_mean_col <- grep("^group-IL.*\\(normalized\\)$", names(x), value = TRUE)
  treatment_mean_col <- setdiff(treatment_mean_col, vehicle_mean_col)[1]
  stopifnot(vehicle_mean_col %in% names(x), treatment_mean_col %in% names(x))

  provider_fc <- as.numeric(x$FC)
  provider_reg <- tolower(as.character(x$Regulation))
  provider_vehicle_vs_treatment_log2FC_from_FC <- dplyr::case_when(
    provider_reg == "up" ~ log2(provider_fc),
    provider_reg == "down" ~ -log2(provider_fc),
    TRUE ~ NA_real_
  )
  treatment_vs_vehicle_log2FC_from_FC <- -provider_vehicle_vs_treatment_log2FC_from_FC
  treatment_vs_vehicle_log2FC_from_group_means <- as.numeric(x[[treatment_mean_col]]) - as.numeric(x[[vehicle_mean_col]])

  out <- x |>
    dplyr::transmute(
      assay = "circRNA microarray",
      analysis_source = "Arraystar_provider_DE_subset_direction_converted_Treatment_minus_Vehicle",
      provider_contrast = sheet,
      contrast = contrast_name,
      feature_id = circRNA,
      probeID,
      alias = Alias,
      gene_symbol = GeneSymbol,
      provider_regulation_original = Regulation,
      provider_FC_original = as.numeric(FC),
      provider_vehicle_vs_treatment_log2FC_from_FC = provider_vehicle_vs_treatment_log2FC_from_FC,
      treatment_vs_vehicle_log2FC_from_FC = treatment_vs_vehicle_log2FC_from_FC,
      treatment_vs_vehicle_log2FC_from_group_means = treatment_vs_vehicle_log2FC_from_group_means,
      log2_fold_change = treatment_vs_vehicle_log2FC_from_group_means,
      vehicle_mean_normalized_log2 = as.numeric(.data[[vehicle_mean_col]]),
      treatment_mean_normalized_log2 = as.numeric(.data[[treatment_mean_col]]),
      raw_p_value = as.numeric(P),
      BH_FDR = as.numeric(FDR),
      source,
      chrom,
      strand,
      txStart,
      txEnd,
      circRNA_type,
      best_transcript,
      circRNA_length,
      label_for_plot = ifelse(is.na(GeneSymbol) | GeneSymbol == "", circRNA, GeneSymbol)
    ) |>
    classify_de() |>
    dplyr::mutate(
      treatment_vs_vehicle_log2FC = log2_fold_change,
      treatment_vs_vehicle_direction = direction
    ) |>
    dplyr::select(
      assay, analysis_source, provider_contrast, contrast, feature_id, probeID, alias, gene_symbol,
      provider_regulation_original, treatment_vs_vehicle_log2FC, treatment_vs_vehicle_direction,
      provider_FC_original, provider_vehicle_vs_treatment_log2FC_from_FC,
      treatment_vs_vehicle_log2FC_from_FC, treatment_vs_vehicle_log2FC_from_group_means,
      vehicle_mean_normalized_log2, treatment_mean_normalized_log2,
      raw_p_value, BH_FDR, direction, nominal_discovery, FDR_significant,
      source, chrom, strand, txStart, txEnd, circRNA_type, best_transcript, circRNA_length, label_for_plot
    )

  stem <- paste0("circRNA_Arraystar_provider_", contrast_name)
  write_split_de(out, out_dir, stem)
  out
}

make_concordance <- function(provider, limma, contrast_name) {
  if (!nrow(provider)) return(data.frame())
  provider2 <- provider |>
    dplyr::transmute(
      feature_id,
      gene_symbol_provider = gene_symbol,
      provider_contrast,
      provider_regulation_original,
      provider_treatment_vs_vehicle_log2FC = treatment_vs_vehicle_log2FC,
      provider_direction = treatment_vs_vehicle_direction,
      provider_raw_p_value = raw_p_value,
      provider_BH_FDR = BH_FDR,
      provider_nominal_discovery = nominal_discovery,
      provider_FDR_significant = FDR_significant
    )
  limma2 <- limma |>
    dplyr::transmute(
      feature_id,
      gene_symbol_limma = gene_symbol,
      limma_treatment_vs_vehicle_log2FC = log2_fold_change,
      limma_direction = direction,
      limma_raw_p_value = raw_p_value,
      limma_BH_FDR = BH_FDR,
      limma_nominal_discovery = nominal_discovery,
      limma_FDR_significant = FDR_significant
    )
  dplyr::inner_join(provider2, limma2, by = "feature_id") |>
    dplyr::mutate(
      contrast = contrast_name,
      same_treatment_vs_vehicle_direction = provider_direction == limma_direction,
      both_nominal_discovery = provider_nominal_discovery & limma_nominal_discovery,
      both_FDR_significant = provider_FDR_significant & limma_FDR_significant
    ) |>
    dplyr::select(contrast, dplyr::everything())
}

original_tabs <- list(
  IL4_vs_Vehicle = run_original_contrast("IL4"),
  IL13_vs_Vehicle = run_original_contrast("IL13")
)

limma_tabs <- list(
  IL4_vs_Vehicle = run_limma_check("IL4"),
  IL13_vs_Vehicle = run_limma_check("IL13")
)

for (contrast in names(limma_tabs)) {
  stem <- paste0("circRNA_limma_", contrast)
  write_split_de(limma_tabs[[contrast]], out_dir, stem)
}

provider_tabs <- list(
  IL4_vs_Vehicle = run_arraystar_provider("V_vs_IL4", "IL4_vs_Vehicle"),
  IL13_vs_Vehicle = run_arraystar_provider("V_vs_IL-13", "IL13_vs_Vehicle")
)

concordance <- dplyr::bind_rows(
  make_concordance(provider_tabs$IL4_vs_Vehicle, limma_tabs$IL4_vs_Vehicle, "IL4_vs_Vehicle"),
  make_concordance(provider_tabs$IL13_vs_Vehicle, limma_tabs$IL13_vs_Vehicle, "IL13_vs_Vehicle")
)
write_tsv(concordance, file.path(out_dir, "circRNA_Arraystar_provider_vs_limma_concordance.tsv"))

summary_counts <- dplyr::bind_rows(lapply(names(original_tabs), function(contrast) {
  summarise_counts(original_tabs[[contrast]], contrast, "recovered_original_circRNA_DEG_workbook")
}))
write_tsv(summary_counts, file.path(out_dir, "circRNA_summary_counts.tsv"))
write_tsv(summary_counts, file.path(out_dir, "circRNA_original_workflow_summary_counts.tsv"))

limma_summary <- dplyr::bind_rows(lapply(names(limma_tabs), function(contrast) {
  summarise_counts(limma_tabs[[contrast]], contrast, "internal_check_limma_on_Arraystar_normalized_log2_matrix")
}))
write_tsv(limma_summary, file.path(out_dir, "circRNA_limma_summary_counts.tsv"))

provider_summary <- dplyr::bind_rows(lapply(names(provider_tabs), function(contrast) {
  summarise_counts(provider_tabs[[contrast]], contrast, "Arraystar_provider_DE_subset_direction_converted")
}))
write_tsv(provider_summary, file.path(out_dir, "circRNA_Arraystar_provider_summary_counts.tsv"))

shared_nominal <- dplyr::inner_join(
  dplyr::filter(original_tabs$IL4_vs_Vehicle, nominal_discovery) |>
    dplyr::select(feature_id, alias, gene_symbol, direction_IL4 = direction, log2FC_IL4 = log2_fold_change, P_IL4 = raw_p_value, FDR_IL4 = BH_FDR),
  dplyr::filter(original_tabs$IL13_vs_Vehicle, nominal_discovery) |>
    dplyr::select(feature_id, direction_IL13 = direction, log2FC_IL13 = log2_fold_change, P_IL13 = raw_p_value, FDR_IL13 = BH_FDR),
  by = "feature_id"
) |>
  dplyr::filter(direction_IL4 == direction_IL13) |>
  dplyr::mutate(shared_direction = direction_IL4)
write_tsv(shared_nominal, file.path(out_dir, "circRNA_shared_direction_nominal_discovery.tsv"))

shared_fdr <- dplyr::inner_join(
  dplyr::filter(original_tabs$IL4_vs_Vehicle, FDR_significant) |>
    dplyr::select(feature_id, alias, gene_symbol, direction_IL4 = direction, log2FC_IL4 = log2_fold_change, P_IL4 = raw_p_value, FDR_IL4 = BH_FDR),
  dplyr::filter(original_tabs$IL13_vs_Vehicle, FDR_significant) |>
    dplyr::select(feature_id, direction_IL13 = direction, log2FC_IL13 = log2_fold_change, P_IL13 = raw_p_value, FDR_IL13 = BH_FDR),
  by = "feature_id"
) |>
  dplyr::filter(direction_IL4 == direction_IL13) |>
  dplyr::mutate(shared_direction = direction_IL4)
write_tsv(shared_fdr, file.path(out_dir, "circRNA_shared_direction_FDR_significant.tsv"))

key_features <- dplyr::bind_rows(original_tabs) |>
  dplyr::filter(grepl("KIAA1199|CEMIP|TNC|LIFR|PAPPA", paste(feature_id, alias, gene_symbol), ignore.case = TRUE))
write_tsv(key_features, file.path(out_dir, "circRNA_key_feature_validation.tsv"))

for (contrast in names(original_tabs)) {
  tab <- original_tabs[[contrast]]
  top <- tab |>
    dplyr::filter(nominal_discovery) |>
    dplyr::arrange(raw_p_value) |>
    dplyr::slice_head(n = 50)
  expr_cols <- metadata$source_column
  expr_cols <- expr_cols[expr_cols %in% names(top)]
  if (nrow(top) > 1 && length(expr_cols) > 1) {
    mat <- as.matrix(top[, expr_cols, drop = FALSE])
    mode(mat) <- "numeric"
    rownames(mat) <- make.unique(paste(top$gene_symbol, top$feature_id, sep = "|"))
    write_tsv(cbind(feature_id = rownames(mat), as.data.frame(mat)), file.path(src_dir, paste0("circRNA_", contrast, "_heatmap_source_data.tsv")))
    grDevices::pdf(file.path(fig_dir, paste0("circRNA_", contrast, "_top_nominal_heatmap.pdf")), width = 7, height = 8)
    pheatmap::pheatmap(mat, scale = "row", main = paste("circRNA top nominal", contrast), fontsize_row = 5)
    grDevices::dev.off()
    grDevices::png(file.path(fig_dir, paste0("circRNA_", contrast, "_top_nominal_heatmap.png")), width = 2100, height = 2400, res = 300)
    pheatmap::pheatmap(mat, scale = "row", main = paste("circRNA top nominal", contrast), fontsize_row = 5)
    grDevices::dev.off()
  }
}

message("Final circRNA outputs generated from recovered original DEG tables: ", out_dir)
