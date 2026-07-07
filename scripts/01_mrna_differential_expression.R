#!/usr/bin/env Rscript

# Publication-facing mRNA-seq differential-expression workflow.
# Recomputes edgeR statistics from the compact processed count matrix.

source(file.path(getwd(), "R", "final_helpers.R"))
load_pkgs(c("dplyr", "ggplot2", "ggrepel", "openxlsx", "edgeR", "org.Hs.eg.db", "AnnotationDbi"))

table_dir <- ensure_dir(p("results", "tables"))
figure_dir <- ensure_dir(p("results", "figures"))

collapse_values <- function(x) {
  paste(unique(stats::na.omit(as.character(x))), collapse = ";")
}

annotate_ensembl <- function(feature_ids) {
  raw <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = unique(feature_ids),
    keytype = "ENSEMBL",
    columns = c("SYMBOL", "ENTREZID", "GENENAME")
  )
  raw |>
    dplyr::group_by(ENSEMBL) |>
    dplyr::summarise(
      gene_symbol = collapse_values(SYMBOL),
      entrez_id = collapse_values(ENTREZID),
      gene_name = collapse_values(GENENAME),
      .groups = "drop"
    ) |>
    dplyr::rename(feature_id = ENSEMBL)
}

metadata <- read_tsv(p("data", "metadata", "mrna_samples.tsv"))
required_meta <- c("sample_name", "source_column", "condition", "biological_replicate")
missing_meta <- setdiff(required_meta, names(metadata))
if (length(missing_meta)) stop("Missing mRNA metadata columns: ", paste(missing_meta, collapse = ", "))

counts <- openxlsx::read.xlsx(p("data", "processed", "mrna", "all_counts.xlsx"))
if (!"GeneID" %in% names(counts)) stop("mRNA count matrix must contain GeneID")
if (!all(metadata$source_column %in% names(counts))) {
  stop("mRNA metadata references missing count columns: ", paste(setdiff(metadata$source_column, names(counts)), collapse = ", "))
}

run_contrast <- function(treatment) {
  contrast_name <- paste0(treatment, "_vs_Vehicle")
  meta <- metadata |>
    dplyr::filter(condition %in% c("Vehicle", treatment)) |>
    dplyr::mutate(condition = factor(condition, levels = c("Vehicle", treatment))) |>
    dplyr::arrange(condition, biological_replicate)
  vehicle_cols <- meta$source_column[meta$condition == "Vehicle"]
  treatment_cols <- meta$source_column[meta$condition == treatment]
  if (length(vehicle_cols) != 3 || length(treatment_cols) != 3) {
    stop("Expected three Vehicle and three ", treatment, " biological replicates for mRNA")
  }

  expr <- as.matrix(counts[, meta$source_column, drop = FALSE])
  mode(expr) <- "numeric"
  rownames(expr) <- counts$GeneID

  dge <- edgeR::DGEList(counts = expr)
  keep <- rowSums(edgeR::cpm(dge) > 1) >= 2
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  group <- factor(ifelse(meta$condition == "Vehicle", "Vehicle", "Treatment"), levels = c("Vehicle", "Treatment"))
  design <- stats::model.matrix(~group)
  dge <- edgeR::calcNormFactors(dge)
  dge <- edgeR::estimateDisp(dge, design)
  fit <- edgeR::glmFit(dge, design)
  lrt <- edgeR::glmLRT(fit, coef = "groupTreatment")
  norm_cpm <- edgeR::cpm(dge, normalized.lib.sizes = TRUE)
  raw_counts <- dge$counts

  stats <- edgeR::topTags(lrt, n = Inf, adjust.method = "BH", sort.by = "none")$table |>
    tibble::rownames_to_column("feature_id") |>
    dplyr::rename(log2_fold_change = logFC, raw_p_value = PValue, BH_FDR = FDR)

  ann <- annotate_ensembl(stats$feature_id)
  ordered_cols <- c(vehicle_cols, treatment_cols)
  out <- stats |>
    dplyr::left_join(ann, by = "feature_id") |>
    dplyr::mutate(
      assay = "mRNA-seq",
      analysis_method = "edgeR_TMM_GLM_LRT",
      contrast = contrast_name,
      mean_raw_count_vehicle = rowMeans(raw_counts[feature_id, vehicle_cols, drop = FALSE]),
      mean_raw_count_treatment = rowMeans(raw_counts[feature_id, treatment_cols, drop = FALSE]),
      mean_CPM_vehicle = rowMeans(norm_cpm[feature_id, vehicle_cols, drop = FALSE]),
      mean_CPM_treatment = rowMeans(norm_cpm[feature_id, treatment_cols, drop = FALSE]),
      label_for_plot = ifelse(is.na(gene_symbol) | gene_symbol == "", feature_id, gene_symbol)
    ) |>
    dplyr::bind_cols(as.data.frame(raw_counts[stats$feature_id, ordered_cols, drop = FALSE]) |>
                       stats::setNames(paste0("raw_count_", ordered_cols))) |>
    dplyr::bind_cols(as.data.frame(norm_cpm[stats$feature_id, ordered_cols, drop = FALSE]) |>
                       stats::setNames(paste0("CPM_", ordered_cols))) |>
    classify_de() |>
    dplyr::select(
      assay, analysis_method, contrast, feature_id, gene_symbol, entrez_id, gene_name,
      log2_fold_change, logCPM, LR, raw_p_value, BH_FDR, direction,
      nominal_discovery, FDR_significant,
      mean_raw_count_vehicle, mean_raw_count_treatment,
      mean_CPM_vehicle, mean_CPM_treatment,
      dplyr::starts_with("raw_count_"), dplyr::starts_with("CPM_"), label_for_plot
    )

  prefix <- if (treatment == "IL4") "Figure_3B_mRNA_IL4_volcano" else "Figure_3C_mRNA_IL13_volcano"
  plot_volcano(out, paste("mRNA", treatment, "vs Vehicle"), file.path(figure_dir, paste0(prefix, ".pdf")))
  out
}

all_tabs <- list(
  IL4_vs_Vehicle = run_contrast("IL4"),
  IL13_vs_Vehicle = run_contrast("IL13")
)

summary_counts <- dplyr::bind_rows(lapply(names(all_tabs), function(contrast) {
  tab <- all_tabs[[contrast]]
  dplyr::tibble(
    assay = "mRNA-seq",
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
    dplyr::select(feature_id, gene_symbol, direction_IL4 = direction),
  dplyr::filter(all_tabs$IL13_vs_Vehicle, nominal_discovery) |>
    dplyr::select(feature_id, direction_IL13 = direction),
  by = "feature_id"
) |>
  dplyr::filter(direction_IL4 == direction_IL13)

expected <- c(IL4_vs_Vehicle = 769, IL13_vs_Vehicle = 697, shared_same_direction = 546)
observed <- c(
  IL4_vs_Vehicle = summary_counts$nominal_discovery[summary_counts$contrast == "IL4_vs_Vehicle"],
  IL13_vs_Vehicle = summary_counts$nominal_discovery[summary_counts$contrast == "IL13_vs_Vehicle"],
  shared_same_direction = nrow(shared_nominal)
)
if (!identical(as.integer(observed), as.integer(expected))) {
  stop("mRNA count mismatch: observed ", paste(names(observed), observed, sep = "=", collapse = "; "))
}

write_tsv(
  dplyr::bind_rows(all_tabs) |>
    dplyr::arrange(contrast, raw_p_value, feature_id),
  file.path(table_dir, "Table_S3_mRNA_differential_expression.tsv")
)

overlap_counts <- dplyr::tibble(
  category = c("IL-4 regulated", "IL-13 regulated", "Shared same direction"),
  count = as.integer(c(observed["IL4_vs_Vehicle"], observed["IL13_vs_Vehicle"], observed["shared_same_direction"]))
)
save_figure(
  plot_overlap_counts(overlap_counts, "mRNA differential expression overlap", "Nominally regulated genes"),
  "Figure_3A_mRNA_overlap",
  width = 6.2,
  height = 4.6
)

message("Table S3 and Figure 3A-C outputs written.")
