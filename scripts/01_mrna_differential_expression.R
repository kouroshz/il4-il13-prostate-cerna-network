#!/usr/bin/env Rscript

# Final mRNA-seq differential-expression workflow.
# Recomputes edgeR statistics from the compact processed count matrix.

source(file.path(getwd(), "R", "final_helpers.R"))
load_pkgs(c("dplyr", "ggplot2", "ggrepel", "pheatmap", "openxlsx", "edgeR", "org.Hs.eg.db", "AnnotationDbi"))

out_dir <- ensure_dir(p("data", "processed", "final", "mrna"))
fig_dir <- ensure_dir(p("results", "final", "figures", "mrna"))
src_dir <- ensure_dir(p("results", "final", "source_data", "mrna"))

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

summarise_counts <- function(tab, contrast) {
  dplyr::tibble(
    assay = "mRNA-seq",
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
      analysis_source = "edgeR_GLM_from_compact_count_matrix_Treatment_minus_Vehicle",
      contrast = contrast_name,
      mean_raw_count_vehicle = rowMeans(raw_counts[feature_id, vehicle_cols, drop = FALSE]),
      mean_raw_count_treatment = rowMeans(raw_counts[feature_id, treatment_cols, drop = FALSE]),
      mean_CPM_vehicle = rowMeans(norm_cpm[feature_id, vehicle_cols, drop = FALSE]),
      mean_CPM_treatment = rowMeans(norm_cpm[feature_id, treatment_cols, drop = FALSE]),
      manual_log2FC_from_mean_CPM = log2((mean_CPM_treatment + 1e-8) / (mean_CPM_vehicle + 1e-8)),
      label_for_plot = ifelse(is.na(gene_symbol) | gene_symbol == "", feature_id, gene_symbol)
    ) |>
    dplyr::bind_cols(as.data.frame(raw_counts[stats$feature_id, ordered_cols, drop = FALSE]) |>
                       stats::setNames(paste0("raw_count_", ordered_cols))) |>
    dplyr::bind_cols(as.data.frame(norm_cpm[stats$feature_id, ordered_cols, drop = FALSE]) |>
                       stats::setNames(paste0("CPM_", ordered_cols))) |>
    classify_de() |>
    dplyr::select(
      assay, analysis_source, contrast, feature_id, gene_symbol, entrez_id, gene_name,
      log2_fold_change, logCPM, LR, raw_p_value, BH_FDR, direction,
      nominal_discovery, FDR_significant,
      mean_raw_count_vehicle, mean_raw_count_treatment,
      mean_CPM_vehicle, mean_CPM_treatment, manual_log2FC_from_mean_CPM,
      dplyr::starts_with("raw_count_"), dplyr::starts_with("CPM_"), label_for_plot
    )

  stem <- paste0("mRNA_", contrast_name)
  write_split_de(out, out_dir, stem)
  write_tsv(
    out |>
      dplyr::transmute(
        feature_id,
        gene_symbol,
        display_label = label_for_plot,
        log2_fold_change,
        raw_p_value,
        BH_FDR,
        nominal_discovery,
        FDR_significant,
        direction
      ),
    file.path(src_dir, paste0(stem, "_volcano_source_data.tsv"))
  )
  plot_volcano(out, paste("mRNA", contrast_name), file.path(fig_dir, paste0(stem, "_volcano.pdf")))
  out
}

all_tabs <- list(
  IL4_vs_Vehicle = run_contrast("IL4"),
  IL13_vs_Vehicle = run_contrast("IL13")
)

summary_counts <- dplyr::bind_rows(lapply(names(all_tabs), function(contrast) summarise_counts(all_tabs[[contrast]], contrast)))
write_tsv(summary_counts, file.path(out_dir, "mRNA_summary_counts.tsv"))

shared_nominal <- dplyr::inner_join(
  dplyr::filter(all_tabs$IL4_vs_Vehicle, nominal_discovery) |>
    dplyr::select(feature_id, gene_symbol, direction_IL4 = direction, log2FC_IL4 = log2_fold_change, P_IL4 = raw_p_value, FDR_IL4 = BH_FDR),
  dplyr::filter(all_tabs$IL13_vs_Vehicle, nominal_discovery) |>
    dplyr::select(feature_id, direction_IL13 = direction, log2FC_IL13 = log2_fold_change, P_IL13 = raw_p_value, FDR_IL13 = BH_FDR),
  by = "feature_id"
) |>
  dplyr::filter(direction_IL4 == direction_IL13) |>
  dplyr::mutate(shared_direction = direction_IL4)
write_tsv(shared_nominal, file.path(out_dir, "mRNA_shared_direction_nominal_discovery.tsv"))

shared_fdr <- dplyr::inner_join(
  dplyr::filter(all_tabs$IL4_vs_Vehicle, FDR_significant) |>
    dplyr::select(feature_id, gene_symbol, direction_IL4 = direction, log2FC_IL4 = log2_fold_change, P_IL4 = raw_p_value, FDR_IL4 = BH_FDR),
  dplyr::filter(all_tabs$IL13_vs_Vehicle, FDR_significant) |>
    dplyr::select(feature_id, direction_IL13 = direction, log2FC_IL13 = log2_fold_change, P_IL13 = raw_p_value, FDR_IL13 = BH_FDR),
  by = "feature_id"
) |>
  dplyr::filter(direction_IL4 == direction_IL13) |>
  dplyr::mutate(shared_direction = direction_IL4)
write_tsv(shared_fdr, file.path(out_dir, "mRNA_shared_direction_FDR_significant.tsv"))

cpm_cols <- grep("^CPM_", names(all_tabs$IL4_vs_Vehicle), value = TRUE)
cpm_mat <- all_tabs$IL4_vs_Vehicle |>
  dplyr::select(feature_id, dplyr::all_of(cpm_cols))
rownames(cpm_mat) <- cpm_mat$feature_id
cpm_mat <- as.matrix(cpm_mat[, cpm_cols, drop = FALSE])
mode(cpm_mat) <- "numeric"
log_cpm <- log2(cpm_mat + 1)

pca <- stats::prcomp(t(log_cpm), scale. = TRUE)
pca_scores <- as.data.frame(pca$x[, 1:2, drop = FALSE])
pca_scores$sample <- rownames(pca_scores)
pca_scores$source_column <- sub("^CPM_", "", pca_scores$sample)
pca_scores <- pca_scores |>
  dplyr::left_join(metadata |> dplyr::select(source_column, condition), by = "source_column")
write_tsv(pca_scores, file.path(src_dir, "mRNA_PCA_source_data.tsv"))

pca_plot <- ggplot2::ggplot(pca_scores, ggplot2::aes(PC1, PC2, color = condition, label = sample)) +
  ggplot2::geom_point(size = 3) +
  ggrepel::geom_text_repel(size = 2.5, show.legend = FALSE) +
  ggplot2::theme_bw(base_size = 10) +
  ggplot2::labs(title = "mRNA-seq PCA", color = "Condition")
ggplot2::ggsave(file.path(fig_dir, "mRNA_PCA.pdf"), pca_plot, width = 5.5, height = 4.5)
ggplot2::ggsave(file.path(fig_dir, "mRNA_PCA.png"), pca_plot, width = 5.5, height = 4.5, dpi = 300)

cor_mat <- stats::cor(log_cpm, method = "pearson")
write_tsv(as.data.frame(as.table(cor_mat)) |>
            dplyr::rename(sample_1 = Var1, sample_2 = Var2, pearson_correlation = Freq),
          file.path(src_dir, "mRNA_sample_correlation_source_data.tsv"))
grDevices::pdf(file.path(fig_dir, "mRNA_sample_correlation.pdf"), width = 6, height = 5.5)
pheatmap::pheatmap(cor_mat, main = "mRNA sample correlation")
grDevices::dev.off()
grDevices::png(file.path(fig_dir, "mRNA_sample_correlation.png"), width = 1800, height = 1600, res = 300)
pheatmap::pheatmap(cor_mat, main = "mRNA sample correlation")
grDevices::dev.off()

key_genes <- dplyr::bind_rows(all_tabs) |>
  dplyr::filter(gene_symbol %in% c("CEMIP", "TNC", "LIFR"))
write_tsv(key_genes, file.path(out_dir, "mRNA_key_gene_validation.tsv"))

message("Final mRNA outputs recomputed from compact count matrix: ", out_dir)
