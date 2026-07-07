#!/usr/bin/env Rscript

# Final miRNA-seq differential-expression workflow.
# Recomputes Welch-test statistics from the verified compact CPM matrix.
# The manuscript-aligned exploratory miRNA screen uses nominal P < 0.05
# without forcing the mRNA/circRNA 1.5-fold cutoff. Experimentally evaluated
# candidate miRNAs are highlighted explicitly as a separate plot category.

source(file.path(getwd(), "R", "final_helpers.R"))
load_pkgs(c("dplyr", "ggplot2", "ggrepel", "pheatmap"))

out_dir <- ensure_dir(p("data", "processed", "final", "mirna"))
fig_dir <- ensure_dir(p("results", "final", "figures", "mirna"))
src_dir <- ensure_dir(p("results", "final", "source_data", "mirna"))

candidate_mirnas <- c("hsa-miR-140-3p", "hsa-miR-135b-5p", "hsa-miR-625-3p")

plot_mirna_volcano <- function(tab, title, out_pdf, out_png = sub("[.]pdf$", ".png", out_pdf)) {
  plot_tab <- tab |>
    dplyr::filter(!is.na(raw_p_value), !is.na(log2_fold_change)) |>
    dplyr::mutate(
      neg_log10_p = -log10(pmax(raw_p_value, .Machine$double.xmin)),
      plot_category = dplyr::case_when(
        biologically_prioritized ~ "biologically_prioritized_tested",
        FDR_significant ~ "FDR_significant",
        nominal_discovery ~ "nominal_P_lt_0.05",
        TRUE ~ "not_nominal"
      )
    )
  labels <- dplyr::bind_rows(
    plot_tab |> dplyr::filter(biologically_prioritized),
    plot_tab |>
      dplyr::filter(nominal_discovery, !biologically_prioritized) |>
      dplyr::arrange(raw_p_value) |>
      dplyr::slice_head(n = 9)
  ) |>
    dplyr::distinct(feature_id, .keep_all = TRUE)

  g <- ggplot2::ggplot(plot_tab, ggplot2::aes(log2_fold_change, neg_log10_p, color = plot_category)) +
    ggplot2::geom_point(alpha = 0.65, size = 1.2) +
    ggplot2::geom_point(
      data = plot_tab |> dplyr::filter(biologically_prioritized),
      ggplot2::aes(log2_fold_change, neg_log10_p),
      shape = 21,
      fill = "white",
      color = "#b2182b",
      stroke = 0.9,
      size = 2.4,
      inherit.aes = FALSE
    ) +
    ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
    ggrepel::geom_text_repel(
      data = labels,
      ggplot2::aes(label = label_for_plot),
      size = 2.4,
      max.overlaps = Inf,
      show.legend = FALSE
    ) +
    ggplot2::scale_color_manual(
      values = c(
        biologically_prioritized_tested = "#b2182b",
        FDR_significant = "#1b7837",
        nominal_P_lt_0.05 = "#2166ac",
        not_nominal = "grey70"
      ),
      breaks = c("biologically_prioritized_tested", "FDR_significant", "nominal_P_lt_0.05", "not_nominal")
    ) +
    ggplot2::labs(
      title = title,
      x = "log2 fold change (Treatment - Vehicle)",
      y = "-log10 nominal P",
      color = "Status"
    ) +
    ggplot2::theme_bw(base_size = 10)
  ensure_dir(dirname(out_pdf))
  ggplot2::ggsave(out_pdf, g, width = 6.5, height = 5.2)
  ggplot2::ggsave(out_png, g, width = 6.5, height = 5.2, dpi = 300)
  invisible(g)
}

safe_ttest <- function(x, y) {
  tryCatch({
    if (length(stats::na.omit(x)) < 2 || length(stats::na.omit(y)) < 2) return(NA_real_)
    stats::t.test(x, y, alternative = "two.sided", var.equal = FALSE)$p.value
  }, error = function(e) NA_real_)
}

summarise_counts <- function(tab, contrast) {
  dplyr::tibble(
    assay = "miRNA-seq",
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

metadata <- read_tsv(p("data", "metadata", "mirna_samples.tsv"))
required_meta <- c("sample_name", "source_column", "condition", "biological_replicate")
missing_meta <- setdiff(required_meta, names(metadata))
if (length(missing_meta)) stop("Missing miRNA metadata columns: ", paste(missing_meta, collapse = ", "))

expr <- utils::read.csv(p("data", "processed", "mirna", "Expression_Browser_CPM.csv"), check.names = FALSE)
if (!"Name" %in% names(expr)) stop("miRNA CPM matrix must contain Name")
if (!all(metadata$source_column %in% names(expr))) {
  stop("miRNA metadata references missing CPM columns: ", paste(setdiff(metadata$source_column, names(expr)), collapse = ", "))
}

cpm <- expr[, c("Name", metadata$source_column), drop = FALSE]
names(cpm) <- c("feature_id", metadata$sample_name)
cpm[, metadata$sample_name] <- lapply(cpm[, metadata$sample_name, drop = FALSE], as.numeric)
cpm <- cpm[rowSums(cpm[, metadata$sample_name, drop = FALSE], na.rm = TRUE) != 0, , drop = FALSE]

run_contrast <- function(treatment) {
  contrast_name <- paste0(treatment, "_vs_Vehicle")
  meta <- metadata |>
    dplyr::filter(condition %in% c("Vehicle", treatment)) |>
    dplyr::mutate(condition = factor(condition, levels = c("Vehicle", treatment))) |>
    dplyr::arrange(condition, biological_replicate)
  vehicle_cols <- meta$sample_name[meta$condition == "Vehicle"]
  treatment_cols <- meta$sample_name[meta$condition == treatment]
  if (length(vehicle_cols) != 3 || length(treatment_cols) != 3) {
    stop("Expected three Vehicle and three ", treatment, " biological replicates for miRNA")
  }

  rows <- lapply(seq_len(nrow(cpm)), function(i) {
    vehicle_cpm <- as.numeric(cpm[i, vehicle_cols])
    treatment_cpm <- as.numeric(cpm[i, treatment_cols])
    vehicle_log <- log2(vehicle_cpm + 1)
    treatment_log <- log2(treatment_cpm + 1)
    dplyr::tibble(
      feature_id = cpm$feature_id[i],
      mean_CPM_vehicle = mean(vehicle_cpm, na.rm = TRUE),
      mean_CPM_treatment = mean(treatment_cpm, na.rm = TRUE),
      mean_log2_CPM_plus_1_vehicle = mean(vehicle_log, na.rm = TRUE),
      mean_log2_CPM_plus_1_treatment = mean(treatment_log, na.rm = TRUE),
      log2_fold_change = mean(treatment_log, na.rm = TRUE) - mean(vehicle_log, na.rm = TRUE),
      raw_p_value = safe_ttest(treatment_log, vehicle_log)
    )
  })

  ordered_cols <- c(vehicle_cols, treatment_cols)
  out <- dplyr::bind_rows(rows) |>
    dplyr::mutate(
      BH_FDR = stats::p.adjust(raw_p_value, method = "BH"),
      assay = "miRNA-seq",
      analysis_source = "log2_CPM_plus_1_Welch_two_sided_from_compact_CPM_matrix",
      contrast = contrast_name,
      label_for_plot = feature_id
    )
  replicate_cpm <- cpm[match(out$feature_id, cpm$feature_id), ordered_cols, drop = FALSE]
  out <- out |>
    dplyr::bind_cols(replicate_cpm) |>
    dplyr::mutate(dplyr::across(dplyr::all_of(ordered_cols), as.numeric)) |>
    classify_de(use_log2fc_threshold = FALSE) |>
    dplyr::mutate(
      biologically_prioritized = feature_id %in% candidate_mirnas,
      experimentally_evaluated = feature_id %in% candidate_mirnas,
      candidate_highlight = biologically_prioritized,
      candidate_highlight_reason = dplyr::if_else(
        biologically_prioritized,
        "prioritized/tested from integrated ceRNA network; highlighted separately from nominal volcano coloring",
        ""
      ),
      threshold_interpretation = dplyr::case_when(
        biologically_prioritized & nominal_discovery ~ "biologically prioritized/tested and nominal P < 0.05",
        biologically_prioritized ~ "biologically prioritized/tested; not a nominal discovery in this contrast",
        nominal_discovery ~ "nominal P < 0.05 exploratory miRNA screen",
        TRUE ~ "not nominally significant"
      )
    ) |>
    dplyr::select(
      assay, analysis_source, contrast, feature_id, log2_fold_change, raw_p_value, BH_FDR,
      direction, nominal_discovery, FDR_significant,
      biologically_prioritized, experimentally_evaluated, candidate_highlight,
      candidate_highlight_reason, threshold_interpretation,
      mean_CPM_vehicle, mean_CPM_treatment,
      mean_log2_CPM_plus_1_vehicle, mean_log2_CPM_plus_1_treatment,
      dplyr::all_of(ordered_cols), label_for_plot
    )

  stem <- paste0("miRNA_", contrast_name)
  write_split_de(out, out_dir, stem)
  write_tsv(
    out |>
      dplyr::transmute(
        feature_id,
        display_label = label_for_plot,
        log2_fold_change,
        raw_p_value,
        BH_FDR,
        nominal_discovery,
        FDR_significant,
        direction,
        biologically_prioritized,
        experimentally_evaluated,
        candidate_highlight,
        candidate_highlight_reason,
        threshold_interpretation
      ),
    file.path(src_dir, paste0(stem, "_volcano_source_data.tsv"))
  )
  plot_mirna_volcano(out, paste("miRNA", contrast_name), file.path(fig_dir, paste0(stem, "_volcano.pdf")))
  out
}

all_tabs <- list(
  IL4_vs_Vehicle = run_contrast("IL4"),
  IL13_vs_Vehicle = run_contrast("IL13")
)

summary_counts <- dplyr::bind_rows(lapply(names(all_tabs), function(contrast) summarise_counts(all_tabs[[contrast]], contrast)))
write_tsv(summary_counts, file.path(out_dir, "miRNA_summary_counts.tsv"))

make_shared <- function(flag_col, label) {
  x <- dplyr::inner_join(
    dplyr::filter(all_tabs$IL4_vs_Vehicle, .data[[flag_col]]) |>
      dplyr::select(feature_id, direction_IL4 = direction, log2FC_IL4 = log2_fold_change, P_IL4 = raw_p_value, FDR_IL4 = BH_FDR),
    dplyr::filter(all_tabs$IL13_vs_Vehicle, .data[[flag_col]]) |>
      dplyr::select(feature_id, direction_IL13 = direction, log2FC_IL13 = log2_fold_change, P_IL13 = raw_p_value, FDR_IL13 = BH_FDR),
    by = "feature_id"
  ) |>
    dplyr::filter(direction_IL4 == direction_IL13) |>
    dplyr::mutate(shared_direction = direction_IL4)
  write_tsv(x, file.path(out_dir, paste0("miRNA_shared_direction_", label, ".tsv")))
}
make_shared("nominal_discovery", "nominal_discovery")
make_shared("FDR_significant", "FDR_significant")

all_mirnas <- dplyr::full_join(
  all_tabs$IL4_vs_Vehicle |>
    dplyr::select(feature_id, IL4_direction = direction, IL4_nominal = nominal_discovery, IL4_FDR_significant = FDR_significant, IL4_log2FC = log2_fold_change, IL4_P = raw_p_value, IL4_FDR = BH_FDR),
  all_tabs$IL13_vs_Vehicle |>
    dplyr::select(feature_id, IL13_direction = direction, IL13_nominal = nominal_discovery, IL13_FDR_significant = FDR_significant, IL13_log2FC = log2_fold_change, IL13_P = raw_p_value, IL13_FDR = BH_FDR),
  by = "feature_id"
)
write_tsv(all_mirnas, file.path(out_dir, "miRNA_exact_overlap_table.tsv"))

candidate_values <- dplyr::bind_rows(all_tabs) |>
  dplyr::filter(feature_id %in% candidate_mirnas) |>
  dplyr::mutate(
    candidate_axis = dplyr::case_when(
      feature_id == "hsa-miR-140-3p" ~ "miR-140-3p / CEMIP / circPAPPA",
      feature_id == "hsa-miR-135b-5p" ~ "miR-135b-5p / TNC / circLIFR",
      feature_id == "hsa-miR-625-3p" ~ "miR-625-3p / LIFR",
      TRUE ~ ""
    ),
    experimental_status = dplyr::case_when(
      feature_id == "hsa-miR-140-3p" ~ "strongest functional support; pursued further",
      TRUE ~ "experimentally evaluated; negative or less supportive result"
    )
  )
write_tsv(candidate_values, file.path(out_dir, "candidate_miRNA_final_values.tsv"))
write_tsv(
  candidate_values |> dplyr::filter(feature_id == "hsa-miR-140-3p"),
  file.path(out_dir, "miR_140_3p_final_values.tsv")
)

message("Final miRNA outputs recomputed from compact CPM matrix: ", out_dir)
