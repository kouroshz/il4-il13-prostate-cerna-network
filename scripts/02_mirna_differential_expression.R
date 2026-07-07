#!/usr/bin/env Rscript

# Publication-facing miRNA-seq differential-expression workflow.
# Uses log2(CPM + 1) values with two-sided Welch tests.

source(file.path(getwd(), "R", "final_helpers.R"))
load_pkgs(c("dplyr", "ggplot2", "ggrepel"))

table_dir <- ensure_dir(p("results", "tables"))
figure_dir <- ensure_dir(p("results", "figures"))

priority_mirnas <- c("hsa-miR-140-3p", "hsa-miR-135b-5p", "hsa-miR-625-3p")

plot_mirna_volcano <- function(tab, title, out_pdf, out_png = sub("[.]pdf$", ".png", out_pdf)) {
  plot_tab <- tab |>
    dplyr::filter(!is.na(raw_p_value), !is.na(log2_fold_change)) |>
    dplyr::mutate(
      neg_log10_p = -log10(pmax(raw_p_value, .Machine$double.xmin)),
      plot_category = dplyr::case_when(
        biologically_prioritized ~ "Prioritized",
        nominal_discovery ~ "Nominal P < 0.05",
        TRUE ~ "Not nominal"
      )
    )
  labels <- dplyr::bind_rows(
    plot_tab |> dplyr::filter(biologically_prioritized),
    plot_tab |>
      dplyr::filter(nominal_discovery, !biologically_prioritized) |>
      dplyr::arrange(raw_p_value) |>
      dplyr::slice_head(n = 6)
  ) |>
    dplyr::distinct(feature_id, .keep_all = TRUE)

  g <- ggplot2::ggplot(plot_tab, ggplot2::aes(log2_fold_change, neg_log10_p, color = plot_category)) +
    ggplot2::geom_point(alpha = 0.65, size = 1.2) +
    ggplot2::geom_point(
      data = plot_tab |> dplyr::filter(biologically_prioritized),
      ggplot2::aes(log2_fold_change, neg_log10_p),
      shape = 21,
      fill = "white",
      color = "#a33f3f",
      stroke = 0.9,
      size = 2.4,
      inherit.aes = FALSE
    ) +
    ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
    ggrepel::geom_text_repel(
      data = labels,
      ggplot2::aes(label = label_for_plot),
      size = 2.4,
      box.padding = 0.35,
      point.padding = 0.25,
      max.overlaps = Inf,
      show.legend = FALSE
    ) +
    ggplot2::scale_color_manual(
      values = c("Prioritized" = "#a33f3f", "Nominal P < 0.05" = "#416f9f", "Not nominal" = "grey70"),
      breaks = c("Prioritized", "Nominal P < 0.05", "Not nominal")
    ) +
    ggplot2::labs(title = title, x = "log2 fold change (Treatment \u2212 Vehicle)", y = "\u2212log10 nominal P", color = "Status") +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(plot.margin = ggplot2::margin(8, 18, 8, 8))
  ggplot2::ggsave(out_pdf, g, width = 6.5, height = 5.2, device = grDevices::cairo_pdf)
  ggplot2::ggsave(out_png, g, width = 6.5, height = 5.2, dpi = 300)
  invisible(g)
}

safe_ttest <- function(x, y) {
  tryCatch({
    if (length(stats::na.omit(x)) < 2 || length(stats::na.omit(y)) < 2) return(NA_real_)
    stats::t.test(x, y, alternative = "two.sided", var.equal = FALSE)$p.value
  }, error = function(e) NA_real_)
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

priority_axis <- function(feature_id) {
  dplyr::case_when(
    feature_id == "hsa-miR-140-3p" ~ "miR-140-3p / CEMIP / circPAPPA",
    feature_id == "hsa-miR-135b-5p" ~ "miR-135b-5p / TNC / circLIFR",
    feature_id == "hsa-miR-625-3p" ~ "miR-625-3p / LIFR",
    TRUE ~ ""
  )
}

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
      analysis_method = "log2_CPM_plus_1_two_sided_Welch_test",
      contrast = contrast_name,
      label_for_plot = feature_id
    )
  replicate_cpm <- cpm[match(out$feature_id, cpm$feature_id), ordered_cols, drop = FALSE]
  out <- out |>
    dplyr::bind_cols(replicate_cpm) |>
    dplyr::mutate(dplyr::across(dplyr::all_of(ordered_cols), as.numeric)) |>
    classify_de(use_log2fc_threshold = FALSE) |>
    dplyr::mutate(
      biologically_prioritized = feature_id %in% priority_mirnas,
      experimentally_evaluated = feature_id %in% priority_mirnas,
      prioritized_axis = priority_axis(feature_id),
      interpretation_class = dplyr::case_when(
        biologically_prioritized & nominal_discovery ~ "prioritized axis and nominal P < 0.05",
        biologically_prioritized ~ "prioritized axis; not nominal in this contrast",
        nominal_discovery ~ "nominal P < 0.05",
        TRUE ~ "not nominal"
      )
    ) |>
    dplyr::select(
      assay, analysis_method, contrast, feature_id, log2_fold_change, raw_p_value, BH_FDR,
      direction, nominal_discovery, FDR_significant,
      biologically_prioritized, experimentally_evaluated, prioritized_axis, interpretation_class,
      mean_CPM_vehicle, mean_CPM_treatment,
      mean_log2_CPM_plus_1_vehicle, mean_log2_CPM_plus_1_treatment,
      dplyr::all_of(ordered_cols), label_for_plot
    )

  prefix <- if (treatment == "IL4") "Figure_1B_miRNA_IL4_volcano" else "Figure_1C_miRNA_IL13_volcano"
  plot_mirna_volcano(out, paste("miRNA", treatment, "vs Vehicle"), file.path(figure_dir, paste0(prefix, ".pdf")))
  out
}

all_tabs <- list(
  IL4_vs_Vehicle = run_contrast("IL4"),
  IL13_vs_Vehicle = run_contrast("IL13")
)

summary_counts <- dplyr::bind_rows(lapply(names(all_tabs), function(contrast) {
  tab <- all_tabs[[contrast]]
  dplyr::tibble(
    assay = "miRNA-seq",
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

expected <- c(IL4_vs_Vehicle = 45, IL13_vs_Vehicle = 28, shared_same_direction = 10)
observed <- c(
  IL4_vs_Vehicle = summary_counts$nominal_discovery[summary_counts$contrast == "IL4_vs_Vehicle"],
  IL13_vs_Vehicle = summary_counts$nominal_discovery[summary_counts$contrast == "IL13_vs_Vehicle"],
  shared_same_direction = nrow(shared_nominal)
)
if (!identical(as.integer(observed), as.integer(expected))) {
  stop("miRNA count mismatch: observed ", paste(names(observed), observed, sep = "=", collapse = "; "))
}

write_tsv(
  dplyr::bind_rows(all_tabs) |>
    dplyr::arrange(contrast, raw_p_value, feature_id),
  file.path(table_dir, "Table_S1_miRNA_differential_expression.tsv")
)

save_figure(
  plot_two_set_venn(
    title = "miRNA differential-expression summary",
    left_label = "IL-4 nominal P < 0.05",
    right_label = "IL-13 nominal P < 0.05",
    left_only = as.integer(observed["IL4_vs_Vehicle"] - observed["shared_same_direction"]),
    shared = as.integer(observed["shared_same_direction"]),
    right_only = as.integer(observed["IL13_vs_Vehicle"] - observed["shared_same_direction"]),
    shared_note = "Shared = same direction"
  ),
  "Figure_1A_miRNA_overlap",
  width = 6.2,
  height = 4.6
)

message("Table S1 and Figure 1A-C outputs written.")
