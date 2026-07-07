#!/usr/bin/env Rscript

# Publication-facing enrichment workflow for mRNA and circRNA results.

source(file.path(getwd(), "R", "final_helpers.R"))
load_pkgs(c("dplyr", "ggplot2"), required = TRUE)

table_dir <- ensure_dir(p("results", "tables"))
figure_dir <- ensure_dir(p("results", "figures"))

sources <- c("GO:MF", "GO:BP", "GO:CC", "KEGG")

empty_enrichment <- function(assay, contrast, analysis_id, query_count, background_count, note) {
  dplyr::tibble(
    assay = assay,
    contrast = contrast,
    analysis_id = analysis_id,
    source = NA_character_,
    term_id = NA_character_,
    term_name = note,
    p_value = NA_real_,
    query_size = query_count,
    intersection_size = NA_integer_,
    term_size = NA_integer_,
    effective_domain_size = background_count,
    intersection = NA_character_
  )
}

run_enrichment <- function(query, background, assay, contrast, analysis_id) {
  query <- unique(query[!is.na(query) & nzchar(query)])
  background <- unique(background[!is.na(background) & nzchar(background)])
  if (!requireNamespace("gprofiler2", quietly = TRUE)) {
    return(empty_enrichment(assay, contrast, analysis_id, length(query), length(background), "gprofiler2 package not installed"))
  }
  if (!length(query) || !length(background)) {
    return(empty_enrichment(assay, contrast, analysis_id, length(query), length(background), "empty enrichment query or background"))
  }
  res <- tryCatch(
    gprofiler2::gost(
      query = query,
      organism = "hsapiens",
      ordered_query = FALSE,
      multi_query = FALSE,
      significant = FALSE,
      exclude_iea = FALSE,
      measure_underrepresentation = FALSE,
      evcodes = FALSE,
      user_threshold = 0.05,
      correction_method = "g_SCS",
      domain_scope = "custom",
      custom_bg = background,
      sources = sources
    ),
    error = function(e) e
  )
  if (inherits(res, "error") || is.null(res$result) || !nrow(res$result)) {
    note <- if (inherits(res, "error")) conditionMessage(res) else "no enrichment results returned"
    return(empty_enrichment(assay, contrast, analysis_id, length(query), length(background), note))
  }
  result <- res$result
  optional_cols <- c("intersection", "intersection_size", "term_size", "effective_domain_size")
  for (col in setdiff(optional_cols, names(result))) result[[col]] <- NA
  result |>
    dplyr::mutate(
      assay = assay,
      contrast = contrast,
      analysis_id = analysis_id,
      query_size = length(query),
      effective_domain_size = length(background)
    ) |>
    dplyr::select(
      assay, contrast, analysis_id, source, term_id, term_name, p_value,
      query_size, intersection_size, term_size, effective_domain_size, intersection
    )
}

plot_enrichment <- function(tab, title, stem) {
  plot_tab <- tab |>
    dplyr::filter(!is.na(p_value), !is.na(term_name), nzchar(term_name)) |>
    dplyr::arrange(p_value) |>
    dplyr::group_by(source) |>
    dplyr::slice_head(n = 6) |>
    dplyr::ungroup()
  if (!nrow(plot_tab)) {
    plot_tab <- dplyr::tibble(term_name = "No enrichment terms available", source = "Result", p_value = 1)
  }
  plot_tab$term_label <- factor(plot_tab$term_name, levels = rev(unique(plot_tab$term_name)))
  g <- ggplot2::ggplot(plot_tab, ggplot2::aes(x = -log10(pmax(p_value, .Machine$double.xmin)), y = term_label, fill = source)) +
    ggplot2::geom_col(width = 0.68) +
    ggplot2::geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "grey35") +
    ggplot2::facet_wrap(~source, scales = "free_y") +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank(), legend.position = "none") +
    ggplot2::labs(title = title, x = "-log10 adjusted pathway P", y = NULL)
  save_figure(g, stem, width = 8.2, height = 6.4)
}

mrna <- read_tsv(p("results", "tables", "Table_S3_mRNA_differential_expression.tsv"))
circrna <- read_tsv(p("results", "tables", "Table_S2_circRNA_differential_expression.tsv"))

analyses <- list(
  list(
    assay = "mRNA-seq",
    contrast = "IL4_vs_Vehicle",
    analysis_id = "mRNA_IL4_up_nominal",
    query = mrna |> dplyr::filter(contrast == "IL4_vs_Vehicle", nominal_discovery, direction == "up") |> dplyr::pull(feature_id),
    background = mrna |> dplyr::filter(contrast == "IL4_vs_Vehicle") |> dplyr::pull(feature_id),
    figure = "Figure_3D_mRNA_IL4_enrichment",
    title = "mRNA IL-4 upregulated gene enrichment"
  ),
  list(
    assay = "mRNA-seq",
    contrast = "IL13_vs_Vehicle",
    analysis_id = "mRNA_IL13_up_nominal",
    query = mrna |> dplyr::filter(contrast == "IL13_vs_Vehicle", nominal_discovery, direction == "up") |> dplyr::pull(feature_id),
    background = mrna |> dplyr::filter(contrast == "IL13_vs_Vehicle") |> dplyr::pull(feature_id),
    figure = "Figure_3E_mRNA_IL13_enrichment",
    title = "mRNA IL-13 upregulated gene enrichment"
  ),
  list(
    assay = "circRNA microarray",
    contrast = "IL4_vs_Vehicle",
    analysis_id = "circRNA_IL4_up_nominal_host_genes",
    query = circrna |> dplyr::filter(contrast == "IL4_vs_Vehicle", nominal_discovery, direction == "up") |> dplyr::pull(gene_symbol),
    background = circrna |> dplyr::filter(contrast == "IL4_vs_Vehicle") |> dplyr::pull(gene_symbol)
  ),
  list(
    assay = "circRNA microarray",
    contrast = "IL13_vs_Vehicle",
    analysis_id = "circRNA_IL13_up_nominal_host_genes",
    query = circrna |> dplyr::filter(contrast == "IL13_vs_Vehicle", nominal_discovery, direction == "up") |> dplyr::pull(gene_symbol),
    background = circrna |> dplyr::filter(contrast == "IL13_vs_Vehicle") |> dplyr::pull(gene_symbol)
  )
)

results <- lapply(analyses, function(x) {
  tab <- run_enrichment(x$query, x$background, x$assay, x$contrast, x$analysis_id)
  if (!is.null(x$figure)) plot_enrichment(tab, x$title, x$figure)
  tab
})

write_tsv(
  dplyr::bind_rows(results) |>
    dplyr::arrange(assay, contrast, p_value, term_name),
  file.path(table_dir, "Table_S4_enrichment_results.tsv")
)

message("Table S4 and Figure 3D-E outputs written.")
