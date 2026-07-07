project_root <- function() {
  env <- Sys.getenv("IL4_IL13_PROJECT_ROOT", unset = "")
  if (nzchar(env)) return(normalizePath(env, mustWork = FALSE))
  wd <- normalizePath(getwd(), mustWork = TRUE)
  cur <- wd
  repeat {
    if (file.exists(file.path(cur, "config", "analysis.yml"))) return(cur)
    parent <- dirname(cur)
    if (identical(parent, cur)) break
    cur <- parent
  }
  wd
}

ROOT <- project_root()
p <- function(...) file.path(ROOT, ...)

load_pkgs <- function(pkgs, required = TRUE) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) && required) stop("Missing required packages: ", paste(missing, collapse = ", "))
  invisible(lapply(pkgs, function(pkg) suppressPackageStartupMessages(library(pkg, character.only = TRUE))))
}

ensure_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  path
}

write_tsv <- function(x, path) {
  ensure_dir(dirname(path))
  x <- as.data.frame(x, stringsAsFactors = FALSE)
  x[] <- lapply(x, function(col) {
    if (is.list(col)) {
      vapply(col, function(item) paste(as.character(item), collapse = ";"), character(1))
    } else {
      col
    }
  })
  utils::write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}

read_tsv <- function(path) {
  utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
}

copy_if_newer <- function(from, to) {
  ensure_dir(dirname(to))
  invisible(file.copy(from, to, overwrite = TRUE))
}

run_gprofiler <- function(query, background, analysis_id, out_dir, id_type, sources = c("GO:MF", "GO:BP", "GO:CC", "KEGG")) {
  ensure_dir(out_dir)
  query <- unique(query[!is.na(query) & nzchar(query)])
  background <- unique(background[!is.na(background) & nzchar(background)])
  writeLines(query, file.path(out_dir, paste0(analysis_id, "_input_list.txt")))
  writeLines(background, file.path(out_dir, paste0(analysis_id, "_background_list.txt")))
  params <- data.frame(
    parameter = c("analysis_id", "organism", "identifier_type", "sources", "correction_method", "ordered_query", "domain_scope", "query_count", "background_count"),
    value = c(analysis_id, "hsapiens", id_type, paste(sources, collapse = ";"), "g_SCS", "FALSE", "custom", length(query), length(background))
  )
  write_tsv(params, file.path(out_dir, paste0(analysis_id, "_parameters.tsv")))
  if (!requireNamespace("gprofiler2", quietly = TRUE)) {
    writeLines("gprofiler2 is not installed; enrichment was not run.", file.path(out_dir, paste0(analysis_id, "_ERROR.txt")))
    return(data.frame())
  }
  if (length(query) == 0 || length(background) == 0) {
    writeLines("Empty query or background; enrichment was not run.", file.path(out_dir, paste0(analysis_id, "_ERROR.txt")))
    return(data.frame())
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
  if (inherits(res, "error")) {
    writeLines(conditionMessage(res), file.path(out_dir, paste0(analysis_id, "_ERROR.txt")))
    return(data.frame())
  }
  result <- res$result
  if (is.null(result)) result <- data.frame()
  write_tsv(result, file.path(out_dir, paste0(analysis_id, "_complete_results.tsv")))
  sig <- if (nrow(result)) dplyr::filter(result, p_value < 0.05) else result
  write_tsv(sig, file.path(out_dir, paste0(analysis_id, "_significant_results.tsv")))
  plot_terms <- if (nrow(sig)) {
    sig |>
      dplyr::arrange(p_value) |>
      dplyr::group_by(source) |>
      dplyr::slice_head(n = 8) |>
      dplyr::ungroup()
  } else if (nrow(result)) {
    result |>
      dplyr::arrange(p_value) |>
      dplyr::group_by(source) |>
      dplyr::slice_head(n = 8) |>
      dplyr::ungroup() |>
      dplyr::mutate(plot_note = "No terms passed adjusted P < 0.05; top adjusted terms shown for transparency.")
  } else {
    sig
  }
  write_tsv(plot_terms, file.path(out_dir, paste0(analysis_id, "_plotted_terms.tsv")))
  write_tsv(plot_terms, file.path(out_dir, paste0(analysis_id, "_figure_source_data.tsv")))
  if (nrow(plot_terms)) {
    plot_terms$term_label <- factor(plot_terms$term_name, levels = rev(unique(plot_terms$term_name)))
    g <- ggplot2::ggplot(plot_terms, ggplot2::aes(x = -log10(p_value), y = term_label, fill = source)) +
      ggplot2::geom_col() +
      ggplot2::facet_wrap(~source, scales = "free_y") +
      ggplot2::theme_bw(base_size = 9) +
      ggplot2::labs(x = "-log10 adjusted pathway P", y = NULL, title = analysis_id)
    if (!nrow(sig)) {
      g <- g +
        ggplot2::geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "grey30") +
        ggplot2::labs(subtitle = "No terms passed adjusted P < 0.05; top adjusted terms shown")
    }
    ggplot2::ggsave(file.path(out_dir, paste0(analysis_id, "_enrichment.pdf")), g, width = 8, height = 7)
    ggplot2::ggsave(file.path(out_dir, paste0(analysis_id, "_enrichment.png")), g, width = 8, height = 7, dpi = 300)
  }
  result
}

direction_from_logfc <- function(x) {
  dplyr::case_when(is.na(x) ~ NA_character_, x > 0 ~ "up", x < 0 ~ "down", TRUE ~ "unchanged")
}

threshold_log2fc <- log2(1.5)

classify_de <- function(tab, use_log2fc_threshold = TRUE) {
  tab |>
    dplyr::mutate(
      direction = direction_from_logfc(log2_fold_change),
      nominal_discovery = !is.na(raw_p_value) & raw_p_value < 0.05 &
        (!use_log2fc_threshold | (!is.na(log2_fold_change) & abs(log2_fold_change) >= threshold_log2fc)),
      FDR_significant = !is.na(BH_FDR) & BH_FDR < 0.05 &
        (!use_log2fc_threshold | (!is.na(log2_fold_change) & abs(log2_fold_change) >= threshold_log2fc))
    )
}

write_split_de <- function(tab, out_dir, stem) {
  ensure_dir(out_dir)
  write_tsv(tab, file.path(out_dir, paste0(stem, "_complete.tsv")))
  write_tsv(dplyr::filter(tab, nominal_discovery), file.path(out_dir, paste0(stem, "_nominal_discovery.tsv")))
  write_tsv(dplyr::filter(tab, nominal_discovery, direction == "up"), file.path(out_dir, paste0(stem, "_nominal_discovery_up.tsv")))
  write_tsv(dplyr::filter(tab, nominal_discovery, direction == "down"), file.path(out_dir, paste0(stem, "_nominal_discovery_down.tsv")))
  write_tsv(dplyr::filter(tab, FDR_significant), file.path(out_dir, paste0(stem, "_FDR_significant.tsv")))
  write_tsv(dplyr::filter(tab, FDR_significant, direction == "up"), file.path(out_dir, paste0(stem, "_FDR_significant_up.tsv")))
  write_tsv(dplyr::filter(tab, FDR_significant, direction == "down"), file.path(out_dir, paste0(stem, "_FDR_significant_down.tsv")))
}

plot_volcano <- function(tab, title, out_pdf, out_png = sub("[.]pdf$", ".png", out_pdf)) {
  if (!"label_for_plot" %in% names(tab)) tab$label_for_plot <- tab$feature_id
  if (!"gene_symbol" %in% names(tab)) tab$gene_symbol <- NA_character_
  plot_tab <- tab |>
    dplyr::filter(!is.na(raw_p_value), !is.na(log2_fold_change)) |>
    dplyr::mutate(
      neg_log10_p = -log10(pmax(raw_p_value, .Machine$double.xmin)),
      status = dplyr::case_when(
        FDR_significant ~ "FDR_significant",
        nominal_discovery ~ "nominal_discovery",
        TRUE ~ "not_threshold"
      ),
      label_for_plot = dplyr::coalesce(label_for_plot, gene_symbol, feature_id)
    )
  labels <- plot_tab |>
    dplyr::filter(status != "not_threshold") |>
    dplyr::arrange(raw_p_value) |>
    dplyr::slice_head(n = 12)
  g <- ggplot2::ggplot(plot_tab, ggplot2::aes(log2_fold_change, neg_log10_p, color = status)) +
    ggplot2::geom_point(alpha = 0.65, size = 1.2) +
    ggplot2::geom_vline(xintercept = c(-threshold_log2fc, threshold_log2fc), linetype = "dashed", color = "grey40") +
    ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
    ggrepel::geom_text_repel(data = labels, ggplot2::aes(label = label_for_plot), size = 2.4, max.overlaps = Inf, show.legend = FALSE) +
    ggplot2::scale_color_manual(values = c(FDR_significant = "#1b7837", nominal_discovery = "#2166ac", not_threshold = "grey70")) +
    ggplot2::labs(title = title, x = "log2 fold change (Treatment - Vehicle)", y = "-log10 nominal P", color = "Status") +
    ggplot2::theme_bw(base_size = 10)
  ensure_dir(dirname(out_pdf))
  ggplot2::ggsave(out_pdf, g, width = 6.5, height = 5.2)
  ggplot2::ggsave(out_png, g, width = 6.5, height = 5.2, dpi = 300)
  invisible(g)
}

write_session_info <- function(path = p("sessionInfo.txt")) {
  ensure_dir(dirname(path))
  writeLines(capture.output(sessionInfo()), path)
}
