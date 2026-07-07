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

clean_output_dir <- function(path) {
  ensure_dir(path)
  entries <- list.files(path, all.files = TRUE, no.. = TRUE, full.names = TRUE)
  entries <- entries[basename(entries) != ".gitkeep"]
  if (length(entries)) unlink(entries, recursive = TRUE, force = TRUE)
  invisible(path)
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

save_figure <- function(plot, stem, width = 6.5, height = 5.0, dpi = 300) {
  figure_dir <- ensure_dir(p("results", "figures"))
  ggplot2::ggsave(file.path(figure_dir, paste0(stem, ".pdf")), plot, width = width, height = height)
  ggplot2::ggsave(file.path(figure_dir, paste0(stem, ".png")), plot, width = width, height = height, dpi = dpi)
}

plot_overlap_counts <- function(counts, title, y_label = "Features") {
  counts$category <- factor(counts$category, levels = counts$category)
  ggplot2::ggplot(counts, ggplot2::aes(category, count)) +
    ggplot2::geom_col(fill = "#4f7f9f", width = 0.62) +
    ggplot2::geom_text(ggplot2::aes(label = count), vjust = -0.35, size = 4) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.12))) +
    ggplot2::labs(title = title, x = NULL, y = y_label) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(size = 10)
    )
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
