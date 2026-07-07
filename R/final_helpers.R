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
  ggplot2::ggsave(file.path(figure_dir, paste0(stem, ".pdf")), plot, width = width, height = height, device = grDevices::cairo_pdf, bg = "white")
  ggplot2::ggsave(file.path(figure_dir, paste0(stem, ".png")), plot, width = width, height = height, dpi = dpi, bg = "white")
}

plot_two_set_venn <- function(
    title,
    left_label,
    right_label,
    left_only,
    shared,
    right_only,
    shared_note = NULL,
    bottom_note = NULL) {
  circle_points <- function(center_x, set_name) {
    theta <- seq(0, 2 * pi, length.out = 240)
    data.frame(
      x = center_x + 1.35 * cos(theta),
      y = 1.35 * sin(theta),
      set_name = set_name
    )
  }
  circles <- rbind(circle_points(-0.75, "left"), circle_points(0.75, "right"))
  annotations <- data.frame(
    x = c(-1.25, 0, 1.25),
    y = c(0, 0, 0),
    label = c(left_only, shared, right_only)
  )
  g <- ggplot2::ggplot() +
    ggplot2::geom_polygon(
      data = circles,
      ggplot2::aes(x, y, group = set_name, fill = set_name),
      alpha = 0.42,
      color = "grey25",
      linewidth = 0.7
    ) +
    ggplot2::geom_text(data = annotations, ggplot2::aes(x, y, label = label), size = 7, fontface = "bold") +
    ggplot2::annotate("text", x = -1.2, y = 1.48, label = left_label, size = 3.5, fontface = "bold") +
    ggplot2::annotate("text", x = 1.2, y = 1.48, label = right_label, size = 3.5, fontface = "bold") +
    ggplot2::scale_fill_manual(values = c(left = "#78a6d8", right = "#8fcb8f")) +
    ggplot2::coord_equal(xlim = c(-2.5, 2.5), ylim = c(-2.05, 1.75), clip = "off") +
    ggplot2::labs(title = title) +
    ggplot2::theme_void(base_size = 11) +
    ggplot2::theme(
      legend.position = "none",
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 13),
      plot.margin = ggplot2::margin(10, 18, 28, 18)
    )
  if (!is.null(shared_note) && nzchar(shared_note)) {
    g <- g + ggplot2::annotate("text", x = 0, y = -1.52, label = shared_note, size = 3.1, color = "grey25")
  }
  if (!is.null(bottom_note) && nzchar(bottom_note)) {
    note_y <- if (!is.null(shared_note) && nzchar(shared_note)) -1.80 else -1.58
    g <- g + ggplot2::annotate("text", x = 0, y = note_y, label = bottom_note, size = 3.0, color = "grey30")
  }
  g
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

plot_volcano <- function(
    tab,
    title,
    out_pdf,
    out_png = sub("[.]pdf$", ".png", out_pdf),
    y_cap = Inf,
    subtitle = NULL,
    label_n = 12,
    label_features = NULL) {
  if (!"label_for_plot" %in% names(tab)) tab$label_for_plot <- tab$feature_id
  if (!"gene_symbol" %in% names(tab)) tab$gene_symbol <- NA_character_
  plot_tab <- tab |>
    dplyr::filter(!is.na(raw_p_value), !is.na(log2_fold_change)) |>
    dplyr::mutate(
      neg_log10_p = -log10(pmax(raw_p_value, .Machine$double.xmin)),
      plotted_neg_log10_p = pmin(neg_log10_p, y_cap),
      status = dplyr::case_when(
        FDR_significant ~ "FDR < 0.05",
        nominal_discovery ~ "Nominal P < 0.05",
        TRUE ~ "Not nominal"
      ),
      label_for_plot = dplyr::coalesce(label_for_plot, gene_symbol, feature_id)
    )
  if (is.null(label_features)) {
    labels <- plot_tab |>
      dplyr::filter(status != "Not nominal", !is.na(label_for_plot), nzchar(label_for_plot)) |>
      dplyr::arrange(raw_p_value) |>
      dplyr::distinct(label_for_plot, .keep_all = TRUE) |>
      dplyr::slice_head(n = label_n)
  } else {
    labels <- plot_tab |>
      dplyr::filter(
        status != "Not nominal",
        !is.na(label_for_plot),
        nzchar(label_for_plot),
        label_for_plot %in% label_features | gene_symbol %in% label_features | feature_id %in% label_features
      ) |>
      dplyr::arrange(factor(label_for_plot, levels = label_features), raw_p_value) |>
      dplyr::distinct(label_for_plot, .keep_all = TRUE)
  }
  y_label <- if (is.finite(y_cap)) paste0("\u2212log10 nominal P, capped at ", y_cap) else "\u2212log10 nominal P"
  g <- ggplot2::ggplot(plot_tab, ggplot2::aes(log2_fold_change, plotted_neg_log10_p, color = status)) +
    ggplot2::geom_point(alpha = 0.65, size = 1.2) +
    ggplot2::geom_vline(xintercept = c(-threshold_log2fc, threshold_log2fc), linetype = "dashed", color = "grey40") +
    ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
    ggrepel::geom_text_repel(data = labels, ggplot2::aes(label = label_for_plot), size = 2.4, max.overlaps = Inf, show.legend = FALSE) +
    ggplot2::scale_color_manual(
      values = c("FDR < 0.05" = "#1b7837", "Nominal P < 0.05" = "#2166ac", "Not nominal" = "grey70"),
      breaks = c("FDR < 0.05", "Nominal P < 0.05", "Not nominal")
    ) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = "log2 fold change (Treatment \u2212 Vehicle)",
      y = y_label,
      color = "Status"
    ) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme(plot.margin = ggplot2::margin(8, 18, 8, 8))
  if (is.finite(y_cap)) {
    g <- g + ggplot2::scale_y_continuous(limits = c(0, y_cap), expand = ggplot2::expansion(mult = c(0, 0.03)))
  }
  ensure_dir(dirname(out_pdf))
  ggplot2::ggsave(out_pdf, g, width = 6.5, height = 5.2, device = grDevices::cairo_pdf, bg = "white")
  ggplot2::ggsave(out_png, g, width = 6.5, height = 5.2, dpi = 300, bg = "white")
  invisible(g)
}

write_session_info <- function(path = p("sessionInfo.txt")) {
  ensure_dir(dirname(path))
  writeLines(capture.output(sessionInfo()), path)
}
