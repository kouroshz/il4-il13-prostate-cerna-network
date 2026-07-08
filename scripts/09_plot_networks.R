#!/usr/bin/env Rscript

# Publication-facing ceRNA network figures.

source(file.path(getwd(), "R", "final_helpers.R"))
load_pkgs(c("dplyr", "igraph", "ggplot2"))

figure_dir <- ensure_dir(p("results", "figures"))

parse_miranda_file <- function(path, comparison) {
  x <- utils::read.table(path, sep = "\t", header = FALSE, fill = TRUE, quote = "", comment.char = "")
  x <- x[grepl("^>", x$V1) & !grepl("^>>", x$V1), , drop = FALSE]
  if (!nrow(x)) return(data.frame())
  x$V1 <- sub("^>", "", x$V1)
  dplyr::tibble(
    comparison = comparison,
    miRNA = as.character(x$V1),
    target_raw = as.character(x$V2),
    target_primary = sub("[|:].*", "", target_raw),
    target_symbol = ifelse(grepl("[|:]", target_raw), sub(".*[|:]", "", target_raw), target_raw),
    score = as.numeric(x$V3),
    energy = as.numeric(x$V4)
  )
}

miranda_files <- c(
  IL4_down = "IL4_vs_V_miRNA_down_against_3utr_circ_filtered_targets.txt",
  IL13_down = "IL13_vs_V_miRNA_down_against_3utr_circ_filtered_targets.txt"
)
miranda_pairs <- dplyr::bind_rows(lapply(names(miranda_files), function(nm) {
  parse_miranda_file(p("data", "processed", "miranda_outputs", miranda_files[[nm]]), nm)
})) |>
  dplyr::group_by(comparison, miRNA, target_raw, target_primary, target_symbol) |>
  dplyr::summarise(
    miRanda_score = max(score, na.rm = TRUE),
    hybridization_energy = min(energy, na.rm = TRUE),
    site_count = dplyr::n(),
    .groups = "drop"
  )

mirna <- read_tsv(p("results", "tables", "Table_S1_miRNA_differential_expression.tsv"))
circrna <- read_tsv(p("results", "tables", "Table_S2_circRNA_differential_expression.tsv"))
mrna <- read_tsv(p("results", "tables", "Table_S3_mRNA_differential_expression.tsv"))
integrated <- read_tsv(p("results", "tables", "Table_S5_integrated_ceRNA_network.tsv"))

make_edges <- function(treatment) {
  contrast <- paste0(treatment, "_vs_Vehicle")
  pairs <- miranda_pairs |> dplyr::filter(comparison == paste0(treatment, "_down"))
  mir_down <- mirna |>
    dplyr::filter(contrast == !!contrast, nominal_discovery, direction == "down") |>
    dplyr::select(miRNA = feature_id, miRNA_log2FC = log2_fold_change, miRNA_P = raw_p_value, miRNA_FDR = BH_FDR)
  mrna_up <- mrna |>
    dplyr::filter(contrast == !!contrast, nominal_discovery, direction == "up") |>
    dplyr::select(target_symbol = gene_symbol, target = gene_symbol, target_id = feature_id, target_type = assay, target_log2FC = log2_fold_change, target_P = raw_p_value)
  circ_up <- circrna |>
    dplyr::filter(contrast == !!contrast, nominal_discovery, direction == "up") |>
    dplyr::select(target_symbol = gene_symbol, target = alias, target_id = feature_id, target_type = assay, target_log2FC = log2_fold_change, target_P = raw_p_value)
  dplyr::bind_rows(
    pairs |> dplyr::inner_join(mir_down, by = "miRNA", relationship = "many-to-many") |>
      dplyr::inner_join(mrna_up, by = "target_symbol", relationship = "many-to-many"),
    pairs |> dplyr::inner_join(mir_down, by = "miRNA", relationship = "many-to-many") |>
      dplyr::inner_join(circ_up, by = "target_symbol", relationship = "many-to-many")
  ) |>
    dplyr::transmute(
      treatment = treatment,
      source = miRNA,
      target = dplyr::coalesce(target, target_symbol, target_id),
      target_symbol,
      target_type = dplyr::if_else(grepl("circRNA", target_type), "circRNA", "mRNA"),
      miRanda_score,
      hybridization_energy,
      site_count,
      miRNA_log2FC,
      miRNA_P,
      miRNA_FDR,
      target_log2FC,
      target_P
    ) |>
    dplyr::filter(!is.na(source), nzchar(source), !is.na(target), nzchar(target)) |>
    dplyr::distinct()
}

edges_il4 <- make_edges("IL4")
edges_il13 <- make_edges("IL13")
shared_keys <- dplyr::inner_join(
  edges_il4 |> dplyr::select(source, target, target_symbol, target_type),
  edges_il13 |> dplyr::select(source, target, target_symbol, target_type),
  by = c("source", "target", "target_symbol", "target_type")
) |>
  dplyr::distinct()
edges_shared <- dplyr::bind_rows(edges_il4, edges_il13) |>
  dplyr::inner_join(shared_keys, by = c("source", "target", "target_symbol", "target_type")) |>
  dplyr::distinct()

edges_focused <- integrated |>
  dplyr::filter(figure4D_focus == "yes" | (!is.na(prioritized_axis) & nzchar(prioritized_axis))) |>
  dplyr::transmute(
    treatment = "Integrated",
    source = miRNA,
    target = target,
    target_symbol = gene_symbol,
    target_type,
    miRanda_score,
    hybridization_energy = miRanda_energy,
    site_count,
    miRNA_log2FC = NA_real_,
    miRNA_P = NA_real_,
    miRNA_FDR = NA_real_,
    target_log2FC = NA_real_,
    target_P = NA_real_
  ) |>
  dplyr::filter(!is.na(source), nzchar(source), !is.na(target), nzchar(target)) |>
  dplyr::distinct()

plot_placeholder <- function(title, stem) {
  g <- ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0, y = 0, label = "No edges available", size = 5) +
    ggplot2::theme_void() +
    ggplot2::labs(title = title)
  save_figure(g, stem, width = 7, height = 5)
}

plot_network <- function(edges, title, stem, focus_nodes = character()) {
  if (!nrow(edges)) {
    plot_placeholder(title, stem)
    return(invisible(NULL))
  }
  graph_edges <- edges |>
    dplyr::transmute(from = source, to = target, target_type) |>
    dplyr::filter(!is.na(from), nzchar(from), !is.na(to), nzchar(to)) |>
    dplyr::distinct()
  if (!nrow(graph_edges)) {
    plot_placeholder(title, stem)
    return(invisible(NULL))
  }
  node_types <- dplyr::bind_rows(
    graph_edges |> dplyr::transmute(name = from, node_type = "miRNA"),
    graph_edges |> dplyr::transmute(name = to, node_type = target_type)
  ) |>
    dplyr::group_by(name) |>
    dplyr::summarise(node_type = dplyr::first(node_type), .groups = "drop")
  g <- igraph::graph_from_data_frame(graph_edges |> dplyr::select(from, to), directed = TRUE, vertices = node_types)
  colors <- c(miRNA = "#2f6db3", mRNA = "#2e8b57", circRNA = "#e69f00")
  vertex_colors <- colors[igraph::V(g)$node_type]
  vertex_colors[is.na(vertex_colors)] <- "grey75"
  focus_present <- igraph::V(g)$name %in% focus_nodes
  set.seed(14013)
  layout <- igraph::layout_with_fr(g)
  plot_one <- function() {
    plot(
      g,
      layout = layout,
      vertex.size = ifelse(focus_present, 12, 8),
      vertex.label.cex = ifelse(focus_present, 0.82, 0.62),
      vertex.label.color = "grey10",
      vertex.color = vertex_colors,
      vertex.frame.color = ifelse(focus_present, "grey10", "white"),
      edge.arrow.size = 0.28,
      edge.color = "grey60",
      main = title
    )
    legend("topleft", legend = names(colors), fill = colors, bty = "n", cex = 0.8)
  }
  grDevices::pdf(file.path(figure_dir, paste0(stem, ".pdf")), width = 9.5, height = 7.5)
  plot_one()
  grDevices::dev.off()
  grDevices::png(file.path(figure_dir, paste0(stem, ".png")), width = 2850, height = 2250, res = 300, type = "cairo-png")
  plot_one()
  grDevices::dev.off()
}

focus_nodes <- c("CEMIP", "TNC", "LIFR", "circPAPPA", "circTNC", "circLIFR", "miR-140-3p", "miR-135b-5p", "miR-625-3p")

plot_network(edges_il4, "IL-4 ceRNA network", "Figure_4A_IL4_ceRNA_network")
plot_network(edges_il13, "IL-13 ceRNA network", "Figure_4B_IL13_ceRNA_network")
plot_network(edges_shared, "Shared ceRNA network", "Figure_4C_shared_ceRNA_network")
plot_network(edges_focused, "Focused integrated ceRNA network", "Figure_4D_focused_ceRNA_network", focus_nodes = focus_nodes)

message("Figure 4A-D outputs written.")
