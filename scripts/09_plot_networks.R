#!/usr/bin/env Rscript

# Plot simple validation layouts from documented edge tables.
# These plots are reproducible checks, not the submitted Cytoscape layout.

source(file.path(getwd(), "R", "final_helpers.R"))
load_pkgs(c("dplyr", "igraph"))

in_dir <- p("data", "processed", "final", "networks")
fig_dir <- ensure_dir(p("results", "final", "networks"))

plot_network <- function(edge_file, prefix) {
  edges <- read_tsv(file.path(in_dir, edge_file))
  if (!nrow(edges)) {
    writeLines("No edges available for this network.", file.path(fig_dir, paste0(prefix, "_EMPTY.txt")))
    return(invisible(NULL))
  }
  graph_edges <- edges |>
    dplyr::transmute(from = source, to = target, edge_type = interaction_type) |>
    dplyr::filter(!is.na(from), nzchar(from), !is.na(to), nzchar(to))
  if (!nrow(graph_edges)) {
    writeLines("No nonblank edges available for this network.", file.path(fig_dir, paste0(prefix, "_EMPTY.txt")))
    return(invisible(NULL))
  }
  g <- igraph::graph_from_data_frame(graph_edges, directed = TRUE)
  title <- gsub("_", " ", prefix)
  grDevices::pdf(file.path(fig_dir, paste0(prefix, ".pdf")), width = 10, height = 8)
  plot(g, vertex.size = 6, vertex.label.cex = 0.55, edge.arrow.size = 0.25, edge.color = "grey70", main = title)
  grDevices::dev.off()
  grDevices::png(file.path(fig_dir, paste0(prefix, ".png")), width = 2400, height = 1800, res = 250)
  plot(g, vertex.size = 6, vertex.label.cex = 0.55, edge.arrow.size = 0.25, edge.color = "grey70", main = title)
  grDevices::dev.off()
}

plot_network("restricted_historical_candidate_IL4_edges.tsv", "restricted_historical_candidate_IL4_network")
plot_network("restricted_historical_candidate_IL13_edges.tsv", "restricted_historical_candidate_IL13_network")
plot_network("restricted_historical_candidate_strict_shared_edges.tsv", "restricted_historical_candidate_strict_shared_network")
plot_network("restricted_historical_fibrosis_candidate_edges.tsv", "restricted_historical_fibrosis_candidate_network")
plot_network("biologically_prioritized_candidate_edges.tsv", "biologically_prioritized_candidate_network")
plot_network("experimentally_evaluated_candidate_edges.tsv", "experimentally_evaluated_candidate_network")
plot_network("FDR_restricted_provenance_edges.tsv", "FDR_restricted_provenance_network")

message("Network validation figures written to ", fig_dir)
