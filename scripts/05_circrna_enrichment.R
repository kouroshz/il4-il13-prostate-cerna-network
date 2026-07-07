#!/usr/bin/env Rscript

source(file.path(getwd(), "R", "final_helpers.R"))
load_pkgs(c("dplyr", "ggplot2"))

out_dir <- ensure_dir(p("results", "final", "enrichment", "circrna"))
table_dir <- ensure_dir(p("data", "processed", "final", "enrichment_inputs", "circrna"))

contrasts <- c("IL4_vs_Vehicle", "IL13_vs_Vehicle")
all_summary <- list()

for (contrast in contrasts) {
  complete <- read_tsv(p("data", "processed", "final", "circrna", paste0("circRNA_", contrast, "_complete.tsv")))
  up_hosts <- complete |>
    dplyr::filter(nominal_discovery, direction == "up") |>
    dplyr::pull(gene_symbol) |>
    unique()
  background <- complete |>
    dplyr::pull(gene_symbol) |>
    unique()
  analysis_id <- paste0("circRNA_", contrast, "_up_nominal_host_genes")
  writeLines(up_hosts[!is.na(up_hosts) & nzchar(up_hosts)], file.path(table_dir, paste0(analysis_id, "_input_symbols.txt")))
  writeLines(background[!is.na(background) & nzchar(background)], file.path(table_dir, paste0(analysis_id, "_background_symbols.txt")))
  result <- run_gprofiler(up_hosts, background, analysis_id, file.path(out_dir, analysis_id), "HGNC gene symbol")
  all_summary[[analysis_id]] <- dplyr::tibble(
    assay = "circRNA microarray",
    treatment = sub("_vs_Vehicle", "", contrast),
    analysis_id = analysis_id,
    input_gene_count = length(unique(up_hosts[!is.na(up_hosts) & nzchar(up_hosts)])),
    background_gene_count = length(unique(background[!is.na(background) & nzchar(background)])),
    significant_terms = ifelse(nrow(result), sum(result$p_value < 0.05, na.rm = TRUE), 0)
  )
}

write_tsv(dplyr::bind_rows(all_summary), file.path(out_dir, "circRNA_enrichment_summary.tsv"))
message("circRNA enrichment outputs written to ", out_dir)
