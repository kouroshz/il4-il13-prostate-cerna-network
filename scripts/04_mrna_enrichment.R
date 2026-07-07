#!/usr/bin/env Rscript

source(file.path(getwd(), "R", "final_helpers.R"))
load_pkgs(c("dplyr", "ggplot2"))

out_dir <- ensure_dir(p("results", "final", "enrichment", "mrna"))
table_dir <- ensure_dir(p("data", "processed", "final", "enrichment_inputs", "mrna"))

contrasts <- c("IL4_vs_Vehicle", "IL13_vs_Vehicle")
all_summary <- list()

for (contrast in contrasts) {
  complete <- read_tsv(p("data", "processed", "final", "mrna", paste0("mRNA_", contrast, "_complete.tsv")))
  up <- complete |>
    dplyr::filter(nominal_discovery, direction == "up") |>
    dplyr::pull(feature_id) |>
    unique()
  background <- complete$feature_id |> unique()
  analysis_id <- paste0("mRNA_", contrast, "_up_nominal")
  writeLines(up, file.path(table_dir, paste0(analysis_id, "_input_ensembl.txt")))
  writeLines(background, file.path(table_dir, paste0(analysis_id, "_background_ensembl.txt")))
  result <- run_gprofiler(up, background, analysis_id, file.path(out_dir, analysis_id), "Ensembl gene ID")
  all_summary[[analysis_id]] <- dplyr::tibble(
    assay = "mRNA-seq",
    treatment = sub("_vs_Vehicle", "", contrast),
    analysis_id = analysis_id,
    input_gene_count = length(up),
    background_gene_count = length(background),
    significant_terms = ifelse(nrow(result), sum(result$p_value < 0.05, na.rm = TRUE), 0)
  )
}

write_tsv(dplyr::bind_rows(all_summary), file.path(out_dir, "mRNA_enrichment_summary.tsv"))
message("mRNA enrichment outputs written to ", out_dir)
