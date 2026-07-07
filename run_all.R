#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
root <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[1]), mustWork = TRUE)) else getwd()
setwd(root)
Sys.setenv(IL4_IL13_PROJECT_ROOT = root)

scripts <- c(
  "scripts/01_mrna_differential_expression.R",
  "scripts/02_mirna_differential_expression.R",
  "scripts/03_circrna_differential_expression.R",
  "scripts/04_mrna_enrichment.R",
  "scripts/05_circrna_enrichment.R",
  "scripts/06_prepare_miranda_inputs.R",
  "scripts/07_parse_miranda_results.R",
  "scripts/08_construct_networks.R",
  "scripts/09_plot_networks.R",
  "scripts/10_sync_integrated_network_source.R"
)

for (script in scripts) {
  message("\n== Running ", script, " ==")
  status <- system2(file.path(R.home("bin"), "Rscript"), script)
  if (!identical(status, 0L)) stop("Script failed: ", script, call. = FALSE)
}

source(file.path(root, "R", "final_helpers.R"))
write_session_info(file.path(root, "sessionInfo.txt"))
message("\nAll final computational workflow steps completed.")
