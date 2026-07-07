#!/usr/bin/env Rscript

source(file.path(getwd(), "R", "final_helpers.R"))
load_pkgs(c("dplyr"))

out_dir <- ensure_dir(p("data", "processed", "final", "miranda", "inputs"))
log_dir <- ensure_dir(p("results", "final", "miranda"))

inputs <- list.files(p("data", "processed", "miranda_inputs"), full.names = TRUE)
inputs_to_copy <- inputs[!basename(inputs) %in% c("mature.fa", "mature_hsa.fa")]
for (f in inputs_to_copy) copy_if_newer(f, file.path(out_dir, basename(f)))

for (legacy_duplicate in file.path(out_dir, c("mature.fa", "mature_hsa.fa"))) {
  if (file.exists(legacy_duplicate)) unlink(legacy_duplicate)
}

fasta_files <- list.files(out_dir, pattern = "[.]fa|[.]txt$", full.names = TRUE)
checksums <- dplyr::bind_rows(lapply(fasta_files, function(f) {
  seq_lines <- readLines(f, warn = FALSE)
  seq_only <- seq_lines[!grepl("^>", seq_lines)]
  hash <- strsplit(system2("shasum", c("-a", "256", f), stdout = TRUE), "\\s+")[[1]][1]
  dplyr::tibble(
    file = sub(paste0("^", gsub("([.|()\\^{}+$*?]|\\[|\\])", "\\\\\\1", ROOT), "/?"), "", f),
    basename = basename(f),
    bytes = file.info(f)$size,
    sha256 = hash,
    sequence_line_T_count = sum(gregexpr("T", paste(seq_only, collapse = ""), fixed = TRUE)[[1]] > 0)
  )
}))
canonical_mature_hsa <- p("data", "processed", "miranda_inputs", "mature_hsa.fa")
if (!file.exists(canonical_mature_hsa)) {
  stop("Canonical human mature miRNA FASTA not found: ", canonical_mature_hsa)
}
canonical_seq_lines <- readLines(canonical_mature_hsa, warn = FALSE)
canonical_seq_only <- canonical_seq_lines[!grepl("^>", canonical_seq_lines)]
canonical_hash <- strsplit(system2("shasum", c("-a", "256", canonical_mature_hsa), stdout = TRUE), "\\s+")[[1]][1]
checksums <- dplyr::bind_rows(
  checksums,
  dplyr::tibble(
    file = "data/processed/miranda_inputs/mature_hsa.fa",
    basename = "mature_hsa.fa",
    bytes = file.info(canonical_mature_hsa)$size,
    sha256 = canonical_hash,
    sequence_line_T_count = sum(gregexpr("T", paste(canonical_seq_only, collapse = ""), fixed = TRUE)[[1]] > 0)
  )
)
mirna_query_files <- grepl("miRNA_.*_seq[.]fa$|^mature_hsa[.]fa$", checksums$basename)
if (any(checksums$sequence_line_T_count[mirna_query_files] > 0, na.rm = TRUE)) {
  bad <- checksums$basename[mirna_query_files & checksums$sequence_line_T_count > 0]
  stop("Mature miRNA query FASTA files must use RNA bases with U, not T: ", paste(bad, collapse = ", "))
}
checksums$mature_mirna_query_validation <- ifelse(
  mirna_query_files,
  "PASS: mature miRNA query sequence contains U/no T",
  "not a mature miRNA query file"
)
write_tsv(checksums, file.path(log_dir, "miranda_input_checksums.tsv"))

commands <- dplyr::tibble(
  comparison = c("IL4_down", "IL4_up", "IL13_down", "IL13_up"),
  query = c("IL4_vs_V_miRNA_down_seq.fa", "IL4_vs_V_miRNA_up_seq.fa", "IL13_vs_V_miRNA_down_seq.fa", "IL13_vs_V_miRNA_up_seq.fa"),
  target = "all_3putr_circ.fa.txt",
  command = paste("miranda", query, target, "-sc 140", sep = " "),
  score_threshold = 140,
  energy_cutoff = "miRanda default shown in historical output: Energy Threshold 1.000000 kcal/mol"
)
write_tsv(commands, file.path(log_dir, "miranda_commands.tsv"))

version_txt <- if (nzchar(Sys.which("miranda"))) {
  system2("miranda", "--version", stdout = TRUE, stderr = TRUE)
} else {
  "miranda binary not found on PATH; final workflow uses recovered miRanda v3.3a outputs."
}
writeLines(version_txt, file.path(log_dir, "miranda_version.txt"))

message("miRanda inputs and checksums written to ", out_dir)
