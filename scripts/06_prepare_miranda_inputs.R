#!/usr/bin/env Rscript

# Validate miRanda input FASTA files and record run metadata.

source(file.path(getwd(), "R", "final_helpers.R"))
load_pkgs(c("dplyr"))

log_dir <- ensure_dir(p("results", "logs"))
input_dir <- p("data", "processed", "miranda_inputs")

fasta_files <- list.files(input_dir, pattern = "[.]fa$|[.]txt$", full.names = TRUE)
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
write_tsv(checksums, file.path(log_dir, "miRanda_input_checksums.tsv"))

commands <- dplyr::tibble(
  comparison = c("IL4_down", "IL4_up", "IL13_down", "IL13_up"),
  query = c("IL4_vs_V_miRNA_down_seq.fa", "IL4_vs_V_miRNA_up_seq.fa", "IL13_vs_V_miRNA_down_seq.fa", "IL13_vs_V_miRNA_up_seq.fa"),
  target = "all_3putr_circ.fa.txt",
  recovered_output_source = c(
    "IL4_vs_V_miRNA_down_against_3utr_circ.txt",
    "IL4_vs_V_miRNA_up_against_3utr_circ.txt",
    "IL13_vs_V_miRNA_down_against_3utr_circ.txt",
    "IL13_vs_V_miRNA_up_against_3utr_circ.txt"
  ),
  miranda_version = "v3.3a",
  command = paste("miranda", query, target, "-sc 140 -en 1.0 -scale 4 -go -9 -ge -4", sep = " "),
  score_rule = "minimum miRanda alignment score \u2265 140",
  score_threshold = 140,
  energy_threshold_kcal_per_mol = 1.0,
  gap_open_penalty = -9,
  gap_extend_penalty = -4,
  scaling_parameter = 4,
  mature_miRNA_sequence_rule = "RNA sequence with U bases; T bases are not allowed in mature miRNA query FASTA"
)
write_tsv(commands, file.path(log_dir, "miRanda_commands.tsv"))

version_txt <- if (nzchar(Sys.which("miranda"))) {
  system2("miranda", "--version", stdout = TRUE, stderr = TRUE)
} else {
  "miRanda binary not found on PATH; workflow uses recovered miRanda v3.3a outputs."
}
writeLines(version_txt, file.path(log_dir, "miRanda_version.txt"))

message("miRanda input checks and run metadata written to ", log_dir)
