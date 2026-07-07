#!/usr/bin/env Rscript

source(file.path(getwd(), "R", "final_helpers.R"))
load_pkgs(c("dplyr"))

out_dir <- ensure_dir(p("data", "processed", "final", "miranda"))
src_dir <- p("data", "processed", "miranda_outputs")

parse_file <- function(path, comparison) {
  x <- utils::read.table(path, sep = "\t", header = FALSE, fill = TRUE, quote = "", comment.char = "")
  x <- x[grepl("^>", x$V1) & !grepl("^>>", x$V1), , drop = FALSE]
  if (!nrow(x)) return(data.frame())
  x$V1 <- sub("^>", "", x$V1)
  mir_pos <- strsplit(as.character(x$V5), " ")
  target_pos <- strsplit(as.character(x$V6), " ")
  dplyr::tibble(
    comparison = comparison,
    miRNA = as.character(x$V1),
    target_raw = as.character(x$V2),
    target_primary = sub("[|:].*", "", target_raw),
    target_symbol = ifelse(grepl("[|:]", target_raw), sub(".*[|:]", "", target_raw), target_raw),
    score = as.numeric(x$V3),
    energy = as.numeric(x$V4),
    mirna_start = as.numeric(vapply(mir_pos, `[`, character(1), 1)),
    mirna_end = as.numeric(vapply(mir_pos, `[`, character(1), 2)),
    target_start = as.numeric(vapply(target_pos, `[`, character(1), 1)),
    target_end = as.numeric(vapply(target_pos, `[`, character(1), 2)),
    source_file = basename(path)
  )
}

files <- c(
  IL4_down = "IL4_vs_V_miRNA_down_against_3utr_circ_filtered_targets.txt",
  IL4_up = "IL4_vs_V_miRNA_up_against_3utr_circ_filtered_targets.txt",
  IL13_down = "IL13_vs_V_miRNA_down_against_3utr_circ_filtered_targets.txt",
  IL13_up = "IL13_vs_V_miRNA_up_against_3utr_circ_filtered_targets.txt"
)

site_level <- dplyr::bind_rows(lapply(names(files), function(nm) parse_file(file.path(src_dir, files[[nm]]), nm)))
write_tsv(site_level, file.path(out_dir, "miranda_site_level_predictions.tsv"))

collapsed <- site_level |>
  dplyr::group_by(comparison, miRNA, target_raw, target_primary, target_symbol) |>
  dplyr::arrange(energy, dplyr::desc(score), .by_group = TRUE) |>
  dplyr::summarise(
    best_score = max(score, na.rm = TRUE),
    best_energy = min(energy, na.rm = TRUE),
    site_count = dplyr::n(),
    best_target_start = target_start[which.min(energy)][1],
    best_target_end = target_end[which.min(energy)][1],
    source_files = paste(unique(source_file), collapse = ";"),
    .groups = "drop"
  )
write_tsv(collapsed, file.path(out_dir, "miranda_collapsed_mirna_target_pairs.tsv"))

message("Parsed miRanda outputs to ", out_dir)
