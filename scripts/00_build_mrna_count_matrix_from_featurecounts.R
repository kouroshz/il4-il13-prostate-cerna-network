#!/usr/bin/env Rscript

# Optional source-tracking utility.
# Reconstructs a compact mRNA count matrix from individual featureCounts files
# when those files are available locally. The public workflow starts from the
# compact processed mRNA, miRNA, and circRNA matrices.

source(file.path(getwd(), "R", "final_helpers.R"))

metadata_file <- p("data", "metadata", "mrna_sample_metadata.tsv")
if (!file.exists(metadata_file)) {
  stop("Sample metadata not found: ", metadata_file)
}

metadata <- read_tsv(metadata_file)
required <- c("sample_id", "group", "count_file")
missing <- setdiff(required, names(metadata))
if (length(missing)) {
  stop("Sample metadata is missing required columns: ", paste(missing, collapse = ", "))
}

read_featurecounts <- function(path, sample_id) {
  if (!file.exists(path)) stop("FeatureCounts file not found: ", path)
  x <- utils::read.delim(path, comment.char = "#", check.names = FALSE, stringsAsFactors = FALSE)
  if (!"Geneid" %in% names(x)) stop("FeatureCounts file lacks Geneid column: ", path)
  count_col <- names(x)[ncol(x)]
  data.frame(Geneid = x$Geneid, sample_count = x[[count_col]], stringsAsFactors = FALSE) |>
    stats::setNames(c("Geneid", paste0("sample_", sample_id)))
}

count_tabs <- lapply(seq_len(nrow(metadata)), function(i) {
  read_featurecounts(p(metadata$count_file[i]), metadata$sample_id[i])
})
count_matrix <- Reduce(function(x, y) merge(x, y, by = "Geneid", all = TRUE), count_tabs)

out_dir <- ensure_dir(p("data", "processed", "mrna"))
out_tsv <- file.path(out_dir, "all_counts_reconstructed_from_featurecounts.tsv")
write_tsv(count_matrix, out_tsv)

checksum <- dplyr::bind_rows(lapply(metadata$count_file, function(path) {
  full_path <- p(path)
  hash <- if (file.exists(full_path)) {
    strsplit(system2("shasum", c("-a", "256", full_path), stdout = TRUE), "\\s+")[[1]][1]
  } else {
    NA_character_
  }
  dplyr::tibble(file = path, sha256 = hash)
}))
write_tsv(checksum, file.path(out_dir, "featurecounts_input_checksums.tsv"))

message("Reconstructed compact mRNA count matrix: ", out_tsv)
