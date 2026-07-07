#!/usr/bin/env Rscript

# Reconstruct documented network tables from valid inputs.
# The curated integrated network TSV is the Figure 4D/Supplementary Table X
# source. These recovered miRanda-derived R tables are retained as restricted
# candidate/provenance layers rather than complete treatment-wide networks.

source(file.path(getwd(), "R", "final_helpers.R"))
load_pkgs(c("dplyr"))

out_dir <- ensure_dir(p("data", "processed", "final", "networks"))

miranda_provenance <- paste(
  "recovered historical miRanda outputs parsed from selected fibrosis-related",
  "mRNA 3-prime UTR and circRNA target FASTA; score threshold >=140"
)
coverage_note <- paste(
  "restricted historical target coverage; not complete transcriptome-wide or",
  "complete treatment-wide discovery-list coverage"
)

read_complete <- function(assay, treatment) {
  stem <- switch(
    assay,
    mirna = paste0("miRNA_", treatment, "_vs_Vehicle_complete.tsv"),
    mrna = paste0("mRNA_", treatment, "_vs_Vehicle_complete.tsv"),
    circrna = paste0("circRNA_", treatment, "_vs_Vehicle_complete.tsv")
  )
  read_tsv(p("data", "processed", "final", assay, stem))
}

feature_status <- function(nominal, fdr) {
  dplyr::case_when(
    isTRUE(fdr) ~ "FDR_significant",
    isTRUE(nominal) ~ "nominal_discovery",
    TRUE ~ "not_threshold_defined"
  )
}

make_nodes <- function(edges) {
  if (!nrow(edges)) {
    return(data.frame(
      node_id = character(),
      node_type = character(),
      label = character(),
      treatment_membership = character(),
      discovery_status = character(),
      FDR_status = character(),
      network_scope = character()
    ))
  }
  dplyr::bind_rows(
    edges |>
      dplyr::transmute(
        node_id = source,
        node_type = "miRNA",
        label = source,
        treatment_membership = treatment,
        discovery_status = miRNA_discovery_status,
        FDR_status = ifelse(miRNA_FDR_significant, "FDR_significant", "not_FDR_significant"),
        network_scope
      ),
    edges |>
      dplyr::transmute(
        node_id = target,
        node_type = target_type,
        label = target_symbol,
        treatment_membership = treatment,
        discovery_status = target_discovery_status,
        FDR_status = ifelse(target_FDR_significant, "FDR_significant", "not_FDR_significant"),
        network_scope
      )
  ) |>
    dplyr::filter(!is.na(node_id), nzchar(node_id), !is.na(label), nzchar(label)) |>
    dplyr::distinct()
}

make_restricted_network <- function(treatment) {
  mir <- read_complete("mirna", treatment)
  mrna <- read_complete("mrna", treatment)
  circ <- read_complete("circrna", treatment)
  pairs <- read_tsv(p("data", "processed", "final", "miranda", "miranda_collapsed_mirna_target_pairs.tsv"))

  mir_down <- mir |>
    dplyr::filter(nominal_discovery, direction == "down") |>
    dplyr::transmute(
      miRNA = feature_id,
      miRNA_log2FC = log2_fold_change,
      miRNA_P = raw_p_value,
      miRNA_FDR = BH_FDR,
      miRNA_nominal_discovery = nominal_discovery,
      miRNA_FDR_significant = FDR_significant,
      miRNA_direction = direction,
      miRNA_discovery_status = dplyr::if_else(FDR_significant, "FDR_significant", "nominal_discovery")
    )

  mrna_up <- mrna |>
    dplyr::filter(nominal_discovery, direction == "up") |>
    dplyr::transmute(
      target_symbol = gene_symbol,
      target = feature_id,
      target_type = "mRNA",
      target_alias = feature_id,
      target_log2FC = log2_fold_change,
      target_P = raw_p_value,
      target_FDR = BH_FDR,
      target_nominal_discovery = nominal_discovery,
      target_FDR_significant = FDR_significant,
      target_direction = direction,
      target_discovery_status = dplyr::if_else(FDR_significant, "FDR_significant", "nominal_discovery")
    )

  circ_up <- circ |>
    dplyr::filter(nominal_discovery, direction == "up") |>
    dplyr::transmute(
      target_symbol = gene_symbol,
      target = feature_id,
      target_type = "circRNA",
      target_alias = alias,
      target_log2FC = log2_fold_change,
      target_P = raw_p_value,
      target_FDR = BH_FDR,
      target_nominal_discovery = nominal_discovery,
      target_FDR_significant = FDR_significant,
      target_direction = direction,
      target_discovery_status = dplyr::if_else(FDR_significant, "FDR_significant", "nominal_discovery")
    )

  p2 <- pairs |> dplyr::filter(comparison == paste0(treatment, "_down"))

  bind_rows(
    p2 |>
      dplyr::inner_join(mir_down, by = "miRNA", relationship = "many-to-many") |>
      dplyr::inner_join(mrna_up, by = "target_symbol", relationship = "many-to-many"),
    p2 |>
      dplyr::inner_join(mir_down, by = "miRNA", relationship = "many-to-many") |>
      dplyr::inner_join(circ_up, by = "target_symbol", relationship = "many-to-many")
  ) |>
    dplyr::transmute(
      layer = paste0("restricted_historical_candidate_", treatment),
      network_scope = "restricted_historical_fibrosis_candidate",
      treatment,
      source = miRNA,
      target,
      target_symbol,
      target_type,
      interaction_type = paste0("miRNA-", target_type),
      miRanda_score = best_score,
      hybridization_energy = best_energy,
      site_count,
      miRNA_log2FC,
      miRNA_P,
      miRNA_FDR,
      miRNA_nominal_discovery,
      miRNA_FDR_significant,
      miRNA_discovery_status,
      target_log2FC,
      target_P,
      target_FDR,
      target_nominal_discovery,
      target_FDR_significant,
      target_discovery_status,
      discovery_classification = "nominal down miRNA plus nominal up target",
      prediction_provenance = miranda_provenance,
      coverage_limitation = coverage_note,
      target_alias = dplyr::coalesce(target_alias, target_primary),
      target_raw,
      biologically_prioritized = FALSE,
      experimentally_tested = FALSE,
      experimental_result = "",
      fibrosis_rationale = "",
      notes = "Use as restricted candidate-network evidence only."
    ) |>
    dplyr::filter(!is.na(source), nzchar(source), !is.na(target), nzchar(target)) |>
    dplyr::distinct()
}

net_il4 <- make_restricted_network("IL4")
net_il13 <- make_restricted_network("IL13")
all_edges <- dplyr::bind_rows(net_il4, net_il13)

write_tsv(net_il4, file.path(out_dir, "restricted_historical_candidate_IL4_edges.tsv"))
write_tsv(make_nodes(net_il4), file.path(out_dir, "restricted_historical_candidate_IL4_nodes.tsv"))
write_tsv(net_il13, file.path(out_dir, "restricted_historical_candidate_IL13_edges.tsv"))
write_tsv(make_nodes(net_il13), file.path(out_dir, "restricted_historical_candidate_IL13_nodes.tsv"))

shared_keys <- dplyr::inner_join(
  net_il4 |> dplyr::select(source, target, target_symbol, target_type, interaction_type),
  net_il13 |> dplyr::select(source, target, target_symbol, target_type, interaction_type),
  by = c("source", "target", "target_symbol", "target_type", "interaction_type")
) |>
  dplyr::distinct()

strict_shared <- all_edges |>
  dplyr::inner_join(shared_keys, by = c("source", "target", "target_symbol", "target_type", "interaction_type")) |>
  dplyr::mutate(layer = "restricted_historical_candidate_strict_shared")

write_tsv(strict_shared, file.path(out_dir, "restricted_historical_candidate_strict_shared_edges.tsv"))
write_tsv(make_nodes(strict_shared), file.path(out_dir, "restricted_historical_candidate_strict_shared_nodes.tsv"))

fibrosis_symbols <- c("CEMIP", "KIAA1199", "TNC", "LIFR", "PAPPA")
fibrosis_edges <- all_edges |>
  dplyr::filter(target_symbol %in% fibrosis_symbols | grepl("PAPPA", target_symbol, ignore.case = TRUE)) |>
  dplyr::mutate(layer = "restricted_historical_fibrosis_candidate")

write_tsv(fibrosis_edges, file.path(out_dir, "restricted_historical_fibrosis_candidate_edges.tsv"))
write_tsv(make_nodes(fibrosis_edges), file.path(out_dir, "restricted_historical_fibrosis_candidate_nodes.tsv"))

fdr_edges <- all_edges |>
  dplyr::filter(miRNA_FDR_significant, target_FDR_significant) |>
  dplyr::mutate(layer = "FDR_restricted_provenance")

write_tsv(fdr_edges, file.path(out_dir, "FDR_restricted_provenance_edges.tsv"))
write_tsv(make_nodes(fdr_edges), file.path(out_dir, "FDR_restricted_provenance_nodes.tsv"))

lookup_pair <- function(treatment, mirna, symbol) {
  pairs <- read_tsv(p("data", "processed", "final", "miranda", "miranda_collapsed_mirna_target_pairs.tsv"))
  pairs |>
    dplyr::filter(comparison == paste0(treatment, "_down"), miRNA == mirna, target_symbol == symbol) |>
    dplyr::arrange(best_energy) |>
    dplyr::slice_head(n = 1)
}

candidate_base <- dplyr::tibble(
  source = c("hsa-miR-140-3p", "hsa-miR-140-3p", "hsa-miR-135b-5p", "hsa-miR-135b-5p", "hsa-miR-625-3p"),
  target_symbol = c("CEMIP", "CEMIP", "TNC", "TNC", "LIFR"),
  target_type = c("mRNA", "mRNA", "mRNA", "mRNA", "mRNA"),
  treatment = c("IL4", "IL13", "IL4", "IL13", "IL4"),
  fibrosis_rationale = c(
    "central CEMIP axis; selected for functional follow-up",
    "central CEMIP axis; selected for functional follow-up",
    "TNC fibrosis matrix candidate",
    "TNC fibrosis matrix candidate",
    "LIFR fibrosis signaling candidate"
  ),
  experimentally_tested = TRUE,
  experimental_result = c(
    "experimentally_supported",
    "experimentally_supported",
    "experimentally_not_supported_or_inconclusive_pending_author_confirmation",
    "experimentally_not_supported_or_inconclusive_pending_author_confirmation",
    "experimentally_not_supported_or_inconclusive_pending_author_confirmation"
  )
)

make_candidate_row <- function(row_index) {
  row <- candidate_base[row_index, ]
  mir <- read_complete("mirna", row$treatment) |>
    dplyr::filter(feature_id == row$source)
  target_row <- read_complete("mrna", row$treatment) |>
    dplyr::filter(gene_symbol == row$target_symbol) |>
    dplyr::slice_head(n = 1)
  pair <- lookup_pair(row$treatment, row$source, row$target_symbol)

  has_pair <- nrow(pair) > 0
  has_mir <- nrow(mir) > 0
  has_target <- nrow(target_row) > 0
  target_id <- if (has_target) target_row$feature_id[1] else row$target_symbol
  miRNA_status <- if (has_mir) feature_status(mir$nominal_discovery[1], mir$FDR_significant[1]) else "not_tested"
  target_status <- if (has_target) feature_status(target_row$nominal_discovery[1], target_row$FDR_significant[1]) else "not_tested"
  dplyr::tibble(
    layer = "biologically_prioritized_candidate",
    network_scope = "biologically_prioritized_candidate",
    treatment = row$treatment,
    source = row$source,
    target = target_id,
    target_symbol = row$target_symbol,
    target_type = row$target_type,
    interaction_type = paste0("miRNA-", row$target_type),
    miRanda_score = ifelse(has_pair, pair$best_score[1], NA_real_),
    hybridization_energy = ifelse(has_pair, pair$best_energy[1], NA_real_),
    site_count = ifelse(has_pair, pair$site_count[1], NA_integer_),
    miRNA_log2FC = ifelse(has_mir, mir$log2_fold_change[1], NA_real_),
    miRNA_P = ifelse(has_mir, mir$raw_p_value[1], NA_real_),
    miRNA_FDR = ifelse(has_mir, mir$BH_FDR[1], NA_real_),
    miRNA_nominal_discovery = ifelse(has_mir, mir$nominal_discovery[1], FALSE),
    miRNA_FDR_significant = ifelse(has_mir, mir$FDR_significant[1], FALSE),
    miRNA_discovery_status = miRNA_status,
    target_log2FC = ifelse(has_target, target_row$log2_fold_change[1], NA_real_),
    target_P = ifelse(has_target, target_row$raw_p_value[1], NA_real_),
    target_FDR = ifelse(has_target, target_row$BH_FDR[1], NA_real_),
    target_nominal_discovery = ifelse(has_target, target_row$nominal_discovery[1], FALSE),
    target_FDR_significant = ifelse(has_target, target_row$FDR_significant[1], FALSE),
    target_discovery_status = target_status,
    discovery_classification = paste(
      "biologically prioritized; formal miRNA/target classes:",
      miRNA_status,
      "miRNA,",
      target_status,
      "target"
    ),
    prediction_provenance = ifelse(
      has_pair,
      miranda_provenance,
      "reported historical prediction; exact miRanda output not recovered; focused rerun pending"
    ),
    coverage_limitation = ifelse(has_pair, coverage_note, "computational_provenance_pending"),
    target_alias = ifelse(has_target, target_row$gene_symbol[1], row$target_symbol),
    target_raw = ifelse(has_pair, pair$target_raw[1], ""),
    biologically_prioritized = TRUE,
    experimentally_tested = row$experimentally_tested,
    experimental_result = row$experimental_result,
    fibrosis_rationale = row$fibrosis_rationale,
    notes = ifelse(
      has_pair,
      "Candidate edge has recovered restricted historical miRanda support.",
      "Do not invent miRanda score or energy; exact computational output is pending."
    )
  )
}

candidate_edges <- dplyr::bind_rows(lapply(seq_len(nrow(candidate_base)), make_candidate_row)) |>
  dplyr::filter(!is.na(source), nzchar(source), !is.na(target), nzchar(target)) |>
  dplyr::distinct()

write_tsv(candidate_edges, file.path(out_dir, "biologically_prioritized_candidate_edges.tsv"))
write_tsv(make_nodes(candidate_edges), file.path(out_dir, "biologically_prioritized_candidate_nodes.tsv"))
write_tsv(candidate_edges |> dplyr::filter(experimentally_tested), file.path(out_dir, "experimentally_evaluated_candidate_edges.tsv"))
write_tsv(make_nodes(candidate_edges |> dplyr::filter(experimentally_tested)), file.path(out_dir, "experimentally_evaluated_candidate_nodes.tsv"))

write_tsv(all_edges, file.path(out_dir, "restricted_historical_directionally_concordant_candidate_edges.tsv"))
write_tsv(make_nodes(all_edges), file.path(out_dir, "restricted_historical_directionally_concordant_candidate_nodes.tsv"))

crosswalk <- dplyr::tibble(
  superseded_or_legacy_name = c(
    "IL4_exploratory_*",
    "IL13_exploratory_*",
    "strict_shared_*",
    "fibrosis_focused_*",
    "FDR_provenance_*"
  ),
  current_name = c(
    "restricted_historical_candidate_IL4_*",
    "restricted_historical_candidate_IL13_*",
    "restricted_historical_candidate_strict_shared_*",
    "restricted_historical_fibrosis_candidate_*",
    "FDR_restricted_provenance_*"
  ),
  reason = "Current names explicitly state restricted historical miRanda coverage."
)
write_tsv(crosswalk, file.path(out_dir, "network_table_name_crosswalk.tsv"))

message("Restricted network tables written to ", out_dir)
