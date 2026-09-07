#!/usr/bin/env Rscript
# R port of analyze_clusters.py -- analyzes an mmseqs easy-cluster _cluster.tsv output
# for cross-assembly redundancy. Parses cluster membership, attributes each member back
# to its source assembly group via its header prefix (headers in this project are
# already globally unique and group-prefixed), and reports how much of a participant's
# own co-assembly + single-sample contigs/genes are redundant with each other.
#
# Usage: Rscript R/analyze_clusters.R --cluster-tsv X --fasta Y --groups g1 g2 ...
#          --tier contigs|orfs --participant pNNN --out cluster_membership.tsv --summary summary.txt

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
})

# Minimal CLI parser (optparse isn't installed in this env): single-value flags take
# exactly the next token; --groups is the one multi-value flag (like argparse's
# nargs="+") and consumes every token up to the next --flag, matching how
# run_pilot.sh/run_sweep.sh invoke this with an unquoted variable expansion.
parse_args <- function(args) {
  single_flags <- c("cluster-tsv", "fasta", "tier", "participant", "out", "summary")
  out <- list()
  i <- 1
  while (i <= length(args)) {
    flag <- sub("^--", "", args[i])
    if (flag == "groups") {
      j <- i + 1
      vals <- character()
      while (j <= length(args) && !str_starts(args[j], "--")) {
        vals <- c(vals, args[j]); j <- j + 1
      }
      out$groups <- vals
      i <- j
    } else if (flag %in% single_flags) {
      out[[flag]] <- args[i + 1]
      i <- i + 2
    } else {
      i <- i + 1
    }
  }
  out
}

read_fasta_lengths <- function(path) {
  lines <- read_lines(path, progress = FALSE)
  header_idx <- which(str_starts(lines, ">"))
  names <- str_extract(lines[header_idx], "^>(\\S+)", group = 1)
  ends <- c(header_idx[-1] - 1, length(lines))
  seq_lines_per_header <- ends - header_idx  # number of sequence lines following each header
  # sum of nchar for the sequence lines belonging to each header, vectorized via a repeated group index
  group_id <- rep(seq_along(header_idx), times = seq_lines_per_header)
  seq_line_idx <- setdiff(seq_len(length(lines)), header_idx)
  seq_lens <- nchar(lines[seq_line_idx])
  lengths <- as.integer(tapply(seq_lens, group_id, sum))
  # tapply drops groups with zero sequence lines (empty record) -- restore as 0
  full <- integer(length(header_idx))
  full[as.integer(names(table(group_id)))] <- lengths
  setNames(full, names)
}

# Every real group name in this project matches ^(mh_p\d+|spaS\d+)$ -- a single
# vectorized regex extracts the owning group directly (same fix already applied in
# cross_participant_analysis.R, which found the naive loop-over-groups approach far too
# slow at millions of rows).
group_of <- function(seq_ids, groups_sorted_desc) {
  str_match(seq_ids, "^(mh_p\\d+|spaS\\d+)")[, 2]
}

opt <- parse_args(commandArgs(trailingOnly = TRUE))
opt$`cluster-tsv` -> cluster_tsv_path  # keep the hyphenated flag name readable below
groups <- unique(opt$groups)

lengths <- read_fasta_lengths(opt$fasta)

clusters <- read_tsv(cluster_tsv_path, col_names = c("rep_id", "member_id"), col_types = "cc", progress = FALSE) %>%
  mutate(source_group = group_of(member_id, groups),
         length = unname(lengths[member_id]))

unresolved <- sum(is.na(clusters$source_group))
clusters$source_group[is.na(clusters$source_group)] <- "UNKNOWN"

cluster_summary <- clusters %>%
  group_by(rep_id) %>%
  summarise(n_members = n(), n_distinct_groups = n_distinct(source_group), .groups = "drop") %>%
  mutate(is_cross_group = n_distinct_groups > 1)

clusters <- clusters %>% left_join(cluster_summary %>% select(rep_id, is_cross_group), by = "rep_id")

write_tsv(clusters %>% transmute(cluster_id = rep_id, member_id, source_group, length,
                                  cross_group_cluster = as.integer(is_cross_group)),
          opt$out)

total_seqs <- nrow(clusters)
total_bp <- sum(clusters$length)
n_singleton <- sum(cluster_summary$n_members == 1)
n_multi <- sum(cluster_summary$n_members > 1)
n_cross_group <- sum(cluster_summary$is_cross_group)
redundant <- clusters %>% filter(is_cross_group)
redundant_seqs <- nrow(redundant)
redundant_bp <- sum(redundant$length)

per_group <- clusters %>%
  group_by(source_group) %>%
  summarise(total_count = n(), total_bp = sum(length),
            redundant_count = sum(is_cross_group), redundant_bp = sum(length[is_cross_group]),
            .groups = "drop") %>%
  arrange(source_group)

lines <- c(
  sprintf("=== %s / %s ===", opt$participant, opt$tier),
  sprintf("total sequences pooled: %d  (%d bp/residues)", total_seqs, total_bp),
  sprintf("total clusters: %d  (singleton: %d, multi-member: %d)", nrow(cluster_summary), n_singleton, n_multi),
  sprintf("cross-group clusters (redundant across >=2 assemblies): %d", n_cross_group),
  if (total_seqs > 0) sprintf("redundant sequences: %d/%d (%.1f%% by count, %.1f%% by bp)",
                               redundant_seqs, total_seqs, 100 * redundant_seqs / total_seqs, 100 * redundant_bp / total_bp) else NULL,
  "per-group breakdown (redundant/total, count and bp):",
  sprintf("  %s: %d/%d (%.1f%%) by count, %d/%d (%.1f%%) by bp",
          per_group$source_group, per_group$redundant_count, per_group$total_count,
          100 * per_group$redundant_count / per_group$total_count,
          per_group$redundant_bp, per_group$total_bp,
          100 * per_group$redundant_bp / per_group$total_bp)
)
if (unresolved > 0) {
  lines <- c(lines, sprintf("WARNING: %d member IDs did not match any known group prefix", unresolved))
}

summary_text <- paste(lines, collapse = "\n")
writeLines(summary_text, opt$summary)
cat(summary_text, "\n")

if (unresolved > 0) quit(status = 1)
