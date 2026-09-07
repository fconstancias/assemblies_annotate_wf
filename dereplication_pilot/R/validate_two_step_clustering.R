# Direct empirical test of the two-step (within-participant, then cross-participant)
# clustering's known theoretical limitation: at the second level, only one representative
# per within-participant cluster gets compared, so a real cross-participant match can be
# missed if that specific representative wasn't quite close enough (even though another
# member of the same cluster would have matched fine). This bias is directionally
# understood (conservative -- can only cause undercounting, never inflate sharing), but
# not previously measured on real data.
#
# Method: find AMR calls (aro_call) where RGI agrees it's the same named gene/allele but
# the whole-cohort clustering still produced >1 distinct cluster -- candidates for a
# missed merge. Pull one representative protein sequence per candidate cluster, pool them
# ALL (across every aro_call) into one small fasta, and re-cluster that pool directly
# (single-step, no hierarchy) with the same settings as the original ORF-tier pass. If
# candidates that were split in the original two-step approach merge back together here,
# that's real, direct evidence of a missed merge -- not just a theoretical possibility.
# Run from dereplication_pilot/.

source("R/utils.R")

master <- read_tsv(file.path(PILOT_DIR, "master_orf_table_amr.tsv.gz"), col_types = cols(.default = "c"))

multi_cluster_calls <- master %>%
  distinct(aro_call, orf_whole_cohort_id) %>%
  count(aro_call, name = "n_original_clusters") %>%
  filter(n_original_clusters > 1)

cat(sprintf("%d aro_call values have >1 distinct whole-cohort cluster (candidates)\n", nrow(multi_cluster_calls)))

# aro_call <-> orf_whole_cohort_id is a real many-to-many mapping (the same whole-cohort
# ORF cluster can carry >1 distinct aro_call if RGI named different member instances
# slightly differently -- that's the ambiguity this test is partly measuring, not a bug).
# Keep that full mapping for the analysis join, but slice the representative SEQUENCE
# once per orf_whole_cohort_id only -- slicing per (aro_call, orf_whole_cohort_id) instead
# duplicated ~230 clusters into the FASTA under an identical header (same
# orf_whole_cohort_id, written once per aro_call it happened to carry), corrupting the
# re-clustering step with literal duplicate input records under one header.
call_membership <- master %>%
  semi_join(multi_cluster_calls, by = "aro_call") %>%
  distinct(aro_call, orf_whole_cohort_id)
write_tsv(call_membership, file.path(PILOT_DIR, "validate_two_step_call_membership.tsv"))

candidates <- master %>%
  semi_join(multi_cluster_calls, by = "aro_call") %>%
  distinct(orf_whole_cohort_id, orf_id, sample) %>%
  group_by(orf_whole_cohort_id) %>%
  slice(1) %>%  # one representative instance per candidate cluster
  ungroup()

cat(sprintf("%d distinct candidate representative sequences to test\n", nrow(candidates)))

# Extract each candidate's real protein sequence from its group_export/*.faa (grouped
# by source group so each .faa is only opened once, not once per candidate).
extract_sequences <- function(candidates) {
  groups <- unique(candidates$sample)
  seqs <- list()
  for (g in groups) {
    ids_needed <- candidates$orf_id[candidates$sample == g]
    faa_path <- if (is_coassembly(g)) file.path(COASSEMBLY_GENE_EXPORT, paste0(g, ".faa"))
                else file.path(SINGLE_GENE_EXPORT, paste0(g, ".faa"))
    if (!file.exists(faa_path)) next
    lines <- read_lines(faa_path, progress = FALSE)
    header_idx <- which(str_starts(lines, ">"))
    ids <- str_remove(lines[header_idx], "^>")
    want <- which(ids %in% ids_needed)
    for (w in want) {
      start <- header_idx[w] + 1
      end <- if (w < length(header_idx)) header_idx[w + 1] - 1 else length(lines)
      seqs[[ids[w]]] <- paste(lines[start:end], collapse = "")
    }
  }
  seqs
}

cat("extracting real protein sequences...\n")
seqs <- extract_sequences(candidates)
cat(sprintf("extracted %d of %d requested sequences\n", length(seqs), nrow(candidates)))

candidates <- candidates %>% mutate(seq = unlist(seqs[orf_id]))
candidates <- candidates %>% filter(!is.na(seq))

fasta_path <- file.path(PILOT_DIR, "validate_two_step_candidates.faa")
fasta_lines <- character(2 * nrow(candidates))
fasta_lines[seq(1, length(fasta_lines), 2)] <- paste0(">", candidates$orf_whole_cohort_id)
fasta_lines[seq(2, length(fasta_lines), 2)] <- candidates$seq
writeLines(fasta_lines, fasta_path)
cat(sprintf("wrote %d candidate sequences to %s\n", nrow(candidates), fasta_path))
write_tsv(candidates %>% select(-seq), file.path(PILOT_DIR, "validate_two_step_candidates.tsv"))
