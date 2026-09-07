# R port of cross_participant_analysis.py -- within- vs. between-participant sharing at
# the second clustering level (cross_participant/ or cross_participant_orfs/, which
# cluster each participant's own already-deduplicated representatives against everyone
# else's). Run from dereplication_pilot/.

source("R/utils.R")

load_group_to_participant <- function() {
  purrr::map_dfr(all_participants(), function(p) {
    load_cluster_membership(p, "contigs") %>% distinct(source_group) %>% mutate(participant = p)
  }) %>% distinct(source_group, .keep_all = TRUE)
}

# member_id (a first-level representative ID) is group-prefixed, e.g.
# mh_p110_000000000011 or mh_p110___0 -- every real group name in this project matches
# ^(mh_p\d+|spaS\d+)$, so a single vectorized regex extracts the owning group directly
# without looping over ~300 known groups x millions of rows (that loop-based approach,
# tried first, would mean up to ~1.5B string comparisons -- switched to this instead).
group_of_member <- function(member_ids) {
  str_match(member_ids, "^(mh_p\\d+|spaS\\d+)")[, 2]
}

analyze_cross_tier <- function(tier, cluster_tsv_path, g2p) {
  second_level <- read_tsv(cluster_tsv_path, col_names = c("rep_id", "member_id"),
                            col_types = "cc", progress = FALSE) %>%
    mutate(owning_group = group_of_member(member_id)) %>%
    left_join(g2p, by = c("owning_group" = "source_group"))

  per_cluster <- second_level %>%
    group_by(rep_id) %>%
    summarise(n_participants_in_cluster = n_distinct(participant), .groups = "drop")

  hist_tbl <- per_cluster %>% count(n_participants_in_cluster, name = "n_clusters") %>% mutate(tier = tier)
  list(hist = hist_tbl, second_level = second_level, per_cluster = per_cluster)
}

if (TRUE) {  # always run main logic on source() or Rscript alike -- sys.nframe()==0 broke under nested source()
  g2p <- load_group_to_participant()

  contigs_result <- analyze_cross_tier("contigs", file.path(PILOT_DIR, "cross_participant/clu_cluster.tsv"), g2p)
  orfs_result <- analyze_cross_tier("orfs", file.path(PILOT_DIR, "cross_participant_orfs/clu_cluster.tsv"), g2p)

  bind_rows(contigs_result$hist, orfs_result$hist) %>%
    write_tsv(file.path(PILOT_DIR, "R_report_cross_participant.tsv"))

  for (r in list(contigs = contigs_result, orfs = orfs_result)) {
    total <- sum(r$hist$n_clusters)
    within <- sum(r$hist$n_clusters[r$hist$n_participants_in_cluster == 1])
    cat(sprintf("total: %d, within-only: %d (%.1f%%), between: %d (%.1f%%)\n",
                total, within, 100 * within / total, total - within, 100 * (total - within) / total))
  }

  saveRDS(list(contigs = contigs_result, orfs = orfs_result), file.path(PILOT_DIR, "R_cross_participant_intermediate.rds"))
}
