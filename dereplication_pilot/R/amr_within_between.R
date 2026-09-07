# R port of amr_within_between.py -- AMR gene redundancy at both levels: within-
# participant (summed) vs. whole-cohort (after the second-level/cross-participant
# clustering collapses sharing between different people too). Run from
# dereplication_pilot/.

source("R/utils.R")

if (TRUE) {  # always run main logic on source() or Rscript alike -- sys.nframe()==0 broke under nested source()
  # 1. Per participant: which first-level cluster_ids carry >=1 AMR hit, owned by whom
  amr_cluster_owner <- list()
  for (participant in all_participants()) {
    membership <- load_cluster_membership(participant, "orfs")
    groups <- unique(membership$source_group)
    calls <- purrr::map_dfr(groups, load_gene_calls_rgi)

    amr_clusters <- membership %>%
      inner_join(calls, by = c("member_id" = "gene_id")) %>%
      distinct(cluster_id) %>%
      mutate(participant = participant)
    amr_cluster_owner[[participant]] <- amr_clusters
  }
  amr_cluster_owner <- bind_rows(amr_cluster_owner)
  within_total <- nrow(amr_cluster_owner)

  # 2. Second-level clustering: first-level cluster_id (as a member) -> second-level cluster
  second_level <- read_tsv(file.path(PILOT_DIR, "cross_participant_orfs/clu_cluster.tsv"),
                            col_names = c("rep_id", "member_id"), col_types = "cc", progress = FALSE)

  resolved <- amr_cluster_owner %>%
    inner_join(second_level, by = c("cluster_id" = "member_id"))
  unresolved <- nrow(amr_cluster_owner) - nrow(resolved)

  per_whole_gene <- resolved %>%
    group_by(rep_id) %>%
    summarise(n_participants_sharing = n_distinct(participant), .groups = "drop")

  whole_cohort_distinct <- nrow(per_whole_gene)
  shared <- sum(per_whole_gene$n_participants_sharing > 1)

  cat(sprintf("within-participant distinct AMR genes: %d\n", within_total))
  cat(sprintf("whole-cohort distinct AMR genes: %d\n", whole_cohort_distinct))
  cat(sprintf("collapse from cross-participant sharing: %d (%.1f%%)\n",
              within_total - whole_cohort_distinct,
              100 * (within_total - whole_cohort_distinct) / within_total))
  cat(sprintf("shared across >=2 participants: %d (%.1f%%)\n",
              shared, 100 * shared / whole_cohort_distinct))
  if (unresolved > 0) cat(sprintf("WARNING: %d unresolved\n", unresolved))

  hist_tbl <- per_whole_gene %>% count(n_participants_sharing, name = "n_distinct_amr_genes")
  write_tsv(hist_tbl, file.path(PILOT_DIR, "R_report_amr_within_between.tsv"))
}
