# R port of amr_plasmid_within_between.py -- same within-vs-between logic as
# amr_within_between.R, restricted to contig clusters that are both geNomad-plasmid-
# classified and AMR-carrying (uses R_plasmid_amr_join.tsv from plasmid_amr_join.R).
# Run from dereplication_pilot/.

source("R/utils.R")

if (TRUE) {  # always run main logic on source() or Rscript alike -- sys.nframe()==0 broke under nested source()
  plasmid_amr <- read_tsv(file.path(PILOT_DIR, "R_plasmid_amr_join.tsv"),
                            col_types = cols(.default = "c")) %>%
    mutate(is_plasmid = is_plasmid == "TRUE", has_amr = has_amr == "TRUE")

  amr_plasmid_owner <- plasmid_amr %>% filter(is_plasmid, has_amr) %>% select(participant, cluster_id)
  within_total <- nrow(amr_plasmid_owner)

  second_level <- read_tsv(file.path(PILOT_DIR, "cross_participant/clu_cluster.tsv"),
                            col_names = c("rep_id", "member_id"), col_types = "cc", progress = FALSE)

  resolved <- amr_plasmid_owner %>% inner_join(second_level, by = c("cluster_id" = "member_id"))
  unresolved <- nrow(amr_plasmid_owner) - nrow(resolved)

  per_cluster <- resolved %>%
    group_by(rep_id) %>%
    summarise(n_participants_sharing = n_distinct(participant), .groups = "drop")

  whole_cohort_distinct <- nrow(per_cluster)
  shared <- sum(per_cluster$n_participants_sharing > 1)

  cat(sprintf("within-participant distinct AMR+plasmid clusters: %d\n", within_total))
  cat(sprintf("whole-cohort distinct: %d\n", whole_cohort_distinct))
  cat(sprintf("shared across >=2 participants: %d (%.1f%%)\n",
              shared, 100 * shared / whole_cohort_distinct))
  cat(sprintf("max participants sharing one: %d\n", max(per_cluster$n_participants_sharing)))
  if (unresolved > 0) cat(sprintf("WARNING: %d unresolved\n", unresolved))

  hist_tbl <- per_cluster %>% count(n_participants_sharing, name = "n_distinct_clusters")
  write_tsv(hist_tbl, file.path(PILOT_DIR, "R_report_amr_plasmid_within_between.tsv"))
}
