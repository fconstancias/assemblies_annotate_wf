# Analysis half of the two-step clustering validation test (extraction done by
# validate_two_step_clustering.R). Reads the direct single-step re-clustering of the
# candidate pool (validate_two_step_clu_cluster.tsv, one representative sequence per
# original whole-cohort ORF cluster, fasta-headered by that original orf_whole_cohort_id)
# and checks: for each aro_call that the two-step pipeline split into >1 distinct
# orf_whole_cohort_id, do any of those original clusters now fall under the SAME
# single-step cluster? That's direct, real evidence of a merge the two-step hierarchy
# missed (only one representative per within-participant cluster is compared at the
# cross-participant level, so a real match can be missed if that specific representative
# wasn't quite close enough -- even though this whole exercise is symmetric proof, not
# just theory). Written up as a report table + a Rmd-ready summary object.
# Run from dereplication_pilot/.

source("R/utils.R")

# call_membership.tsv carries the real many-to-many aro_call <-> orf_whole_cohort_id
# mapping; candidates.tsv itself is deduplicated to one row per orf_whole_cohort_id
# (the sequence actually tested), so the aro_call join has to come from the former.
candidates <- read_tsv(file.path(PILOT_DIR, "validate_two_step_candidates.tsv"),
                         col_types = cols(.default = "c"))
call_membership <- read_tsv(file.path(PILOT_DIR, "validate_two_step_call_membership.tsv"),
                              col_types = cols(.default = "c"))

reclu <- read_tsv(file.path(PILOT_DIR, "validate_two_step_clu_cluster.tsv"),
                    col_names = c("new_cluster_rep", "orf_whole_cohort_id"),
                    col_types = "cc", progress = FALSE)

cat(sprintf("candidates tested: %d representative sequences, %d distinct aro_call values\n",
            nrow(candidates), n_distinct(call_membership$aro_call)))
cat(sprintf("single-step re-clustering produced %d new clusters over these sequences\n",
            n_distinct(reclu$new_cluster_rep)))

# Per aro_call: did the direct single-step pass merge >=2 of its original (two-step)
# whole-cohort clusters back into one? If every original cluster keeps its own distinct
# new_cluster_rep, the two-step split is confirmed -- these ORFs really were distinct
# enough even under direct comparison, not an artifact of the hierarchy.
per_call <- call_membership %>%
  left_join(reclu, by = "orf_whole_cohort_id") %>%
  group_by(aro_call) %>%
  summarise(
    n_original_clusters = n_distinct(orf_whole_cohort_id),
    n_new_clusters = n_distinct(new_cluster_rep),
    merged = n_new_clusters < n_original_clusters,
    .groups = "drop"
  )

n_confirmed_split <- sum(!per_call$merged)
n_missed_merge <- sum(per_call$merged)

cat(sprintf("\n%d of %d aro_call values: two-step split CONFIRMED (every original cluster stayed distinct under direct re-clustering)\n",
            n_confirmed_split, nrow(per_call)))
cat(sprintf("%d of %d aro_call values: two-step split shows >=1 MISSED MERGE (direct re-clustering pulled originally-distinct clusters together)\n",
            n_missed_merge, nrow(per_call)))

# Pair-level detail for the missed-merge cases, for the report table: which specific
# original clusters merged, and how many total ORF instances (across all samples/
# participants, from the full master table -- not just the one representative tested)
# does each side represent, so the reader can judge real-world weight, not just count.
master <- read_tsv(file.path(PILOT_DIR, "master_orf_table_amr.tsv.gz"), col_types = cols(.default = "c"))
cluster_sizes <- master %>% count(orf_whole_cohort_id, name = "n_orf_instances")

missed_detail <- call_membership %>%
  left_join(reclu, by = "orf_whole_cohort_id") %>%
  semi_join(per_call %>% filter(merged), by = "aro_call") %>%
  left_join(cluster_sizes, by = "orf_whole_cohort_id") %>%
  arrange(aro_call, new_cluster_rep)

write_tsv(per_call, file.path(PILOT_DIR, "R_two_step_validation_summary.tsv"))
write_tsv(missed_detail, file.path(PILOT_DIR, "R_two_step_validation_missed_merges.tsv"))

cat(sprintf("\nwrote R_two_step_validation_summary.tsv (%d aro_call rows)\n", nrow(per_call)))
cat(sprintf("wrote R_two_step_validation_missed_merges.tsv (%d rows, %d aro_call values affected)\n",
            nrow(missed_detail), n_distinct(missed_detail$aro_call)))
