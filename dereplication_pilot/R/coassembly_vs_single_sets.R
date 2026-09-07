# Set breakdown (unique-to-single, unique-to-coassembly, shared) per participant/tier --
# the data behind the UpSet-style comparison requested for report section 2/3. Extends
# report_coassembly_unique.tsv (which only reported the co-assembly side) with the
# single-sample side and the explicit 3-way breakdown, for both tiers.
# Run from dereplication_pilot/.

source("R/utils.R")

out_rows <- list()

for (participant in all_participants()) {
  for (tier in c("contigs", "orfs")) {
    membership <- load_cluster_membership(participant, tier)
    coassembly_group <- membership$source_group[is_coassembly(membership$source_group)][1]
    if (is.na(coassembly_group)) next  # one of the 6 no-co-assembly participants

    cluster_kind <- membership %>%
      distinct(cluster_id, source_group) %>%
      group_by(cluster_id) %>%
      summarise(has_coassembly = any(is_coassembly(source_group)),
                has_single = any(!is_coassembly(source_group)), .groups = "drop") %>%
      mutate(set = case_when(
        has_coassembly & has_single ~ "shared",
        has_coassembly ~ "coassembly_only",
        TRUE ~ "single_only"
      ))

    counts <- cluster_kind %>% count(set, name = "n_clusters")
    out_rows[[paste(participant, tier)]] <- counts %>%
      mutate(participant = participant, tier = tier) %>%
      select(participant, tier, set, n_clusters)
  }
}

result <- bind_rows(out_rows)
write_tsv(result, file.path(PILOT_DIR, "report_coassembly_vs_single_sets.tsv"))

cat("=== aggregated across all 19 co-assembly participants ===\n")
result %>% group_by(tier, set) %>% summarise(n = sum(n_clusters), .groups = "drop") %>% print(n = 20)
