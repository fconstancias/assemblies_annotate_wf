# R port of amr_and_coassembly_reports.py:
#  1. AMR redundancy: raw RGI hits vs. distinct AMR gene clusters, per participant.
#  2. Co-assembly-unique contribution: contigs/genes in a cluster with NO single-sample
#     member at all, per participant/tier, cross-referenced against AMR hits.
# Run from dereplication_pilot/.

source("R/utils.R")

amr_redundancy_report <- function() {
  participants <- all_participants()
  per_participant <- list()

  for (participant in participants) {
    membership <- load_cluster_membership(participant, "orfs")
    groups <- unique(membership$source_group)
    calls <- purrr::map_dfr(groups, load_gene_calls_rgi)

    joined <- membership %>% inner_join(calls, by = c("member_id" = "gene_id"))
    per_participant[[participant]] <- tibble(
      participant = participant,
      raw_hits = nrow(joined),
      distinct_clusters = n_distinct(joined$cluster_id)
    )
  }

  result <- bind_rows(per_participant)
  write_tsv(result %>% select(participant, raw_hits, distinct_clusters),
            file.path(PILOT_DIR, "R_report_amr_redundancy.tsv"))

  total_raw <- sum(result$raw_hits)
  total_distinct <- sum(result$distinct_clusters)
  cat(sprintf("total raw RGI hits: %d\ntotal distinct AMR gene clusters: %d\nredundancy: %.1f%%\n",
              total_raw, total_distinct, 100 * (total_raw - total_distinct) / total_raw))
  result
}

coassembly_unique_report <- function() {
  participants <- all_participants()
  unique_rows <- list()
  amr_unique_rows <- list()

  for (participant in participants) {
    for (tier in c("contigs", "orfs")) {
      membership <- load_cluster_membership(participant, tier)
      coassembly_group <- membership$source_group[is_coassembly(membership$source_group)][1]
      if (is.na(coassembly_group)) next  # one of the 6 no-co-assembly participants

      cluster_groups <- membership %>%
        distinct(cluster_id, source_group) %>%
        group_by(cluster_id) %>%
        summarise(n_distinct_groups = n_distinct(source_group),
                  has_coassembly = any(is_coassembly(source_group)), .groups = "drop")

      ca_members <- membership %>% filter(source_group == coassembly_group)
      ca_total <- nrow(ca_members)

      unique_cluster_ids <- cluster_groups %>%
        filter(has_coassembly, n_distinct_groups == 1) %>% pull(cluster_id)
      ca_unique_members <- ca_members %>% filter(cluster_id %in% unique_cluster_ids)
      ca_unique <- nrow(ca_unique_members)

      unique_rows[[paste(participant, tier)]] <- tibble(
        participant = participant, coassembly_group = coassembly_group, tier = tier,
        unique_count = ca_unique, total_count = ca_total,
        unique_pct = round(100 * ca_unique / ca_total, 1)
      )

      if (tier == "orfs" && ca_unique > 0) {
        calls <- load_gene_calls_rgi(coassembly_group)
        hits <- ca_unique_members %>% inner_join(calls, by = c("member_id" = "gene_id"))
        if (nrow(hits) > 0) {
          amr_unique_rows[[paste(participant, tier)]] <- hits %>%
            transmute(participant, gene_id = member_id, aro_call)
        }
      }
    }
  }

  list(
    unique = bind_rows(unique_rows),
    amr_unique = bind_rows(amr_unique_rows)
  )
}

if (TRUE) {  # always run main logic on source() or Rscript alike -- sys.nframe()==0 broke under nested source()
  amr_redundancy_report()
  result <- coassembly_unique_report()
  write_tsv(result$unique, file.path(PILOT_DIR, "R_report_coassembly_unique.tsv"))
  write_tsv(result$amr_unique, file.path(PILOT_DIR, "R_report_amr_coassembly_unique.tsv"))
  cat(sprintf("co-assembly-unique: %d rows, AMR-unique hits: %d\n",
              nrow(result$unique), nrow(result$amr_unique)))
}
