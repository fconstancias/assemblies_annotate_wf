# R port of validate_propagation.py -- for clusters with >=2 members that independently
# received a real hit from RGI (genes) or geNomad (contigs), do they agree? Direct test
# of whether propagating one cluster member's annotation to the rest would be correct.
#
# Run from dereplication_pilot/ (source("R/utils.R") assumes that working directory).

source("R/utils.R")
library(purrr)

# Note: member_id (gene_id/contig_id) is already globally unique across every group in
# this project (group-prefixed headers, e.g. mh_p110___0/spaS30___0) -- no need to carry
# source_group through the calls tables just to join; it already lives on `membership`.
validate_tier <- function(tier, load_calls_fn) {
  participants <- all_participants()
  results <- list()
  disagreements <- list()

  for (participant in participants) {
    membership <- load_cluster_membership(participant, tier)
    groups <- unique(membership$source_group)

    calls_long <- map_dfr(groups, load_calls_fn)  # columns: (gene_id|contig_id), (aro_call|call)
    id_field <- names(calls_long)[1]

    joined <- membership %>%
      inner_join(calls_long, by = setNames(id_field, "member_id")) %>%
      rename(call = !!names(calls_long)[2])

    cluster_calls <- joined %>%
      group_by(cluster_id) %>%
      summarise(n_hits = n(), n_distinct_calls = n_distinct(call), .groups = "drop") %>%
      filter(n_hits >= 2)

    testable <- nrow(cluster_calls)
    consistent <- sum(cluster_calls$n_distinct_calls == 1)
    inconsistent <- testable - consistent

    results[[participant]] <- tibble(participant = participant, testable = testable,
                                      consistent = consistent, inconsistent = inconsistent)

    if (inconsistent > 0) {
      bad_clusters <- cluster_calls %>% filter(n_distinct_calls > 1) %>% pull(cluster_id)
      disagreements[[participant]] <- joined %>%
        filter(cluster_id %in% bad_clusters) %>%
        select(cluster_id, source_group, member_id, call) %>%
        mutate(participant = participant)
    }
  }

  list(summary = bind_rows(results), disagreements = bind_rows(disagreements))
}

if (TRUE) {  # always run main logic on source() or Rscript alike -- sys.nframe()==0 broke under nested source()
  orfs_result <- validate_tier("orfs", load_gene_calls_rgi)
  contigs_result <- validate_tier("contigs", load_contig_calls_genomad)

  bind_rows(
    orfs_result$summary %>% mutate(tool = "RGI"),
    contigs_result$summary %>% mutate(tool = "geNomad")
  ) %>%
    select(tool, participant, testable, consistent, inconsistent) %>%
    write_tsv(file.path(PILOT_DIR, "R_report_validation.tsv"))

  write_tsv(orfs_result$disagreements, file.path(PILOT_DIR, "R_validate_orfs_disagreements.tsv"))
  write_tsv(contigs_result$disagreements, file.path(PILOT_DIR, "R_validate_contigs_disagreements.tsv"))

  cat("=== RGI (orfs) ===\n")
  print(orfs_result$summary %>% summarise(testable = sum(testable), consistent = sum(consistent), inconsistent = sum(inconsistent)))
  cat("=== geNomad (contigs) ===\n")
  print(contigs_result$summary %>% summarise(testable = sum(testable), consistent = sum(consistent), inconsistent = sum(inconsistent)))
}
