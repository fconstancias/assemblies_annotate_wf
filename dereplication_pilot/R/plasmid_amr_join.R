# R port of plasmid_amr_join.py -- for each contig CLUSTER (not raw contig), was any
# member ever classified plasmid/virus by geNomad, and does any member carry a real AMR
# gene hit? Run from dereplication_pilot/.

source("R/utils.R")

load_contig_has_amr <- function(group) {
  calls <- load_gene_calls_rgi(group)
  if (nrow(calls) == 0) return(character())
  gff <- read_gff3(group)
  gff %>% semi_join(calls, by = "gene_id") %>% pull(contig) %>% unique()
}

if (TRUE) {  # always run main logic on source() or Rscript alike -- sys.nframe()==0 broke under nested source()
  out_rows <- list()

  for (participant in all_participants()) {
    membership <- load_cluster_membership(participant, "contigs")
    groups <- unique(membership$source_group)

    genomad_calls <- purrr::map_dfr(groups, load_contig_calls_genomad)
    amr_contigs <- unique(unlist(purrr::map(groups, load_contig_has_amr)))

    cluster_status <- membership %>%
      left_join(genomad_calls, by = c("member_id" = "contig_id")) %>%
      mutate(is_plasmid = str_detect(coalesce(call, ""), "plasmid"),
             is_virus = str_detect(coalesce(call, ""), "virus"),
             has_amr = member_id %in% amr_contigs) %>%
      group_by(cluster_id) %>%
      summarise(is_plasmid = any(is_plasmid), is_virus = any(is_virus),
                has_amr = any(has_amr), n_members = n(), .groups = "drop") %>%
      mutate(participant = participant)

    out_rows[[participant]] <- cluster_status
  }

  result <- bind_rows(out_rows) %>% select(participant, cluster_id, is_plasmid, is_virus, has_amr, n_members)
  write_tsv(result, file.path(PILOT_DIR, "R_plasmid_amr_join.tsv"))

  summary <- result %>%
    mutate(category = case_when(is_plasmid ~ "plasmid", is_virus ~ "virus", TRUE ~ "other")) %>%
    group_by(category) %>%
    summarise(n_clusters = n(), n_amr_carrying = sum(has_amr), .groups = "drop")
  write_tsv(summary, file.path(PILOT_DIR, "R_report_plasmid_amr.tsv"))
  print(summary)
}
