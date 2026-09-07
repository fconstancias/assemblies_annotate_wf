# R port of amr_presence_absence.py -- presence/absence of each whole-cohort-distinct
# AMR gene across every real sample, resolved to the underlying contig's own Tier-1
# cluster ("genomic context") so persistence/mobility over time is visualizable.
# Run from dereplication_pilot/.

source("R/utils.R")
source("R/longitudinal_similarity.R")  # for load_sample_dates()

if (TRUE) {  # always run main logic on source() or Rscript alike -- sys.nframe()==0 broke under nested source()
  dates <- load_sample_dates()

  amr_cluster_info <- list()  # per participant: cluster_id, aro_call, group, gene_id, contig, contig_cluster

  for (participant in all_participants()) {
    membership <- load_cluster_membership(participant, "orfs")
    contig_membership <- load_cluster_membership(participant, "contigs") %>%
      select(contig_id = member_id, contig_cluster = cluster_id, contig_length = length)

    groups <- unique(membership$source_group)
    calls <- purrr::map_dfr(groups, load_gene_calls_rgi)
    gff_all <- purrr::map_dfr(groups, function(g) read_gff3(g) %>% mutate(source_group = g))

    # Which first-level clusters are AMR-defining (>=1 member, in ANY group, has an RGI
    # hit)? This must stay separate from "which specific gene copies have their own RGI
    # hit" -- presence is intentionally the LOOSER "this group has a member in the
    # AMR-defining cluster" definition (matches the Python original's docstring: a
    # cluster-mate counts as present even if RGI didn't independently flag that sample's
    # own copy). An earlier version of this script inner-joined membership straight to
    # calls, which silently switched to the stricter "has its own RGI hit" definition and
    # undercounted real presence (58593 rows/max prevalence 285 in the original Python
    # output vs 50386/268 from that bug) -- caught by re-comparing against Python's
    # already-published numbers before trusting this port.
    amr_cluster_ids <- membership %>%
      inner_join(calls, by = c("member_id" = "gene_id")) %>%
      distinct(cluster_id, .keep_all = TRUE) %>%
      select(cluster_id, aro_call)

    # resolve one representative contig (and its Tier-1 cluster) per (cluster_id, group),
    # from ALL members of that cluster in that group -- not just RGI-hit members.
    resolved <- membership %>%
      semi_join(amr_cluster_ids, by = "cluster_id") %>%
      distinct(cluster_id, source_group, .keep_all = TRUE) %>%
      left_join(amr_cluster_ids, by = "cluster_id") %>%
      rename(gene_id = member_id) %>%
      left_join(gff_all %>% select(gene_id, contig, source_group), by = c("gene_id", "source_group")) %>%
      left_join(contig_membership, by = c("contig" = "contig_id")) %>%
      mutate(participant = participant)

    amr_cluster_info[[participant]] <- resolved
  }
  amr_cluster_info <- bind_rows(amr_cluster_info)

  # Second-level (cross-participant) clustering
  second_level <- read_tsv(file.path(PILOT_DIR, "cross_participant_orfs/clu_cluster.tsv"),
                            col_names = c("amr_gene", "cluster_id"), col_types = "cc", progress = FALSE)

  presence <- amr_cluster_info %>%
    inner_join(second_level, by = "cluster_id") %>%
    group_by(amr_gene, source_group, participant) %>%
    summarise(aro_call = first(aro_call), contig_cluster = first(contig_cluster),
              contig_length = first(contig_length), .groups = "drop") %>%
    rename(group = source_group) %>%
    mutate(is_coassembly = is_coassembly(group)) %>%
    left_join(dates %>% select(group, date), by = "group") %>%
    mutate(date = ifelse(is_coassembly, NA, as.character(date)))

  write_tsv(presence %>% mutate(date = coalesce(date, "")) %>%
              select(amr_gene, aro_call, group, participant, is_coassembly, date, contig_cluster, contig_length),
            file.path(PILOT_DIR, "R_report_amr_presence_absence.tsv"))

  prevalence <- presence %>%
    group_by(amr_gene) %>%
    summarise(aro_call = first(aro_call), n_samples = n_distinct(group), .groups = "drop") %>%
    arrange(desc(n_samples))
  write_tsv(prevalence, file.path(PILOT_DIR, "R_report_amr_prevalence.tsv"))

  write_tsv(dates %>% arrange(participant, date), file.path(PILOT_DIR, "R_report_sample_calendar.tsv"))

  cat(sprintf("wrote %d presence rows, %d distinct whole-cohort AMR genes\n",
              nrow(presence), n_distinct(presence$amr_gene)))
  cat(sprintf("prevalence range: %d-%d samples\n", min(prevalence$n_samples), max(prevalence$n_samples)))
}
