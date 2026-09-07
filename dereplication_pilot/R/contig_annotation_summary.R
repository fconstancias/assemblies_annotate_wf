# Rich per-contig-cluster annotation table for the requested annotation heatmap: length,
# gene count, AMR gene count, sample/participant prevalence, plasmid/virus classification,
# and MAG bin status (via Binette). Restricted to "interesting" contig clusters (AMR-
# carrying and/or plasmid/virus-classified -- the same set section 10 already uses,
# ~58K of the 5.48M total clusters) since a heatmap of every singleton contig would be
# meaningless. Run from dereplication_pilot/.

source("R/utils.R")

plasmid_amr <- read_tsv(file.path(PILOT_DIR, "plasmid_amr_join.tsv"),
                          col_types = cols(.default = "c")) %>%
  # plasmid_amr_join.tsv is now R-generated (canonical since the Python-to-R port), so
  # its booleans are R-native "TRUE"/"FALSE" -- NOT Python's "True"/"False". Got this
  # wrong on the first pass here despite already having hit and fixed the identical
  # mismatch in amr_plasmid_within_between.R earlier in the same port.
  mutate(is_plasmid = is_plasmid == "TRUE", is_virus = is_virus == "TRUE", has_amr = has_amr == "TRUE",
         n_members = as.integer(n_members))

interesting <- plasmid_amr %>% filter(is_plasmid | is_virus | has_amr)
cat(sprintf("interesting first-level clusters: %d of %d total\n", nrow(interesting), nrow(plasmid_amr)))

# Second-level (cross-participant) identity for each interesting first-level cluster
second_level <- read_tsv(file.path(PILOT_DIR, "cross_participant/clu_cluster.tsv"),
                          col_names = c("whole_cohort_id", "cluster_id"), col_types = "cc", progress = FALSE)
interesting <- interesting %>% left_join(second_level, by = "cluster_id")

per_group_rows <- list()

for (participant in unique(interesting$participant)) {
  clusters_here <- interesting %>% filter(participant == !!participant) %>% pull(cluster_id)
  membership <- load_cluster_membership(participant, "contigs") %>%
    filter(cluster_id %in% clusters_here)
  groups <- unique(membership$source_group)

  gff_counts <- purrr::map_dfr(groups, function(g) {
    read_gff3(g) %>% count(contig, name = "n_genes") %>% mutate(source_group = g)
  })
  amr_counts <- purrr::map_dfr(groups, function(g) {
    calls <- load_gene_calls_rgi(g)
    if (nrow(calls) == 0) return(tibble())
    gff <- read_gff3(g)
    gff %>% semi_join(calls, by = "gene_id") %>% count(contig, name = "n_amr_genes") %>% mutate(source_group = g)
  })
  bins <- purrr::map_dfr(groups, function(g) load_contig_bins(g) %>% mutate(source_group = g))

  per_group_rows[[participant]] <- membership %>%
    left_join(gff_counts, by = c("member_id" = "contig", "source_group")) %>%
    left_join(amr_counts, by = c("member_id" = "contig", "source_group")) %>%
    left_join(bins, by = c("member_id" = "contig_id", "source_group")) %>%
    mutate(n_genes = coalesce(n_genes, 0L), n_amr_genes = coalesce(n_amr_genes, 0L),
           participant = participant)
}

contig_detail <- bind_rows(per_group_rows)
write_tsv(contig_detail, file.path(PILOT_DIR, "report_contig_annotation_detail.tsv"))

# Roll up to whole-cohort distinct contig (second-level cluster) for the heatmap.
# Classification (is_plasmid/is_virus/has_amr) comes from `interesting` itself (already
# per first-level cluster_id, joined to whole_cohort_id); per-instance detail (length,
# gene counts, MAG bin status) comes from contig_detail, joined via the same cluster_id.
classification <- interesting %>%
  group_by(whole_cohort_id) %>%
  summarise(is_plasmid = any(is_plasmid), is_virus = any(is_virus), has_amr = any(has_amr), .groups = "drop")

detail_rollup <- contig_detail %>%
  left_join(interesting %>% select(cluster_id, whole_cohort_id), by = "cluster_id") %>%
  group_by(whole_cohort_id) %>%
  summarise(
    n_samples = n_distinct(source_group),
    n_participants = n_distinct(participant),
    length = max(length, na.rm = TRUE),
    n_genes = max(n_genes, na.rm = TRUE),
    n_amr_genes = max(n_amr_genes, na.rm = TRUE),
    n_good_mag = sum(is_good_mag, na.rm = TRUE),
    n_binned = sum(!is.na(bin_name)),
    n_instances = n(),
    .groups = "drop"
  ) %>%
  mutate(
    # MAG binning is done independently per sample/assembly, so it can genuinely
    # succeed in some participants' assemblies and fail in others for the exact same
    # whole-cohort contig -- confirmed on real data (mh_p550_000000053967: a plasmid
    # shared across 11 participants, good-MAG in 4 of them, unbinned in the rest). A
    # flat any()-derived "good MAG"/"unbinned" label collapsed that real heterogeneity
    # into a misleadingly definitive-looking single value; report the actual fraction
    # instead.
    pct_good_mag = round(100 * n_good_mag / n_instances, 1),
    mag_status = case_when(
      n_binned == 0 ~ "never binned",
      pct_good_mag >= 80 ~ "consistently good MAG",
      n_good_mag > 0 ~ "sometimes good MAG",
      TRUE ~ "binned, never good quality"
    )
  )

summary <- classification %>% left_join(detail_rollup, by = "whole_cohort_id")

write_tsv(summary, file.path(PILOT_DIR, "report_contig_annotation_summary.tsv"))
cat(sprintf("wrote %d whole-cohort distinct interesting contigs\n", nrow(summary)))
print(summary %>% count(mag_status))
