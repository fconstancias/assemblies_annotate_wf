# Master table: every ORF, every sample (42.4M rows across all 296 assemblies), with
# every annotation dimension already joined -- contig, contig length, geNomad
# plasmid/virus call, Binette MAG bin (name + completeness/contamination/is_good_mag),
# this ORF's own RGI AMR call if it has one -- AND the dereplication cluster identity
# for BOTH the ORF and its contig, at both levels:
#   - orf_cluster_id / contig_cluster_id: first-level (within-participant) cluster,
#     i.e. which of that participant's own mmseqs clusters this instance belongs to.
#   - orf_whole_cohort_id / contig_whole_cohort_id: second-level (cross-participant)
#     identity -- the real "same gene/contig as which one in every other sample"
#     answer, via cross_participant_orfs/ and cross_participant/'s own clustering.
# Without these, this table would just be a per-sample annotation dump with no link to
# the cross-sample identity model the rest of this pipeline is built around -- added
# after a real gap was caught before the first (incomplete) build finished.
#
# Written incrementally (one real group at a time, ~296 iterations) to stay memory-safe
# at this scale. Output is gzipped TSV (no parquet library installed in this env); still
# large (expect several GB), so NOT tracked in git -- disk + S3 only.
#
# Run from dereplication_pilot/, as a SLURM job (R/run_master_orf_table.sh).

source("R/utils.R")
source("R/longitudinal_similarity.R")  # for load_sample_dates()

all_real_groups <- function() {
  ca <- list.files(COASSEMBLY_GENE_EXPORT, pattern = "\\.gff3$") %>% str_remove("\\.gff3$")
  sa <- list.files(SINGLE_GENE_EXPORT, pattern = "\\.gff3$") %>% str_remove("\\.gff3$")
  c(ca, sa)
}

group_to_participant_map <- function() {
  purrr::map_dfr(all_participants(), function(p) {
    load_cluster_membership(p, "contigs") %>% distinct(source_group) %>% mutate(participant = p)
  }) %>% distinct(source_group, .keep_all = TRUE)
}

#' Real per-contig lengths + first-level cluster_id, reused from the dereplication
#' pilot's own already-computed contigs/cluster_membership.tsv (member_id -> length,
#' cluster_id). Only covers contigs >=1kb (the pilot's own pooling filter); real
#' gene-bearing contigs are virtually always above that.
load_all_contig_info <- function() {
  purrr::map_dfr(all_participants(), function(p) {
    load_cluster_membership(p, "contigs") %>%
      select(contig_id = member_id, contig_length = length, contig_cluster_id = cluster_id)
  }) %>% distinct(contig_id, .keep_all = TRUE)
}

#' orf_id -> first-level cluster_id, same source pattern as contig info above, but from
#' the orfs-tier cluster_membership.tsv (every real ORF that made it into a >=1
#' sequence dereplication cluster -- i.e. essentially all of them, no length filter
#' applies at the ORF tier the way it does for contigs).
load_all_orf_cluster_ids <- function() {
  purrr::map_dfr(all_participants(), function(p) {
    load_cluster_membership(p, "orfs") %>% select(orf_id = member_id, orf_cluster_id = cluster_id)
  }) %>% distinct(orf_id, .keep_all = TRUE)
}

build_one_group <- function(group, participant, contig_info, orf_cluster_ids,
                             contig_second_level, orf_second_level, dates, amr_only = TRUE) {
  rgi <- load_gene_calls_rgi(group)
  # AMR-only first pass (per explicit request: validate the full schema at a much
  # smaller, faster scale before committing to all 42.4M ORFs): skip groups with zero
  # RGI hits entirely, and join gff3 for only this group's AMR-hit genes rather than
  # every gene -- the expensive full-gff3 read/join is exactly what expanding later
  # (amr_only = FALSE) restores.
  if (amr_only && nrow(rgi) == 0) return(invisible(NULL))

  gff <- read_gff3(group)
  if (nrow(gff) == 0) return(invisible(NULL))
  if (amr_only) gff <- gff %>% semi_join(rgi, by = "gene_id")
  if (nrow(gff) == 0) return(invisible(NULL))

  genomad <- load_contig_calls_genomad(group)
  bins <- load_contig_bins(group)
  sample_date <- if (!is_coassembly(group) && group %in% dates$group) {
    as.character(dates$date[dates$group == group][1])
  } else {
    NA_character_
  }

  gff %>%
    left_join(rgi, by = "gene_id") %>%
    left_join(genomad, by = c("contig" = "contig_id")) %>%
    left_join(bins, by = c("contig" = "contig_id")) %>%
    left_join(contig_info, by = c("contig" = "contig_id")) %>%
    left_join(orf_cluster_ids, by = c("gene_id" = "orf_id")) %>%
    left_join(contig_second_level, by = "contig_cluster_id") %>%
    left_join(orf_second_level, by = "orf_cluster_id") %>%
    transmute(
      sample = group, participant = participant, is_coassembly = is_coassembly(group), date = sample_date,
      orf_id = gene_id, orf_cluster_id, orf_whole_cohort_id,
      contig_id = contig, contig_length, contig_cluster_id, contig_whole_cohort_id,
      is_plasmid = str_detect(coalesce(call, ""), "plasmid"),
      is_virus = str_detect(coalesce(call, ""), "virus"),
      plasmid_score, plasmid_fdr, virus_score, virus_fdr,
      mag_bin = bin_name, mag_completeness = completeness, mag_contamination = contamination,
      is_good_mag = is_good_mag, mag_gtdb_taxonomy = gtdb_taxonomy,
      mag_drep_secondary_cluster = drep_secondary_cluster, aro_call = aro_call
    )
}

if (TRUE) {
  g2p <- group_to_participant_map()
  dates <- load_sample_dates()

  cat("loading real per-contig lengths + first-level cluster IDs...\n")
  contig_info <- load_all_contig_info()
  cat(sprintf("  %s distinct contigs\n", format(nrow(contig_info), big.mark = ",")))

  cat("loading first-level ORF cluster IDs...\n")
  orf_cluster_ids <- load_all_orf_cluster_ids()
  cat(sprintf("  %s distinct ORFs\n", format(nrow(orf_cluster_ids), big.mark = ",")))

  cat("loading second-level (cross-participant) contig identity...\n")
  contig_second_level <- read_tsv(file.path(PILOT_DIR, "cross_participant/clu_cluster.tsv"),
                                    col_names = c("contig_whole_cohort_id", "contig_cluster_id"),
                                    col_types = "cc", progress = FALSE)
  cat(sprintf("  %s rows\n", format(nrow(contig_second_level), big.mark = ",")))

  cat("loading second-level (cross-participant) ORF identity...\n")
  orf_second_level <- read_tsv(file.path(PILOT_DIR, "cross_participant_orfs/clu_cluster.tsv"),
                                 col_names = c("orf_whole_cohort_id", "orf_cluster_id"),
                                 col_types = "cc", progress = FALSE)
  cat(sprintf("  %s rows\n", format(nrow(orf_second_level), big.mark = ",")))

  # AMR_ONLY: first pass per explicit request ("let's start with AMR genes first then
  # we will expand from there") -- validates the full schema (cluster IDs at both
  # levels, geNomad, Binette, dates) at a much smaller/faster scale. Flip to FALSE and
  # rerun for the full 42.4M-row table once this is confirmed correct.
  AMR_ONLY <- TRUE
  out_path <- file.path(PILOT_DIR, if (AMR_ONLY) "master_orf_table_amr.tsv.gz" else "master_orf_table.tsv.gz")
  if (file.exists(out_path)) file.remove(out_path)

  groups <- all_real_groups()
  cat(sprintf("processing %d real groups (AMR_ONLY=%s)\n", length(groups), AMR_ONLY))

  total_rows <- 0
  for (i in seq_along(groups)) {
    g <- groups[i]
    participant <- g2p$participant[g2p$source_group == g]
    if (length(participant) == 0) participant <- NA_character_
    rows <- build_one_group(g, participant, contig_info, orf_cluster_ids,
                             contig_second_level, orf_second_level, dates, amr_only = AMR_ONLY)
    if (!is.null(rows) && nrow(rows) > 0) {
      write_tsv(rows, out_path, append = file.exists(out_path))
      total_rows <- total_rows + nrow(rows)
    }
    if (i %% 20 == 0) cat(sprintf("[%d/%d] %s -- %s total rows so far\n", i, length(groups), g, format(total_rows, big.mark = ",")))
  }
  cat(sprintf("DONE: %s total ORF rows written to %s\n", format(total_rows, big.mark = ","), out_path))
}
