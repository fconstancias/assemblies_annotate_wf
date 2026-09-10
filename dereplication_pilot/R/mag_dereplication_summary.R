# Real dereplicated MAG count + cross-participant sharing, directly from dRep's own
# per-strand species-level clustering (dRep/data_tables/Cdb.csv, already run on
# Binette's own bins -- see utils.R's load_bin_taxonomy() for how this was found and
# reconciled to Binette's own bin_name scheme). NOT cross-strand: co-assembly and
# single-assembly dRep runs are two separate, independent clusterings that were never
# compared to each other, so the same organism recovered in both a participant's
# co-assembly AND their own single-sample assemblies counts as two separate
# dereplicated MAGs here -- a real limitation of the underlying data, flagged
# explicitly in the report rather than silently reconciled.
# Run from dereplication_pilot/.

source("R/utils.R")
source("R/longitudinal_similarity.R")  # for load_sample_dates()

group_to_participant_map <- function() {
  purrr::map_dfr(all_participants(), function(p) {
    load_cluster_membership(p, "contigs") %>% distinct(source_group) %>% mutate(participant = p)
  }) %>% distinct(source_group, .keep_all = TRUE)
}

load_strand_drep <- function(is_coa, strand_label) {
  base <- dirname(if (is_coa) COASSEMBLY_BINETTE else SINGLE_BINETTE)
  p <- file.path(base, "dRep/data_tables/Cdb.csv")
  read_csv(p, col_types = cols(genome = "c", secondary_cluster = "c", .default = "c"), progress = FALSE) %>%
    transmute(
      is_coassembly = is_coa,
      strand = strand_label,
      dereplicated_mag_id = paste0(strand_label, "_", secondary_cluster),
      group = str_remove(genome, "_binette_bin\\d+\\.fa$"),
      bin_name = str_remove(genome, "\\.fa$")
    )
}

load_strand_winners <- function(is_coa, strand_label) {
  base <- dirname(if (is_coa) COASSEMBLY_BINETTE else SINGLE_BINETTE)
  p <- file.path(base, "dRep/data_tables/Wdb.csv")
  read_csv(p, col_types = cols(genome = "c", cluster = "c", .default = "c"), progress = FALSE) %>%
    transmute(dereplicated_mag_id = paste0(strand_label, "_", cluster),
              winner_bin_name = str_remove(genome, "\\.fa$"))
}

#' Real total bin counts independent of dRep's own quality pre-filter (dRep only ever
#' sees a subset -- 9,094 of 9,996 real Binette final_bins -- so "total bins" and
#' "bins entering dRep" are two genuinely different numbers, not the same one twice).
count_real_bins <- function() {
  quality_files <- c(
    Sys.glob(file.path(dirname(COASSEMBLY_BINETTE), "binette/*/final_bins_quality_reports.tsv")),
    Sys.glob(file.path(dirname(SINGLE_BINETTE), "binette/*/final_bins_quality_reports.tsv"))
  )
  purrr::map_dfr(quality_files, function(p) {
    read_tsv(p, col_types = cols(completeness = "d", contamination = "d", .default = "c"), progress = FALSE)
  }) %>%
    summarise(total_bins = n(), good_mags = sum(completeness >= 50 & contamination < 10, na.rm = TRUE))
}

if (TRUE) {
  bin_totals <- count_real_bins()
  cat(sprintf("real total bins (both strands): %d, good-quality (>=50%% comp, <10%% contam): %d\n",
              bin_totals$total_bins, bin_totals$good_mags))
  write_tsv(bin_totals, file.path(PILOT_DIR, "R_report_mag_totals.tsv"))

  g2p <- group_to_participant_map()
  dates <- load_sample_dates()

  mags <- bind_rows(
    load_strand_drep(TRUE, "coassembly"),
    load_strand_drep(FALSE, "single")
  ) %>%
    left_join(g2p, by = c("group" = "source_group")) %>%
    left_join(dates, by = c("group", "participant"))

  cat(sprintf("total bins entering dRep (both strands): %d\n", nrow(mags)))
  cat(sprintf("bins with no participant match: %d\n", sum(is.na(mags$participant))))

  # Canonical taxonomy per dereplicated MAG: dRep's own cluster "winner" (Wdb.csv, the
  # highest-scoring genome) -- one real, specific genome's own GTDB-Tk call, not a vote
  # across the cluster's members (which should agree anyway, since dRep only clusters
  # genomes within its own ANI threshold, but the winner is the principled choice).
  winners <- bind_rows(load_strand_winners(TRUE, "coassembly"), load_strand_winners(FALSE, "single"))
  all_taxonomy <- bind_rows(load_bin_taxonomy(TRUE), load_bin_taxonomy(FALSE)) %>%
    distinct(bin_name, .keep_all = TRUE)
  winner_tax <- winners %>%
    left_join(all_taxonomy, by = c("winner_bin_name" = "bin_name")) %>%
    select(dereplicated_mag_id, gtdb_taxonomy)

  per_mag <- mags %>%
    group_by(strand, dereplicated_mag_id) %>%
    summarise(n_bins = n(), n_participants = n_distinct(participant), .groups = "drop") %>%
    left_join(winner_tax, by = "dereplicated_mag_id") %>%
    mutate(phylum = str_extract(gtdb_taxonomy, "(?<=p__)[^;]+"))

  cat(sprintf("dereplicated MAGs (species-level genomes): %d\n", nrow(per_mag)))
  cat(sprintf("shared across >=2 participants: %d (%.1f%%)\n",
              sum(per_mag$n_participants > 1), 100 * mean(per_mag$n_participants > 1)))

  write_tsv(per_mag, file.path(PILOT_DIR, "R_report_mag_dereplication_summary.tsv"))

  per_participant <- mags %>%
    filter(!is.na(participant)) %>%
    distinct(participant, dereplicated_mag_id) %>%
    count(participant, name = "n_dereplicated_mags")
  write_tsv(per_participant, file.path(PILOT_DIR, "R_report_mag_per_participant.tsv"))
  cat(sprintf("wrote per-participant dereplicated MAG counts for %d participants\n", nrow(per_participant)))

  # Per-instance table (one row per bin) for a presence/across-time view -- only
  # single-sample instances carry a real date (co-assemblies pool an entire
  # participant's timeline into one assembly, so they have no single date).
  instances <- mags %>%
    select(dereplicated_mag_id, strand, group, participant, is_coassembly, date) %>%
    left_join(per_mag %>% select(dereplicated_mag_id, gtdb_taxonomy, phylum), by = "dereplicated_mag_id")
  write_tsv(instances, file.path(PILOT_DIR, "R_report_mag_instances.tsv"))
  cat(sprintf("wrote %d per-bin-instance rows (%d with a real sample date)\n",
              nrow(instances), sum(!is.na(instances$date))))
}
