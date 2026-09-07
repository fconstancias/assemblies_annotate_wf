# Shared loaders for the dereplication pipeline -- R port of validate_propagation.py's
# path builders and RGI/geNomad joins. Ported (not just translated) using vectorized
# readr/dplyr joins rather than the Python originals' per-row loops, since several of
# these gff3 files run 300K+ lines and a naive R for-loop over that would be far slower
# than the equivalent Python -- vectorized joins are what makes R competitive here.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
})

REPO <- "/maps/projects/hansen_ol-AUDIT/scratch/NILU/metagenomes/assembly_annotation_wf"
COASSEMBLY_GENE_EXPORT <- file.path(REPO, "coassembly_production/gene_export")
SINGLE_GENE_EXPORT <- file.path(REPO, "single_assembly_production/gene_export")
COASSEMBLY_FUNCSCAN <- file.path(REPO, "coassembly_production/funcscan_run/results")
SINGLE_FUNCSCAN <- file.path(REPO, "single_assembly_production/funcscan_run/results")
COASSEMBLY_MAP <- file.path(REPO, "coassembly_production/map_run/results")
SINGLE_MAP <- file.path(REPO, "single_assembly_production/map_run/results")
PILOT_DIR <- file.path(REPO, "dereplication_pilot")

is_coassembly <- function(group) str_starts(group, "mh_p")

gff3_path <- function(group) {
  base <- if (is_coassembly(group)) COASSEMBLY_GENE_EXPORT else SINGLE_GENE_EXPORT
  file.path(base, paste0(group, ".gff3"))
}

rgi_path <- function(group) {
  base <- if (is_coassembly(group)) COASSEMBLY_FUNCSCAN else SINGLE_FUNCSCAN
  file.path(base, "arg/rgi", group, paste0(group, ".txt"))
}

contigid_map_path <- function(group) {
  base <- if (is_coassembly(group)) COASSEMBLY_MAP else SINGLE_MAP
  file.path(base, group, "preprocessing", paste0(group, "_contigID.map"))
}

genomad_summary_paths <- function(group) {
  base <- if (is_coassembly(group)) COASSEMBLY_MAP else SINGLE_MAP
  dir <- file.path(base, group, "prediction/genomad")
  list(
    plasmid = file.path(dir, paste0(group, "_5kb_contigs_plasmid_summary.tsv")),
    virus   = file.path(dir, paste0(group, "_5kb_contigs_virus_summary.tsv"))
  )
}

COASSEMBLY_BINETTE <- "/maps/projects/hansen_ol-AUDIT/scratch/NILU/metagenomes/spa_coassembly_all/results_coassembly_all/binette"
SINGLE_BINETTE <- "/maps/projects/hansen_ol-AUDIT/scratch/NILU/metagenomes/spa_single_all/results_spa_single_all/binette"

binette_paths <- function(group) {
  base <- if (is_coassembly(group)) COASSEMBLY_BINETTE else SINGLE_BINETTE
  dir <- file.path(base, group)
  list(contig_to_bin = file.path(dir, "final_contig_to_bin.tsv"),
       bin_quality = file.path(dir, "final_bins_quality_reports.tsv"))
}

#' Real GTDB-Tk taxonomy + dRep cross-sample species-level cluster for every Binette
#' bin -- NOT a new pass to run: dRep and GTDB-Tk already exist as siblings of binette/
#' itself (spa_coassembly_all/results_coassembly_all/{dRep,gtdbtk_classify}/ and the
#' spa_single_all equivalent), run ONCE per assembly strand directly on Binette's own
#' final_bins/*.fa (via binette_renamed_bins/, genome names "<group>_binette_binN.fa" --
#' this pipeline's own group-prefix convention, already matching load_contig_bins()'s
#' bin_name once ".fa" is stripped). This is a DIFFERENT, unrelated pass from the
#' top-level metagenomes/dRep + metagenomes/gtdbtk_classify (those use concoct-manual/
#' semibin2 bin identities -- a separate, earlier binning strategy that mixes both
#' binners across BOTH assembly strands, confirmed not reconcilable to Binette's own
#' bins by name). GTDB-Tk classified essentially every individual bin here (8286
#' summary rows for 8285 single-assembly final_bins), not just dRep cluster
#' representatives, so taxonomy is a direct per-bin join -- dRep's secondary_cluster is
#' exposed alongside it since it directly answers "is this bin the same organism as
#' that bin from a different sample/participant", the cross-sample MAG identity link
#' load_contig_bins() alone can't provide (Binette itself bins each assembly group
#' independently, never across samples). Cached per strand (one shared file pair, not
#' per-group) since build_one_group() calls load_contig_bins() ~296 times.
.bin_taxonomy_cache <- new.env()
load_bin_taxonomy <- function(is_coa) {
  key <- as.character(is_coa)
  if (!is.null(.bin_taxonomy_cache[[key]])) return(.bin_taxonomy_cache[[key]])

  base <- dirname(if (is_coa) COASSEMBLY_BINETTE else SINGLE_BINETTE)
  gtdb_bac <- file.path(base, "gtdbtk_classify/gtdbtk.bac120.summary.tsv")
  gtdb_ar  <- file.path(base, "gtdbtk_classify/gtdbtk.ar53.summary.tsv")
  drep_cdb <- file.path(base, "dRep/data_tables/Cdb.csv")

  read_gtdb <- function(p) {
    if (!file.exists(p)) return(tibble(bin_name = character(), gtdb_taxonomy = character()))
    read_tsv(p, col_types = cols(user_genome = "c", classification = "c", .default = "c"), progress = FALSE) %>%
      transmute(bin_name = str_remove(user_genome, "\\.fa$"), gtdb_taxonomy = classification)
  }
  gtdb <- bind_rows(read_gtdb(gtdb_bac), read_gtdb(gtdb_ar))

  drep <- if (file.exists(drep_cdb)) {
    read_csv(drep_cdb, col_types = cols(genome = "c", secondary_cluster = "c", .default = "c"), progress = FALSE) %>%
      transmute(bin_name = str_remove(genome, "\\.fa$"), drep_secondary_cluster = secondary_cluster)
  } else tibble(bin_name = character(), drep_secondary_cluster = character())

  result <- full_join(gtdb, drep, by = "bin_name")
  .bin_taxonomy_cache[[key]] <- result
  result
}

#' contig_id -> (bin_name, completeness, contamination, is_good_mag) for one group.
#' "good MAG" follows the common completeness>=50 / contamination<10 threshold.
#' bin_name is prefixed with the group ("mh_p110_binette_bin1", not just "binette_bin1")
#' -- Binette numbers bins independently PER GROUP, so the raw name alone is genuinely
#' ambiguous read in isolation (two different rows both saying "binette_bin1" could be
#' completely unrelated bins from different samples) even though every actual JOIN in
#' this codebase was already safe (always keyed on contig_id + source_group, never on
#' bin_name alone) -- a real output-labeling gap, not a join bug, caught directly.
load_contig_bins <- function(group) {
  paths <- binette_paths(group)
  empty_tax <- tibble(gtdb_taxonomy = character(), drep_secondary_cluster = character())
  if (!file.exists(paths$contig_to_bin)) {
    return(tibble(contig_id = character(), bin_name = character(),
                   completeness = double(), contamination = double(), is_good_mag = logical()) %>%
             bind_cols(empty_tax))
  }
  c2b <- read_tsv(paths$contig_to_bin, col_names = c("contig_id", "bin_name"), col_types = "cc", progress = FALSE) %>%
    mutate(bin_name = paste0(group, "_", bin_name))
  taxonomy <- load_bin_taxonomy(is_coassembly(group))
  if (!file.exists(paths$bin_quality)) {
    return(c2b %>% mutate(completeness = NA_real_, contamination = NA_real_, is_good_mag = NA) %>%
             left_join(taxonomy, by = "bin_name"))
  }
  bq <- read_tsv(paths$bin_quality, col_types = cols(name = "c", completeness = "d", contamination = "d", .default = "c"), progress = FALSE) %>%
    mutate(name = paste0(group, "_", name))
  c2b %>%
    left_join(bq %>% select(bin_name = name, completeness, contamination), by = "bin_name") %>%
    mutate(is_good_mag = completeness >= 50 & contamination < 10) %>%
    left_join(taxonomy, by = "bin_name")
}

#' Read one group's gff3 as a tidy (contig, start, stop, gene_id) table.
read_gff3 <- function(group) {
  p <- gff3_path(group)
  if (!file.exists(p)) return(tibble(contig = character(), start = integer(), stop = integer(), gene_id = character()))
  read_tsv(p, comment = "#", col_names = c("contig", "source", "type", "start", "stop", "score", "strand", "frame", "attrs"),
            col_types = cols(contig = "c", start = "i", stop = "i", attrs = "c", .default = "c"),
            progress = FALSE) %>%
    mutate(gene_id = str_match(attrs, "ID=([^;]+)")[, 2]) %>%
    filter(!is.na(gene_id)) %>%
    select(contig, start, stop, gene_id)
}

#' Read one group's RGI hit table (already tab-separated with a real header).
read_rgi <- function(group) {
  p <- rgi_path(group)
  if (!file.exists(p)) return(tibble())
  read_tsv(p, col_types = cols(Contig = "c", Start = "i", Stop = "i", Best_Hit_ARO = "c", .default = "c"), progress = FALSE)
}

#' gene_id -> Best_Hit_ARO for one group, joined via (contig,start,stop) coordinates --
#' RGI's own ORF_ID uses prodigal-style numbering, not our gene_export ID scheme, so
#' coordinates are the only common key (confirmed against both real formats directly
#' when this join was first built in Python).
load_gene_calls_rgi <- function(group) {
  rgi <- read_rgi(group)
  if (nrow(rgi) == 0) return(tibble(gene_id = character(), aro_call = character()))
  gff <- read_gff3(group)
  rgi %>%
    inner_join(gff, by = c("Contig" = "contig", "Start" = "start", "Stop" = "stop")) %>%
    distinct(gene_id, .keep_all = TRUE) %>%
    transmute(gene_id, aro_call = Best_Hit_ARO)
}

#' contig_id -> classification ("plasmid"/"virus"/"plasmid+virus") PLUS geNomad's own
#' real, continuous scores (plasmid_score/plasmid_fdr, virus_score/virus_fdr) -- these
#' are score-based, not a hard binary truth, and appearing in the summary file at all
#' already reflects geNomad's own internal threshold; the earlier version of this
#' function only recorded that appearance as a plain TRUE/FALSE, discarding exactly the
#' confidence information needed to tell "confidently plasmid" from "barely crossed the
#' threshold" apart -- caught by a real, ambiguous example (a geNomad-plasmid-classified,
#' Binette-binned contig carrying a canonical *chromosomal* AMR gene call, soxR). Via
#' MAP's own contigID.map (geNomad's input contigs get renamed contig_<N> before
#' running; this file resolves that back to our own contig IDs).
load_contig_calls_genomad <- function(group) {
  empty <- tibble(contig_id = character(), call = character(),
                    plasmid_score = double(), plasmid_fdr = double(),
                    virus_score = double(), virus_fdr = double())
  mpath <- contigid_map_path(group)
  if (!file.exists(mpath)) return(empty)
  id_map <- read_tsv(mpath, col_names = c("genomad_id", "contig_id"), col_types = "cc", progress = FALSE) %>%
    mutate(genomad_id = str_remove(genomad_id, "^>"))

  paths <- genomad_summary_paths(group)
  plasmid_tbl <- if (file.exists(paths$plasmid)) {
    read_tsv(paths$plasmid, col_types = cols(seq_name = "c", plasmid_score = "d", fdr = "d", .default = "c"), progress = FALSE) %>%
      select(genomad_id = seq_name, plasmid_score, plasmid_fdr = fdr)
  } else tibble(genomad_id = character(), plasmid_score = double(), plasmid_fdr = double())
  virus_tbl <- if (file.exists(paths$virus)) {
    read_tsv(paths$virus, col_types = cols(seq_name = "c", virus_score = "d", fdr = "d", .default = "c"), progress = FALSE) %>%
      select(genomad_id = seq_name, virus_score, virus_fdr = fdr)
  } else tibble(genomad_id = character(), virus_score = double(), virus_fdr = double())

  full_join(plasmid_tbl, virus_tbl, by = "genomad_id") %>%
    inner_join(id_map, by = "genomad_id") %>%
    mutate(call = case_when(
      !is.na(plasmid_score) & !is.na(virus_score) ~ "plasmid+virus",
      !is.na(plasmid_score) ~ "plasmid",
      !is.na(virus_score) ~ "virus",
      TRUE ~ NA_character_
    )) %>%
    select(contig_id, call, plasmid_score, plasmid_fdr, virus_score, virus_fdr)
}

#' Read a participant's real Tier-1 or Tier-2 cluster_membership.tsv (already produced
#' by the mmseqs clustering + analyze_clusters.py pass -- these are large real data
#' products on disk, not regenerated here).
load_cluster_membership <- function(participant, tier) {
  p <- file.path(PILOT_DIR, participant, tier, "cluster_membership.tsv")
  read_tsv(p, col_types = cols(.default = "c"), progress = FALSE)
}

all_participants <- function() {
  # "p" + digits, and must be a real directory -- a plain "p*" glob also matched
  # plasmid_amr_join.log (a leftover file from the Python pipeline) the hard way.
  candidates <- Sys.glob(file.path(PILOT_DIR, "p*"))
  candidates <- candidates[dir.exists(candidates)]
  sort(basename(candidates)[str_detect(basename(candidates), "^p[0-9]+$")])
}
