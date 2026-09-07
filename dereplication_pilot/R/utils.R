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

#' contig_id -> classification ("plasmid"/"virus"/"plasmid+virus"), via MAP's own
#' contigID.map (geNomad's input contigs get renamed contig_<N> before running; this
#' file resolves that back to our own contig IDs).
load_contig_calls_genomad <- function(group) {
  mpath <- contigid_map_path(group)
  if (!file.exists(mpath)) return(tibble(contig_id = character(), call = character()))
  id_map <- read_tsv(mpath, col_names = c("genomad_id", "contig_id"), col_types = "cc", progress = FALSE) %>%
    mutate(genomad_id = str_remove(genomad_id, "^>"))

  paths <- genomad_summary_paths(group)
  plasmid_ids <- if (file.exists(paths$plasmid)) {
    read_tsv(paths$plasmid, col_types = cols(seq_name = "c", .default = "c"), progress = FALSE)$seq_name
  } else character()
  virus_ids <- if (file.exists(paths$virus)) {
    read_tsv(paths$virus, col_types = cols(seq_name = "c", .default = "c"), progress = FALSE)$seq_name
  } else character()

  calls <- bind_rows(
    tibble(genomad_id = plasmid_ids, tag = "plasmid"),
    tibble(genomad_id = virus_ids, tag = "virus")
  ) %>%
    inner_join(id_map, by = "genomad_id") %>%
    group_by(contig_id) %>%
    summarise(call = paste(sort(unique(tag)), collapse = "+"), .groups = "drop")
  calls
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
