# R port of consolidate_for_report.py -- parses every participant/tier's summary.txt
# (written by R/analyze_clusters.R) into clean report_overview.tsv/report_per_group.tsv,
# plus the sensitivity-sweep comparison (report_sweep.tsv). Run from dereplication_pilot/.

source("R/utils.R")

parse_summary <- function(path) {
  text <- read_lines(path)
  total_line <- text[str_starts(text, "total sequences pooled")]
  total_seqs <- as.integer(str_match(total_line, "pooled: (\\d+)")[, 2])
  # as.numeric, not as.integer: bp totals for the largest participants exceed R's 32-bit
  # integer range (~2.1B) -- confirmed the hard way, p431's contig total_bp (4.15B)
  # silently became NA under as.integer() with only a generic "NAs introduced by
  # coercion" warning, easy to miss.
  total_bp <- as.numeric(str_match(total_line, "\\((\\d+) bp")[, 2])

  clusters_line <- text[str_starts(text, "total clusters")]
  m <- str_match(clusters_line, "total clusters: (\\d+)\\s+\\(singleton: (\\d+), multi-member: (\\d+)\\)")
  total_clusters <- as.integer(m[2]); singleton <- as.integer(m[3]); multi <- as.integer(m[4])

  cross_line <- text[str_starts(text, "cross-group clusters")]
  cross_group <- as.integer(str_match(cross_line, ": (\\d+)$")[, 2])

  redundant_line <- text[str_starts(text, "redundant sequences")]
  m2 <- str_match(redundant_line, "redundant sequences: (\\d+)/(\\d+) \\(([\\d.]+)% by count, ([\\d.]+)% by bp\\)")
  redundant_count <- as.integer(m2[2]); redundant_pct <- as.numeric(m2[4]); redundant_bp_pct <- as.numeric(m2[5])

  per_group_lines <- text[str_detect(text, "^\\s+\\S+: \\d+/\\d+ \\([\\d.]+%\\) by count")]
  per_group <- tibble()
  if (length(per_group_lines) > 0) {
    gm <- str_match(per_group_lines,
      "^\\s+(\\S+): (\\d+)/(\\d+) \\(([\\d.]+)%\\) by count, (\\d+)/(\\d+) \\(([\\d.]+)%\\) by bp")
    per_group <- tibble(group = gm[, 2], redundant_count = as.integer(gm[, 3]), total_count = as.integer(gm[, 4]),
                         redundant_pct = as.numeric(gm[, 5]), redundant_bp = as.numeric(gm[, 6]),
                         total_bp_g = as.numeric(gm[, 7]), redundant_bp_pct = as.numeric(gm[, 8]))
  }

  list(total_seqs = total_seqs, total_bp = total_bp, total_clusters = total_clusters,
       singleton = singleton, multi = multi, cross_group_clusters = cross_group,
       redundant_count = redundant_count, redundant_pct = redundant_pct,
       redundant_bp_pct = redundant_bp_pct, per_group = per_group)
}

overview_rows <- list()
per_group_rows <- list()

for (tier in c("contigs", "orfs")) {
  for (participant in all_participants()) {
    f <- file.path(PILOT_DIR, participant, tier, "summary.txt")
    if (!file.exists(f)) next
    d <- parse_summary(f)
    overview_rows[[paste(participant, tier)]] <- tibble(
      participant = participant, tier = tier, total_seqs = d$total_seqs, total_bp = d$total_bp,
      total_clusters = d$total_clusters, singleton = d$singleton, multi = d$multi,
      redundant_count = d$redundant_count, redundant_pct = d$redundant_pct, redundant_bp_pct = d$redundant_bp_pct
    )
    if (nrow(d$per_group) > 0) {
      per_group_rows[[paste(participant, tier)]] <- d$per_group %>%
        mutate(participant = participant, tier = tier, is_coassembly = is_coassembly(group)) %>%
        select(participant, tier, group, is_coassembly, redundant_count, total_count,
               redundant_pct, redundant_bp, total_bp = total_bp_g, redundant_bp_pct)
    }
  }
}

write_tsv(bind_rows(overview_rows), file.path(PILOT_DIR, "report_overview.tsv"))
write_tsv(bind_rows(per_group_rows), file.path(PILOT_DIR, "report_per_group.tsv"))

# Sensitivity sweep (p110/p97 only -- the two participants the sweep was run on)
sweep_rows <- list()
for (participant in c("p110", "p97")) {
  d0 <- parse_summary(file.path(PILOT_DIR, participant, "contigs/summary.txt"))
  sweep_rows[[length(sweep_rows) + 1]] <- tibble(participant, tier = "contigs", setting = "cov-mode 2 (original)",
                                                   redundant_pct = d0$redundant_pct, redundant_bp_pct = d0$redundant_bp_pct)
  d1 <- parse_summary(file.path(PILOT_DIR, participant, "contigs/sweep_covmode0/summary.txt"))
  sweep_rows[[length(sweep_rows) + 1]] <- tibble(participant, tier = "contigs", setting = "cov-mode 0",
                                                   redundant_pct = d1$redundant_pct, redundant_bp_pct = d1$redundant_bp_pct)
  d2 <- parse_summary(file.path(PILOT_DIR, participant, "orfs/summary.txt"))
  sweep_rows[[length(sweep_rows) + 1]] <- tibble(participant, tier = "orfs", setting = "min-seq-id 0.95 (original)",
                                                   redundant_pct = d2$redundant_pct, redundant_bp_pct = d2$redundant_bp_pct)
  d3 <- parse_summary(file.path(PILOT_DIR, participant, "orfs/sweep_minid099/summary.txt"))
  sweep_rows[[length(sweep_rows) + 1]] <- tibble(participant, tier = "orfs", setting = "min-seq-id 0.99",
                                                   redundant_pct = d3$redundant_pct, redundant_bp_pct = d3$redundant_bp_pct)
}
write_tsv(bind_rows(sweep_rows), file.path(PILOT_DIR, "report_sweep.tsv"))

cat(sprintf("wrote overview (%d rows), per_group (%d rows), sweep (%d rows)\n",
            length(overview_rows), sum(sapply(per_group_rows, nrow)), length(sweep_rows)))
