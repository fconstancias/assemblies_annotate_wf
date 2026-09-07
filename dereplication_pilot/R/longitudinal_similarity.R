# R port of longitudinal_similarity.py -- cluster-sharing (Jaccard) between every pair
# of a participant's own single-sample timepoints, vs. real calendar days between them.
# Run from dereplication_pilot/.

source("R/utils.R")
library(lubridate)

METADATA <- "/maps/projects/hansen_ol-AUDIT/scratch/NILU/metagenomes/to_bin_sample_host_date.tsv"

load_sample_dates <- function() {
  read_tsv(METADATA, col_types = cols(.default = "c"), progress = FALSE) %>%
    mutate(group = paste0("spaS", str_remove(Sample, "^S_")),
           participant = paste0("p", Subject),
           date = dmy(Time)) %>%
    filter(!is.na(date)) %>%
    select(group, participant, date)
}

run_tier <- function(tier, dates) {
  out_rows <- list()

  for (participant in all_participants()) {
    p <- file.path(PILOT_DIR, participant, tier, "cluster_membership.tsv")
    if (!file.exists(p)) next
    membership <- load_cluster_membership(participant, tier)

    clusters_of_group <- membership %>%
      filter(!is_coassembly(source_group)) %>%
      group_by(source_group) %>%
      summarise(clusters = list(unique(cluster_id)), .groups = "drop")

    groups <- clusters_of_group$source_group
    if (length(groups) < 2) next

    for (i in seq_along(groups)) {
      for (j in seq_len(length(groups))) {
        if (j <= i) next
        gi <- groups[i]; gj <- groups[j]
        di <- dates$date[dates$group == gi]
        dj <- dates$date[dates$group == gj]
        if (length(di) == 0 || length(dj) == 0) next

        ci <- clusters_of_group$clusters[[i]]
        cj <- clusters_of_group$clusters[[j]]
        shared <- length(intersect(ci, cj))
        uni <- length(union(ci, cj))

        out_rows[[length(out_rows) + 1]] <- tibble(
          participant = participant, group_i = gi, group_j = gj,
          days_between = abs(as.integer(dj[1] - di[1])),
          shared_clusters = shared, n_clusters_i = length(ci), n_clusters_j = length(cj),
          jaccard = if (uni > 0) shared / uni else 0
        )
      }
    }
  }

  bind_rows(out_rows)
}

if (TRUE) {  # always run main logic on source() or Rscript alike -- sys.nframe()==0 broke under nested source()
  dates <- load_sample_dates()
  for (tier in c("orfs", "contigs")) {
    result <- run_tier(tier, dates)
    write_tsv(result, file.path(PILOT_DIR, sprintf("R_longitudinal_similarity_%s.tsv", tier)))
    cat(sprintf("[%s] wrote %d pairs, days %d-%d, jaccard %.3f-%.3f\n",
                tier, nrow(result), min(result$days_between), max(result$days_between),
                min(result$jaccard), max(result$jaccard)))
  }
}
