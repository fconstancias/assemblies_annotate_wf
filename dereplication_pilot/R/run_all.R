# Driver: runs the full R port of the dereplication analysis pipeline, in dependency
# order, writing R_report_*.tsv outputs (the "R_" prefix keeps them alongside the
# original Python-derived report_*.tsv during side-by-side validation; drop the prefix
# once the Rmd is switched over and the Python originals are retired).
#
# Usage: Rscript R/run_all.R   (run from dereplication_pilot/)

stopifnot("run this from dereplication_pilot/ (Rscript R/run_all.R)" = file.exists("R/utils.R"))

cat("=== 1/9: consolidate_for_report ===\n"); source("R/consolidate_for_report.R")
cat("=== 2/9: validate_propagation ===\n"); source("R/validate_propagation.R")
cat("=== 3/9: amr_and_coassembly_reports ===\n"); source("R/amr_and_coassembly_reports.R")
cat("=== 4/9: longitudinal_similarity ===\n"); source("R/longitudinal_similarity.R")
cat("=== 5/9: plasmid_amr_join ===\n"); source("R/plasmid_amr_join.R")
cat("=== 6/9: cross_participant_analysis ===\n"); source("R/cross_participant_analysis.R")
cat("=== 7/9: amr_within_between ===\n"); source("R/amr_within_between.R")
cat("=== 8/9: amr_plasmid_within_between ===\n"); source("R/amr_plasmid_within_between.R")
cat("=== 9/9: amr_presence_absence ===\n"); source("R/amr_presence_absence.R")
cat("=== all done ===\n")
