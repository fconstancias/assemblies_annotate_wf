# Contig-row presence/absence, analogous to amr_presence_absence.R but pivoted to
# CONTIGS as the row instead of genes -- one row per whole-cohort-distinct "interesting"
# contig (AMR-carrying and/or plasmid/virus-classified), tracked across every real
# sample, with the section-14-style annotations (length, participant/sample counts, AMR
# gene count, MAG status) available as row side-bars.
#
# Cheap by design: reuses report_contig_annotation_detail.tsv (per-instance data,
# already computed by contig_annotation_summary.R -- no gff3/RGI re-parsing needed) and
# just joins in the whole-cohort identity + real sample dates.
# Run from dereplication_pilot/.

source("R/utils.R")
source("R/longitudinal_similarity.R")  # for load_sample_dates()

detail <- read_tsv(file.path(PILOT_DIR, "report_contig_annotation_detail.tsv"),
                    col_types = cols(.default = "c")) %>%
  mutate(length = as.numeric(length), n_genes = as.integer(n_genes), n_amr_genes = as.integer(n_amr_genes))

second_level <- read_tsv(file.path(PILOT_DIR, "cross_participant/clu_cluster.tsv"),
                          col_names = c("whole_cohort_id", "cluster_id"), col_types = "cc", progress = FALSE)

dates <- load_sample_dates()

presence <- detail %>%
  inner_join(second_level, by = "cluster_id") %>%
  transmute(whole_cohort_id, group = source_group, participant,
            is_coassembly = is_coassembly(source_group),
            local_length = length, local_n_amr_genes = n_amr_genes) %>%
  left_join(dates %>% select(group, date), by = "group") %>%
  mutate(date = ifelse(is_coassembly, NA, as.character(date)))

write_tsv(presence %>% mutate(date = coalesce(date, "")),
          file.path(PILOT_DIR, "report_contig_presence_absence.tsv"))

cat(sprintf("wrote %d contig presence rows, %d distinct whole-cohort contigs\n",
            nrow(presence), n_distinct(presence$whole_cohort_id)))
