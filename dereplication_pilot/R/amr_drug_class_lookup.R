# Drug Class / Resistance Mechanism / AMR Gene Family per aro_call (RGI's Best_Hit_ARO) --
# real columns RGI already writes (confirmed in the raw .txt header) but load_gene_calls_rgi()
# never kept, since nothing needed them until now. These are a property of the CARD/ARO
# term itself, not of any one instance, so one representative RGI row per distinct
# aro_call (from any group that carries it) is enough -- no need to touch every instance
# or rebuild the master table. Scoped to only the groups already known to carry >=1 AMR
# hit (from master_orf_table_amr.tsv.gz's own `sample` column), not all 296 real groups.
# Run from dereplication_pilot/.

source("R/utils.R")

master <- read_tsv(file.path(PILOT_DIR, "master_orf_table_amr.tsv.gz"), col_types = cols(.default = "c"))
groups <- unique(master$sample)
cat(sprintf("scanning %d AMR-carrying groups for drug class info...\n", length(groups)))

lookup <- purrr::map_dfr(groups, function(g) {
  read_rgi(g) %>% distinct(Best_Hit_ARO, `Drug Class`, `Resistance Mechanism`, `AMR Gene Family`)
}) %>%
  distinct(Best_Hit_ARO, .keep_all = TRUE) %>%
  transmute(aro_call = Best_Hit_ARO, drug_class = `Drug Class`,
            resistance_mechanism = `Resistance Mechanism`, amr_gene_family = `AMR Gene Family`)

cat(sprintf("resolved drug class for %d of %d distinct aro_call values\n",
            sum(!is.na(lookup$drug_class) & lookup$drug_class != ""), n_distinct(master$aro_call)))

write_tsv(lookup, file.path(PILOT_DIR, "R_report_amr_drug_class_lookup.tsv"))
