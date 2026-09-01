# Single-sample production run — all 277 real single-sample assembly groups

Companion to `../coassembly_production/` — same architecture, same config, same fixes,
scaled to the single-sample strand of the project (`spa_single_all/`, 277 groups across
all 25 participants — see that directory's `SINGLE_ASSEMBLY.md`) instead of the 19
per-participant co-assemblies.

## Scope

All 277 groups in `spa_single_all_anvio_kegg/group_list.txt` (real contigs-dbs already
copied there for the anvi'o-KEGG run — reused directly for gene export here too, read-only,
doesn't conflict with that array's concurrent KOfam-annotation writes to the same files).
Real per-group scale: 132K-321K genes (smaller than co-assembly's 230K-696K, but 277 groups
vs 19 — ~14x the total scope).

## What's running

- **`gene_export/`**: gff3+faa per group, `export_gene_calls_array.sh` (SLURM array,
  `%40` throttle — gene export is lightweight, no HMM search, safe to run more concurrently
  than the heavier KOfam/MAP arrays). Completed for all 277 groups, verified clean (CDS
  count == FAA sequence count for every group, no mismatches).
- **`map_run/`**: MAP, full scope (`skip_sanntis=false`, `skip_gecco=true`,
  `skip_antismash=true` — BGC is SanntiS's job, same decision as the co-assembly run).
  `nextflow.config` is a direct copy of `../coassembly_production/map_run/nextflow.config`
  — reuses every fix already proven there (boolean-CLI regression, VFDB/geNomad DB-download
  bugs, ICEFINDER2_LITE disabled, SANNTIS/INTERPROSCAN memory+time overrides) rather than
  re-discovering them at this scale too. One deliberate change: `executor.queueSize` trimmed
  from 200 to 100, since this run shares the cluster with the still-running single-sample
  anvi'o-KEGG array at the same time, on top of its own 277-group scope.
- **anvi'o-KEGG** (`../../spa_single_all_anvio_kegg/`, outside this repo): already running
  independently before this MAP run started, same architecture as the co-assembly one
  (own contigs-db copies, never the shared upstream). Not duplicated here.

## Not yet started

funcscan on the single-sample groups — MAP was the explicit ask; funcscan would need its
own gene_export-reuse + samplesheet + config-copy, same pattern, if/when wanted.
