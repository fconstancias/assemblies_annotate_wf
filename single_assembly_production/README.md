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

## Status: ✅ complete, 2026-09-03

All 277 groups finished MAP's full pipeline (real SANNTIS/BGC output confirmed for every
group) and are uploaded to S3 — the sync watcher (`s3_sync_watcher.sh`) exited cleanly
having reached "everything complete" (277/277 real `UPLOADED:` log entries, not a crash).

## anvi'o-KEGG BadConstraints incident, 2026-09-03

Job 1851734 (`../../spa_single_all_anvio_kegg/run_kofams_array.sh`, `--array=0-276%20`)
looked complete at "196/196" earlier, but that was actually a partial/stuck state, not
full completion — indices 196-276 (81 groups) were stuck `PENDING`/`BadConstraints`
indefinitely. Confirmed real: `sacct` showed 0-195 genuinely `COMPLETED`, direct file
checks confirmed all 81 groups in the 196-276 range genuinely had no
`*_metabolism_modules.txt` output. Root cause wasn't real resource unavailability —
`sinfo` showed 7 fully idle nodes and 1538/1664 idle CPUs cluster-wide at the time, no
account/QOS submit or job-count limits configured (`sacctmgr show assoc`/`show qos` both
empty) — `scontrol show job` on one of the stuck tasks showed internally inconsistent
resource fields (`NumCPUs=6-6`/`CPUs/Task=8` alongside `MinCPUsNode=32`) pointing to a
stale/corrupted SLURM job record for just that pending range, not a real constraint.
Fixed by `scancel`-ing just the stuck pending range (`1851734_[196-276]`, leaving the 196
genuinely-completed tasks untouched) and resubmitting the same script scoped to just
those indices (`sbatch --array=196-276%20 run_kofams_array.sh`, job 1874380) — the array
task ID → group_list.txt line-number mapping is unchanged, so this fills exactly the gap.
Confirmed fixed: the resubmission scheduled normally (`JobArrayTaskLimit`, not
`BadConstraints`) with real tasks running on real nodes immediately.

## funcscan — started 2026-09-03

Same pattern as `map_run/`: `funcscan_run/nextflow.config` is a direct copy of
`../coassembly_production/funcscan_run/nextflow.config` (reuses every fix already proven
there — AMPIR/RGI_MAIN/RUNDBCAN memory+time overrides, AMPCOMBI2_CLUSTER
`errorStrategy = 'ignore'`, the blanket per-attempt memory scaling), with `queueSize`
trimmed from 200 to 100 since this run shares the cluster with the concurrent
single-sample anvi'o-KEGG array (up to 640 CPUs at `%20` throttle) on top of its own
277-group scope. `samplesheet_single.csv` reuses the already-exported, already-verified
`gene_export/*.{gff3,faa}` plus the same `spa_single_all/.../spades/{group}/final.contigs.reformatted.fa`
assembly-path convention MAP's samplesheet uses — verified zero missing files across all
277 groups before launch. Launched with the same flags as both prior funcscan runs
(`--run_arg_screening --run_cazyme_screening --run_amp_screening --amp_skip_amplify
--amp_skip_macrel`), named run `single_assembly_production_funcscan` for future
`-resume`. Confirmed launched cleanly — real per-sample tasks (ABRICATE, AMPIR, FARGENE,
RUNDBCAN_DATABASE, ...) submitting and running on real nodes within a minute of launch.

`funcscan_s3_sync_watcher.sh` covers S3 upload for this run — same idempotent
`rclone copy --checksum` pattern as every other watcher here, but funcscan's
`results/` layout is tool-first (`results/<category>/<tool>/<sample>/...`), not
sample-first like MAP's, across 8+ tool directories — a granular per-sample-per-tool
upload loop (like the coassembly watcher's) isn't worth the fragility at this scale, so
this one does periodic whole-tree syncs instead, gated by counting how many samples have
reached dbCAN's substrate-prediction step (the last stage of RUNDBCAN, funcscan's
slowest tool chain) as the "worth syncing now" proxy, with a final sync once all 277
groups reach it.
