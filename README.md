# assembly_annotation_wf

Functional annotation and longitudinal quantification for the NCCR fecal-metagenome
cohort's assemblies (single-sample and co-assembly), built on top of
`fconstancias/shotgun_metagenomics_binning`'s contigs/bins/gene-calls. Combines three
independent pipelines against the same upstream anvi'o contigs-dbs: nf-core/funcscan
(AMR/AMP/CAZyme), EBI's mobilome-annotation-pipeline (MAP: MGE/plasmid/BGC), and anvi'o's
own KEGG KOfam + metabolism-estimation.

**Start here, not this file, if you're picking this project back up**: `CLAUDE.md` is the
maintained operational memory (gotchas, fixes, current state) — read it first, every
session. This README is the stable orientation layer; CLAUDE.md is where the fast-moving
detail lives.

## Documentation map

| File | Purpose |
|---|---|
| `CLAUDE.md` | Operational memory: what's been tried, what broke, current production state. Read first. |
| `COMMANDS.md` | Literal commands actually run this project, in order — copy-paste reference. |
| `REPRODUCE.md` | Generic/templated version of the same, for adapting to new samples or a new cluster. |
| `REPRODUCE.md` §8 | Config gotchas worth applying proactively rather than rediscovering. |
| This file | Orientation: what this repo is, where things live, patch inventory. |

## Layout

- `funcscan_run/`, `map_run/`, `funcscan_gff_test/` — early single/curated-sample test runs
  (2 samples: `megaS121`, `co728`) that established the working recipe before production.
- `coassembly_production/` — the real production run against all 19 co-assemblies
  (`spa_coassembly_all/`), started 2026-08-29 under a hard access-deadline (2026-08-31).
  Contains `funcscan_run/`, `map_run/`, `anvio_kegg/`, `gene_export/`, and
  `s3_sync_watcher.sh` (incremental S3 upload as groups finish, see below).
- `gene_export/`, `bgc_comparison/`, `metabolism_comparison/` — analysis/comparison work,
  not pipeline runs themselves.
- `scripts/` — standalone helpers, notably `gff3_to_gbk.py` (real GenBank conversion for
  AMPCOMBI2, see CLAUDE.md's Downstream section for why this exists).

Single-sample-assembly KEGG annotation (`spa_single_all/`'s 277 groups) runs in a **separate
directory outside this repo**: `../spa_single_all_anvio_kegg/` (own `run_kofams_array.sh` +
`s3_sync_watcher.sh`, mirroring `coassembly_production/anvio_kegg/`'s recipe exactly — copies
of both scripts are kept here under `single_assembly_production/anvio_kegg/` for reference/
reproducibility, but the live, running data is not in this repo). Started 2026-08-31, in
progress — see CLAUDE.md's Production run section for current status.

## Patches

Two upstream pipelines needed local patches to work against Prodigal-called genes from our
own anvi'o contigs-dbs (neither pipeline is designed around "genes already called upstream"
as the default case):

| # | File | Pipeline | What it fixes |
|---|---|---|---|
| 1 | `funcscan_gff_test/patches/01_schema_drop_protein_requires_gbk.patch` | funcscan | Input schema wrongly required a `.gbk` column whenever `protein` was supplied; we supply protein without gbk. |
| 2 | `funcscan_gff_test/patches/02_ampcombi2_parsetables_optional_gbk.patch` | funcscan | Makes AMPCOMBI2's `--gbk` optional, generating a Biopython-built placeholder `.gbff` when real GenBank isn't supplied. |
| 3 | `funcscan_gff_test/patches/03_ampcombi_download_dramp_nan_fix.patch` | funcscan | `ampcombi_download.py` crashed on non-string (NaN) `Sequence` rows in the DRAMP reference table. |
| 4 | `coassembly_production/map_run/04_disable_broken_icefinder2_lite.patch` | MAP | ICEFINDER2_LITE's internal join crashes whenever one of its parallel branches produces no match — confirmed recurring across multiple participants, not an isolated case. Disabled the subworkflow entirely (substituting an empty channel the existing `remainder:true` downstream already handles) rather than excluding affected samples one at a time. |

Patches 1-3 are generated from hand-edited files in the staged Nextflow checkout
(`~/.nextflow/assets/nf-core/funcscan/`) via `git diff` — see `COMMANDS.md` §3 for the exact
commands and `CLAUDE.md`'s Downstream section for the two hardest-won specifics (the
placeholder `.gbff` must be Biopython-generated, and the placeholder logic must only fire
when `gbk_input` is genuinely absent). Patch 4 came from the same workflow, applied directly
to `workflows/mobilomeannotation.nf` in the staged MAP checkout.

MAP additionally needed two `nextflow.config`-level `errorStrategy = 'ignore'` overrides
(not code patches) for upstream bugs unrelated to our own input format — `DB_DOWNLOAD_VFDB`'s
missing curl dependency, and geNomad's mismatched tarball folder name. See
`coassembly_production/map_run/nextflow.config`.

Beyond these, a long tail of resource-limit fixes (OOM ceilings, one SLURM-timeout
misdiagnosed as OOM) were needed to get funcscan through at real 19-participant production
scale — not code patches, just `nextflow.config` `withName` resource overrides. Full history
in `git log --oneline`; each commit message states the actual failure mode and fix, not just
"increase memory."

## S3 backup

Every project directory under `../` (this repo, `spa_single_all/`, `spa_coassembly_all/`,
`spa_assembly_comparison/`) is backed up to a UNIL S3 bucket
(`recn-fac-fbm-dbc-slehtine-stool-sampling`), one prefix per project, via `rclone copy
--checksum` (real hash verification, not just size/modtime) followed by an independent
`rclone check` before trusting any upload complete. Known-broken/junk data (a since-fixed
dRep bug's stale output, verbose SLURM logs) is explicitly excluded — see each project's own
upload script for its exact exclude list; this repo excludes `coassembly_production/`'s
still-in-progress pieces (MAP's BGC/SANNTIS step, funcscan's CAZyme substrate step for
whichever groups haven't finished yet) via `coassembly_production/s3_sync_watcher.sh`, which
polls every 15 minutes and uploads each group's output the moment it's genuinely complete
(not just "a file exists" — cross-checked against whether a job is still actively writing to
that group). Uploads run as their own SLURM batch jobs, not head-node background processes,
so they survive session/login disconnects; they do **not** survive a full account
suspension (job-kill), which is the access-cutoff scenario assumed for 2026-08-31 — anything
still mid-computation at that moment is lost and needs rerunning once access resumes, but
everything already finished and synced is safe regardless.
