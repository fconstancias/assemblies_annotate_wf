# Exact commands used this session — assemblies → functional annotation

Literal, copy-paste commands actually run, in order, on the real data/paths used this
session (as opposed to `REPRODUCE.md`'s generic/templated version — use that one when
adapting to different samples/paths, use this one to see exactly what was actually
executed). Paths below are real, from the `hansen_ol-AUDIT` project on the esrum cluster.

## 0. Environment activation (needed before every block below)

```bash
source /opt/software/mamba/23.3.1/etc/profile.d/conda.sh
source /usr/share/Modules/init/bash
```
Non-interactive shells don't source `.bashrc`'s conda-init block (gated behind an
interactive-only check) — always source directly like this first.

## 1. Install: conda environments

```bash
# anvi'o -- assumed pre-existing as `anvio-9` in this project, not built fresh here

# Dedicated Nextflow env
mamba create -n env_nf -c bioconda -c conda-forge nextflow -y
conda activate env_nf
nextflow self-update   # 25.10.2 -> 26.04.6, needed for funcscan 4.0.0's version floor

# gff3->gbk converter env (only needed for real per-participant GenBank / AMPCOMBI2_CLUSTER)
mamba create -n gff2gbk -c conda-forge python=3.10 biopython -y
conda activate gff2gbk
python3 -m pip install --no-deps bcbio-gff   # bioconda's own package pins an incompatible
python3 -m pip install six                   #   old biopython -- don't `conda install` it

# Singularity module (already installed system-wide, just needs loading each session)
module load singularity/3.8.7
```

## 2. Install: pipeline checkouts (first `nextflow run` auto-stages these)

```bash
conda activate env_nf
module load singularity/3.8.7
export TMPDIR=/tmp   # always -- see REPRODUCE.md §8

# Triggers staging at ~/.nextflow/assets/nf-core/funcscan/
nextflow run nf-core/funcscan -r 4.0.0 --help

# Triggers staging at ~/.nextflow/assets/EBI-Metagenomics/mobilome-annotation-pipeline/
nextflow run EBI-Metagenomics/mobilome-annotation-pipeline -r v5.0.0 --help
```

## 3. Patch: funcscan (applied directly to the staged checkout)

```bash
cd ~/.nextflow/assets/nf-core/funcscan

# Patch 1 -- schema_input.json: drop "protein": ["gbk"] from dependentRequired
# (edited assets/schema_input.json by hand, then:)
git diff -- assets/schema_input.json > /path/to/funcscan_gff_test/patches/01_schema_drop_protein_requires_gbk.patch

# Patch 2 -- ampcombi2/parsetables/main.nf: optional --gbk with Biopython-generated placeholder
git diff -- modules/nf-core/ampcombi2/parsetables/main.nf > /path/to/funcscan_gff_test/patches/02_ampcombi2_parsetables_optional_gbk.patch

# Patch 3 -- ampcombi_download.py: skip non-string (NaN) Sequence rows
git diff -- bin/ampcombi_download.py > /path/to/funcscan_gff_test/patches/03_ampcombi_download_dramp_nan_fix.patch

# Verify a patch applies cleanly against the pristine checkout before trusting it:
git stash push -- <file>
git apply --check /path/to/patch
git stash pop
```
The actual file edits these patches came from are described in full in `CLAUDE.md`'s
Downstream section — the two hardest-won specifics: (a) the placeholder `.gbff` must be
generated with Biopython itself, not hand-written (GenBank's LOCUS line is a rigid
fixed-column format Biopython's own parser rejects otherwise); (b) that placeholder logic
must only run when `gbk_input` is genuinely absent — a real supplied gbk file crashes if
treated as a directory to `mkdir`/`ls` into.

MAP additionally needed one code patch, applied the same way against the staged checkout's
`workflows/mobilomeannotation.nf`:

```bash
cd ~/.nextflow/assets/EBI-Metagenomics/mobilome-annotation-pipeline
# Patch 4 -- disable ICEFINDER2_LITE entirely: its internal join logic crashes whenever one
# of its parallel branches produces no match (confirmed recurring across multiple
# participants, not isolated) -- substituted an empty channel the existing remainder:true
# downstream already handles, rather than excluding affected samples one at a time.
git diff -- workflows/mobilomeannotation.nf > /path/to/coassembly_production/map_run/04_disable_broken_icefinder2_lite.patch
```

Beyond that one patch, MAP needed only `nextflow.config`-level `errorStrategy = 'ignore'`
overrides for two upstream bugs (`DB_DOWNLOAD_VFDB`'s missing curl, geNomad's mismatched
tarball folder name) — see `map_run/nextflow.config` for the by-hand DB rescue each one
needs.

## 4. Reference database setup

```bash
conda activate anvio-9

# anvi'o KEGG/KOfam -- one-time, ~10GB, no login needed
mkdir -p /maps/projects/hansen_ol-AUDIT/scratch/NILU/metagenomes  # parent must pre-exist
anvi-setup-kegg-data --mode all \
    --kegg-data-dir /maps/projects/hansen_ol-AUDIT/scratch/NILU/metagenomes/anvio_kegg_data \
    -T 8
# NOTE: do NOT pre-mkdir the kegg-data-dir itself -- anvi'o refuses to touch a directory
# that already exists, as a safety check. Let it create that one.
```

## 5. Per-sample gene export (real commands, real sample)

```bash
conda activate anvio-9
DB=/maps/projects/hansen_ol-AUDIT/scratch/NILU/metagenomes/BINNING/01_ANVIO_DBs/megaS121/megaS121.db

anvi-get-sequences-for-gene-calls -c "$DB" --export-gff3 -o megaS121.gff3

anvi-get-sequences-for-gene-calls -c "$DB" --get-aa-sequences \
    --defline-format "{contigs_db_project_name}___{gene_caller_id}" -o megaS121.faa

anvi-get-sequences-for-gene-calls -c "$DB" \
    --defline-format "{contigs_db_project_name}___{gene_caller_id}" -o megaS121.fna

# Verification that actually caught the ID-mismatch bug in the first place:
grep -c $'\tCDS\t' megaS121.gff3     # must equal...
grep -c "^>" megaS121.faa            # ...this, exactly
```

## 6. Run: funcscan (real launch command, production co-assembly scale)

```bash
conda activate env_nf
module load singularity/3.8.7
export TMPDIR=/tmp
cd coassembly_production/funcscan_run

nextflow run nf-core/funcscan -r 4.0.0 \
    --input samplesheet_coassembly.csv --outdir results \
    --run_arg_screening --run_cazyme_screening \
    --run_amp_screening --amp_skip_amplify --amp_skip_macrel \
    -c nextflow.config -profile singularity -resume coassembly_production
```

## 7. Run: MAP (real launch command)

```bash
conda activate env_nf
module load singularity/3.8.7
export TMPDIR=/tmp
cd coassembly_production/map_run

nextflow run EBI-Metagenomics/mobilome-annotation-pipeline -r v5.0.0 \
    --input samplesheet_coassembly.csv --outdir results \
    -c nextflow.config \
    -c /maps/projects/hansen_ol-AUDIT/scratch/NILU/metagenomes/assembly_to_MGE/my_paths.config \
    -profile singularity -resume coassembly_production_v2
```
`skip_sanntis`/`skip_gecco`/`skip_antismash` are set as real Groovy booleans inside
`nextflow.config`'s `params {}` block, not as CLI flags — see REPRODUCE.md §6 for why.

## 8. Run: anvi'o-KEGG (real SLURM array submission)

```bash
sbatch coassembly_production/anvio_kegg/run_kofams_array.sh
```
That script itself runs, per array task (one per participant):
```bash
conda activate anvio-9
anvi-run-kegg-kofams -c "$DB" \
    --kegg-data-dir /maps/projects/hansen_ol-AUDIT/scratch/NILU/metagenomes/anvio_kegg_data \
    -T "$SLURM_CPUS_PER_TASK"
anvi-estimate-metabolism -c "$DB" \
    --kegg-data-dir /maps/projects/hansen_ol-AUDIT/scratch/NILU/metagenomes/anvio_kegg_data \
    --metagenome-mode -O "${P}_metabolism"
```

## 9. Run: anvi'o-KEGG for single-sample assemblies (277 groups, separate dir)

Same recipe as §8, against `spa_single_all/`'s 277 single-sample contigs-dbs instead of the
19 co-assemblies. Runs in `../spa_single_all_anvio_kegg/` (outside this repo — copies of both
scripts kept here under `single_assembly_production/anvio_kegg/` for reference). Own copy of
each contigs-db first (never mutate the shared upstream `spa_single_all` dbs, same convention
as §8), array throttled to 20 concurrent tasks (`%20` — cluster has 1664 total CPUs across 13
nodes, 277 concurrent x 32 CPUs/task would need the whole cluster to itself).

```bash
sbatch spa_single_all_anvio_kegg/run_kofams_array.sh
sbatch spa_single_all_anvio_kegg/s3_sync_watcher.sh
```
Real per-group runtime for co-assembly (§8) was 7.5-15h each (230K-667K genes/participant);
single-sample groups are smaller (132K-321K genes/group sampled) but still multi-hour, not
minutes — full 277-group batch expected to take multiple days even with throttled
parallelism, not something to wait on interactively.

## 10. S3 incremental sync watchers

Both `coassembly_production/s3_sync_watcher.sh` and
`spa_single_all_anvio_kegg/s3_sync_watcher.sh` poll every 15 minutes for newly-finished work
(cross-checked against `squeue` so a group mid-(re)computation is never uploaded prematurely)
and upload it via `rclone copy --checksum`. Both run as their own long-lived SLURM batch jobs
(`--time=7-00:00:00`), not head-node background processes, so they survive session/login
disconnects (confirmed: a head-node Snakemake driver survived 8+ days untouched, and MAP's
own Nextflow driver is `nohup`'d and reparented to init — a plain SSH disconnect does not
kill either). They do **not** survive a full account suspension (job-kill), the assumed
access-cutoff scenario for 2026-08-31 -- see `README.md`'s S3 backup section.

```bash
sbatch coassembly_production/s3_sync_watcher.sh
sbatch spa_single_all_anvio_kegg/s3_sync_watcher.sh
```

## 11. Operational commands used repeatedly (not install/patch/run, but load-bearing)

```bash
# Clear a stale Nextflow session lock after confirming no process actually holds it
ps aux | grep -i nextflow | grep -v grep    # confirm empty first
rm -f .nextflow/cache/*/db/LOCK

# Check real cluster capacity before trusting a queueSize value
sinfo -p standardqueue -o "%D %C"
squeue -u $USER | wc -l
```
