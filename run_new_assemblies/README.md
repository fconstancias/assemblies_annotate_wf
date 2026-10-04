# Annotate new assemblies: anvi'o contigs DB → funcscan + MAP

AMR, virulence/toxins, AMPs, CAZymes, BGCs, plasmids/phages/IS/integrons for any set of anvi'o contigs DBs (single
assemblies, co-assemblies, a MAG catalogue). funcscan and MAP run on **anvi'o's own gene calls**, so every hit carries
the anvi'o gene ID (`<project>___<gene_callers_id>`) and joins back to anvi'o (coverage, bins, collections).
Tested end to end on matph (2026-10). Background and history: `../REPRODUCE.md`, `../CLAUDE.md`.

## Quick start (esrum)

```bash
D=/maps/projects/hansen_ol-AUDIT/scratch/NILU/metagenomes/assembly_annotation_wf/run_new_assemblies
RUN=/path/to/my_run                     # one folder per run; everything is written inside
mkdir -p $RUN && cd $RUN

# 0. once per user: patch the pipelines (idempotent)
$D/scripts/01_apply_patches.sh

# 1. list the contigs DBs (tab-separated, header line)
printf "sample\tcontigs_db\n" > samples.tsv
printf "S1\t/path/to/S1.db\n" >> samples.tsv          # one line per DB; sample = a short unique name

# 2. export contigs + genes from anvi'o (SLURM array, ~5 min per DB)
sbatch --array=1-$(tail -n +2 samples.tsv | wc -l) $D/scripts/02_export_from_anvio.sbatch $RUN

# 3. samplesheets
$D/scripts/03_make_samplesheets.sh $RUN

# 4. launch both pipelines (drivers run in the background with nohup; tasks go to SLURM)
$D/scripts/04_launch_funcscan.sh $RUN        # log: $RUN/funcscan/nohup_funcscan.log
$D/scripts/05_launch_map.sh $RUN             # log: $RUN/map/nohup_map.log

# 5. when both are finished: checks + SignalP fill
sbatch $D/scripts/06_post_run.sbatch $RUN    # status 0 = everything accounted for

# 6. optional, per sample: all annotations into a COPY of the contigs DB
sbatch $D/../scripts/import_annotations_into_anvio.sbatch /path/to/S1.db $RUN/gene_export/S1.gff3 S1 \
       $RUN/map/results/S1 $RUN/funcscan/results $RUN/annotated/S1
```
Re-running the same launch command **resumes** (finished tasks are reused). Add `-preview` as second argument to a
launch script to only check samplesheet, config and patches (nothing runs): `$D/scripts/05_launch_map.sh $RUN -preview`.

## Install (once)

| What | How |
|---|---|
| Nextflow env | `mamba create -n env_nf -c conda-forge -c bioconda nextflow` then `conda activate env_nf && nextflow self-update` (needs ≥ 25.10.4 for funcscan 4.0.0). One normal conda env, not nested. |
| anvi'o | env `anvio-9` (https://anvio.org/install/), only for steps 2 and 6 |
| Singularity | `module load singularity/3.8.7` (done by `scripts/00_env.sh`); containers are pulled by Nextflow into the cache set in the configs |
| Pipelines + patches | `scripts/01_apply_patches.sh` (pulls funcscan 4.0.0 and MAP v5.0.0 if missing, applies `patches/`) |
| MAP databases | already on esrum, paths in `config/map_db_paths.config`; new cluster: `../REPRODUCE.md` §4 |
| funcscan databases | downloaded by funcscan itself into the run's `work/` |

## Options

**Input**: any anvi'o contigs DB with gene calls (`anvi-gen-contigs-database`; pyrodigal-gv / prodigal). Several DBs
= several lines in `samples.tsv`; they run in parallel in one funcscan + one MAP run.

**funcscan** (`scripts/04_launch_funcscan.sh`): runs ARG (ABRicate, AMRFinderPlus, DeepARG, fARGene, RGI →
hAMRonization, argNorm), CAZymes (dbCAN), AMPs (ampir → AMPcombi). Edit the launch line to change screens
(nf-core/funcscan parameters). BGC screening is **not** used: funcscan's BGC tools need real GenBank input; BGCs come
from MAP (SanntiS).

**MAP** (`scripts/05_launch_map.sh`): geNomad + CheckV, ISEScan, IntegronFinder, compositional outliers, PathoFact2
(virulence/toxins, VFDB), AMRFinderPlus/RGI/DeepARG, SanntiS + InterProScan. In `map/nextflow.config` (`params {}`):
`skip_sanntis = true` drops SanntiS **and** InterProScan (saves ~15–25 h per sample, but no BGCs, Pfam or SignalP);
`skip_gecco` / `skip_antismash` stay `true`. ICEfinder2 is disabled by patch 04.

**Resources**: `../templates/*.nextflow.config` (copied into the run at first launch; edit the run's copy for that run,
the template for all future runs). Typical durations per sample of 300–500k genes: InterProScan 9–17 h, SanntiS 4–7 h,
dbCAN ~5 h, the rest < 1 h. Time and memory already set for those.

**Other cluster / user**: edit in `../templates/*.nextflow.config` `clusterOptions = '--account=…'`, `queue`,
`singularity.cacheDir`; in `scripts/00_env.sh` the conda path, env name and singularity module; in
`config/map_db_paths.config` the DB paths.

## Outputs (in `$RUN`)

| Path | Content |
|---|---|
| `gene_export/<S>.{contigs.fa,gff3,faa,fna}` | input given to both pipelines (anvi'o IDs) |
| `funcscan/results/` | `reports/hamronization_summarize/` (all ARG tools), `reports/ampcombi2/`, `cazyme/dbcan/…/<S>_overview.tsv`, per-tool folders |
| `map/results/<S>/` | `<S>_combined_report.signalp.tsv` (per protein: VFDB, toxin/virulence prob, AMR, MGE, BGC, SignalP), `gff/<S>_mobilome.gff.gz`, `prediction/` (geNomad, CheckV, ISEScan, IntegronFinder, AMR tools, PathoFact2, BGC, InterProScan), `preprocessing/<S>_contigID.map` (MAP renames contigs `contig_N`) |
| `annotated/<S>/` (step 6) | `<S>_annotated.db` + `functions/` (functions-txt, `evidence_amr.tsv`, `evidence_vf_toxin.tsv`, `contig_genomad.tsv`) |

## Know before interpreting (checked by step 5)

- MAP's mobilome keeps geNomad calls with score > 0.8, complete IS elements and complete integrons only.
- MAP's combined report lists only PathoFact2 + AMR proteins; its BGC / MGE columns are filled for those only.
- hAMRonization drops AMRFinderPlus hits with an internal stop codon (listed as WARN).
- MAP's `signalP` column is always empty in v5.0.0 → filled in `<S>_combined_report.signalp.tsv` (step 5).
- Compositional outliers are weak evidence (composition only), not mobility.

## Layout

`scripts/` 00–06 above · `patches/funcscan/` 01–03, `patches/map/` 04–06 (what each fixes: `../README.md`
"Patches") · `config/map_db_paths.config` (MAP DB paths, comments explain each DB workaround) · configs:
`../templates/` · checks used by 06: `../scripts/check_map_outputs.sh`, `../scripts/check_annotation_completeness.py`,
`../scripts/add_signalp_to_report.py`.

**When a run hits a new problem**: fix it in `../templates/` (config) or as a new patch in `patches/` (+ a line in
`../README.md` "Patches"), and in the checks if it was silent, so the next run starts from it.
