#!/bin/bash
#SBATCH --account=cbmr
#SBATCH --partition=standardqueue
#SBATCH --job-name=kofams_coassembly
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --array=0-18
#SBATCH --output=/maps/projects/hansen_ol-AUDIT/scratch/NILU/metagenomes/assembly_annotation_wf/coassembly_production/logs/kofams_%a.log

source /opt/software/mamba/23.3.1/etc/profile.d/conda.sh
conda activate anvio-9

PARTICIPANTS=(mh_p110 mh_p470 mh_p341 mh_p601 mh_p894 mh_p728 mh_p408 mh_p550 mh_p722 mh_p644 mh_p97 mh_p398 mh_p125 mh_p813 mh_p292 mh_p860 mh_p523 mh_p584 mh_p789)
P=${PARTICIPANTS[$SLURM_ARRAY_TASK_ID]}
BASE=/maps/projects/hansen_ol-AUDIT/scratch/NILU/metagenomes/spa_coassembly_all/results_coassembly_all/concoct

anvi-run-kegg-kofams \
  -c "$BASE/$P/$P.db" \
  --kegg-data-dir /maps/projects/hansen_ol-AUDIT/scratch/NILU/metagenomes/anvio_kegg_data \
  -T "$SLURM_CPUS_PER_TASK"

echo "=== $P kofams done, running estimate-metabolism ==="

anvi-estimate-metabolism \
  -c "$BASE/$P/$P.db" \
  --kegg-data-dir /maps/projects/hansen_ol-AUDIT/scratch/NILU/metagenomes/anvio_kegg_data \
  --metagenome-mode \
  -O /maps/projects/hansen_ol-AUDIT/scratch/NILU/metagenomes/assembly_annotation_wf/coassembly_production/anvio_kegg/${P}_metabolism

echo "=== $P fully done ==="
