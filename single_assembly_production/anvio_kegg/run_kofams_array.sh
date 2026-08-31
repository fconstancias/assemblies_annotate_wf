#!/bin/bash
#SBATCH --account=cbmr
#SBATCH --partition=standardqueue
#SBATCH --job-name=kofams_single
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --array=0-276%20
#SBATCH --output=/maps/projects/hansen_ol-AUDIT/scratch/NILU/metagenomes/spa_single_all_anvio_kegg/logs/kofams_%a.log

source /opt/software/mamba/23.3.1/etc/profile.d/conda.sh
conda activate anvio-9

# Same recipe as the co-assembly run (assembly_annotation_wf/coassembly_production/anvio_kegg/):
# operate on our own copy under contigs_dbs/, never the shared upstream spa_single_all
# contigs-dbs, since anvi-run-kegg-kofams mutates the db it's given in place.
# %20 throttle: cluster has 1664 total CPUs across 13 nodes -- 277 concurrent x 32 CPUs
# would need the whole cluster to itself; 20 concurrent is a considerate ceiling that still
# makes real progress (co-assembly's 19-task array ran without issue at similar scale).
GROUP_LIST=/maps/projects/hansen_ol-AUDIT/scratch/NILU/metagenomes/spa_single_all_anvio_kegg/group_list.txt
P=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$GROUP_LIST")
BASE=/maps/projects/hansen_ol-AUDIT/scratch/NILU/metagenomes/spa_single_all_anvio_kegg/contigs_dbs

anvi-run-kegg-kofams \
  -c "$BASE/$P.db" \
  --kegg-data-dir /maps/projects/hansen_ol-AUDIT/scratch/NILU/metagenomes/anvio_kegg_data \
  -T "$SLURM_CPUS_PER_TASK"

echo "=== $P kofams done, running estimate-metabolism ==="

anvi-estimate-metabolism \
  -c "$BASE/$P.db" \
  --kegg-data-dir /maps/projects/hansen_ol-AUDIT/scratch/NILU/metagenomes/anvio_kegg_data \
  --metagenome-mode \
  -O /maps/projects/hansen_ol-AUDIT/scratch/NILU/metagenomes/spa_single_all_anvio_kegg/${P}_metabolism

echo "=== $P fully done ==="
