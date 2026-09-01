#!/bin/bash
#SBATCH --account=cbmr
#SBATCH --partition=standardqueue
#SBATCH --job-name=export_single
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --array=0-276%40
#SBATCH --output=/maps/projects/hansen_ol-AUDIT/scratch/NILU/metagenomes/assembly_annotation_wf/single_assembly_production/logs/export_%a.log

source /opt/software/mamba/23.3.1/etc/profile.d/conda.sh
conda activate anvio-9

# Gene export is lightweight (I/O + a quick anvi'o call, no HMM search) unlike kofams'
# 32-CPU tasks -- %40 throttle is safely higher, still considerate of shared cluster use.
GROUP_LIST=/maps/projects/hansen_ol-AUDIT/scratch/NILU/metagenomes/assembly_annotation_wf/single_assembly_production/group_list.txt
P=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$GROUP_LIST")

# Reuse the same contigs-db copies already made for the anvi'o-KEGG run (never the shared
# upstream spa_single_all dbs) -- read-only access here, doesn't conflict with that array's
# concurrent KOfam annotation writes to the same files.
DB=/maps/projects/hansen_ol-AUDIT/scratch/NILU/metagenomes/spa_single_all_anvio_kegg/contigs_dbs/$P.db
OUT=/maps/projects/hansen_ol-AUDIT/scratch/NILU/metagenomes/assembly_annotation_wf/single_assembly_production/gene_export

anvi-get-sequences-for-gene-calls -c "$DB" --export-gff3 -o "$OUT/${P}.gff3"
anvi-get-sequences-for-gene-calls -c "$DB" --get-aa-sequences \
    --defline-format "{contigs_db_project_name}___{gene_caller_id}" -o "$OUT/${P}.faa"

n_gff=$(grep -c $'\tCDS\t' "$OUT/${P}.gff3")
n_faa=$(grep -c "^>" "$OUT/${P}.faa")
echo "$P: gff3_CDS=$n_gff faa=$n_faa"
if [ "$n_gff" != "$n_faa" ]; then
    echo "$P: MISMATCH -- gff3 CDS count != faa sequence count" >&2
    exit 1
fi
