#!/bin/bash
#SBATCH --job-name=derep_cross_participant_orfs
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --partition=standardqueue
#SBATCH --account=cbmr
#SBATCH --output=/maps/projects/hansen_ol-AUDIT/scratch/NILU/metagenomes/assembly_annotation_wf/dereplication_pilot/logs/%x_%j.log

# ORF-tier counterpart to run_cross_participant.sh -- pools every participant's
# own Tier-2 (gene/protein) cluster REPRESENTATIVES and clusters those across
# all 25 participants, to test within- vs. between-participant gene/AMR-gene
# sharing (not just contigs). ~8.8M representative proteins / 2.84GB pooled
# input -- smaller total bytes than the contig pass despite more sequences,
# since proteins are much shorter than full contigs.

set -euo pipefail

PILOT_DIR=/maps/projects/hansen_ol-AUDIT/scratch/NILU/metagenomes/assembly_annotation_wf/dereplication_pilot
OUT=$PILOT_DIR/cross_participant_orfs
mkdir -p "$OUT/tmp"

source /opt/software/mamba/23.3.1/etc/profile.d/conda.sh

echo "=== pooling all participants' representative ORFs ==="
: > "$OUT/pooled_reps.faa"
for f in "$PILOT_DIR"/p*/orfs/clu_rep_seq.fasta; do
  cat "$f" >> "$OUT/pooled_reps.faa"
done
grep -c "^>" "$OUT/pooled_reps.faa"

conda activate genomad
echo "=== mmseqs easy-cluster across all participants (same settings as the within-participant pass) ==="
mmseqs easy-cluster "$OUT/pooled_reps.faa" "$OUT/clu" "$OUT/tmp" \
  --min-seq-id 0.95 -c 0.9 --cov-mode 1 --cluster-mode 2 --threads "$SLURM_CPUS_PER_TASK"
conda deactivate

echo "=== DONE ==="
wc -l "$OUT/clu_cluster.tsv"
