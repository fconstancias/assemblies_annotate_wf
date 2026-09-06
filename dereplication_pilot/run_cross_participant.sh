#!/bin/bash
#SBATCH --job-name=derep_cross_participant
#SBATCH --cpus-per-task=32
#SBATCH --mem=256G
#SBATCH --time=24:00:00
#SBATCH --partition=standardqueue
#SBATCH --account=cbmr
#SBATCH --output=/maps/projects/hansen_ol-AUDIT/scratch/NILU/metagenomes/assembly_annotation_wf/dereplication_pilot/logs/%x_%j.log

# Second-level clustering: pool every participant's own Tier-1 (contig) cluster
# REPRESENTATIVES (already deduplicated within-participant by the first-pass
# run_pilot.sh) and cluster those across all 25 participants. This tests
# within- vs. between-participant redundancy without the much larger cost of
# re-clustering all 296 assemblies' raw contigs together from scratch --
# same two-level design real gene-catalog studies use (dedupe locally, merge
# at the next level). ~5.48M representative contigs / 7.9GB pooled input.

set -euo pipefail

PILOT_DIR=/maps/projects/hansen_ol-AUDIT/scratch/NILU/metagenomes/assembly_annotation_wf/dereplication_pilot
OUT=$PILOT_DIR/cross_participant
mkdir -p "$OUT/tmp"

source /opt/software/mamba/23.3.1/etc/profile.d/conda.sh

echo "=== pooling all participants' representative contigs ==="
: > "$OUT/pooled_reps.fa"
for f in "$PILOT_DIR"/p*/contigs/clu_rep_seq.fasta; do
  cat "$f" >> "$OUT/pooled_reps.fa"
done
grep -c "^>" "$OUT/pooled_reps.fa"

conda activate genomad
echo "=== mmseqs easy-cluster across all participants (same settings as the within-participant pass) ==="
mmseqs easy-cluster "$OUT/pooled_reps.fa" "$OUT/clu" "$OUT/tmp" \
  --min-seq-id 0.95 -c 0.8 --cov-mode 2 --cluster-mode 2 --threads "$SLURM_CPUS_PER_TASK"
conda deactivate

echo "=== DONE ==="
wc -l "$OUT/clu_cluster.tsv"
