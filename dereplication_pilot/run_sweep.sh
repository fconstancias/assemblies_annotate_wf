#!/bin/bash
#SBATCH --job-name=derep_sweep
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=2:00:00
#SBATCH --partition=standardqueue
#SBATCH --account=cbmr
#SBATCH --output=/maps/projects/hansen_ol-AUDIT/scratch/NILU/metagenomes/assembly_annotation_wf/dereplication_pilot/logs/sweep_%x_%j.log

# Parameter-sensitivity sweep for the dereplication pilot -- reuses the already-pooled
# fasta files from run_pilot.sh's first pass (no need to re-pool), just re-clusters
# with different mmseqs settings to see how much the redundancy numbers depend on the
# specific --min-seq-id/-c/--cov-mode/--cluster-mode choices made in the first pass.
#
# Usage: sbatch run_sweep.sh <participant> <tier: contigs|orfs> <variant_tag> \
#          <input_fasta> "<mmseqs extra args, one quoted string>" "<participant groups, one quoted string>"
#
# e.g.: sbatch run_sweep.sh p110 contigs covmode0 \
#         .../p110/contigs/pooled.min1kb.fa "--min-seq-id 0.95 -c 0.8 --cov-mode 0 --cluster-mode 2" \
#         "mh_p110 spaS30 spaS37 spaS222 spaS223"
#
# NOTE: do NOT name the participant-groups variable "GROUPS" -- that collides with
# bash's own special GROUPS variable (process's Unix group IDs) and silently breaks
# cluster->source-group attribution. Bit us for real in run_pilot.sh; using
# PARTICIPANT_GROUPS_STR here instead.

set -euo pipefail

PARTICIPANT="$1"
TIER="$2"
VARIANT_TAG="$3"
INPUT_FASTA="$4"
MMSEQS_ARGS_STR="$5"
PARTICIPANT_GROUPS_STR="$6"

read -ra MMSEQS_ARGS <<< "$MMSEQS_ARGS_STR"

BASE=/maps/projects/hansen_ol-AUDIT/scratch/NILU/metagenomes/assembly_annotation_wf/dereplication_pilot/$PARTICIPANT/$TIER/sweep_$VARIANT_TAG
mkdir -p "$BASE/tmp"

source /opt/software/mamba/23.3.1/etc/profile.d/conda.sh
conda activate genomad
echo "=== [$PARTICIPANT/$TIER/$VARIANT_TAG] mmseqs easy-cluster ${MMSEQS_ARGS[*]} ==="
mmseqs easy-cluster "$INPUT_FASTA" "$BASE/clu" "$BASE/tmp" "${MMSEQS_ARGS[@]}" --threads "$SLURM_CPUS_PER_TASK"
conda deactivate

echo "=== [$PARTICIPANT/$TIER/$VARIANT_TAG] analyzing ==="
conda activate r-binner-compare
Rscript /maps/projects/hansen_ol-AUDIT/scratch/NILU/metagenomes/assembly_annotation_wf/dereplication_pilot/R/analyze_clusters.R \
  --cluster-tsv "$BASE/clu_cluster.tsv" --fasta "$INPUT_FASTA" \
  --groups $PARTICIPANT_GROUPS_STR --tier "$TIER" --participant "${PARTICIPANT}_${VARIANT_TAG}" \
  --out "$BASE/cluster_membership.tsv" --summary "$BASE/summary.txt"

echo "=== [$PARTICIPANT/$TIER/$VARIANT_TAG] DONE ==="
cat "$BASE/summary.txt"
