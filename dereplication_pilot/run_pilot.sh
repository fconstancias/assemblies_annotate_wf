#!/bin/bash
#SBATCH --job-name=derep_pilot
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --partition=standardqueue
#SBATCH --account=cbmr
#SBATCH --output=/maps/projects/hansen_ol-AUDIT/scratch/NILU/metagenomes/assembly_annotation_wf/dereplication_pilot/logs/%x_%j.log

# Dereplication pilot -- see ../dereplication_brainstorm.md for the design rationale
# and /home/ljc444/.claude/plans/happy-moseying-clarke.md for the pilot plan itself.
# Resources bumped (8->16 CPU, 16->64GB, 4h->12h) for the full 25-participant run --
# the original 2-participant pilot (4/8 single-samples) was smaller than the real
# range (up to participant 813's 35 single-sample assemblies), so sizing generously
# for the largest participants rather than re-tuning per participant.
#
# Usage: sbatch run_pilot.sh <participant_tag> <coassembly_group_or_empty> <single_group1> [single_group2 ...]
# e.g.:  sbatch run_pilot.sh p110 mh_p110 spaS30 spaS37 spaS222 spaS223
# For the 6 participants with no co-assembly (142,151,373,431,654,729 -- excluded from
# co-assembly per CO_ASSEMBLY.md, "more samples expected"), pass "" for coassembly_group
# to pool single-sample assemblies only.

set -euo pipefail

PARTICIPANT="$1"
COASSEMBLY="$2"
shift 2
SINGLE_GROUPS=("$@")

BASE=/maps/projects/hansen_ol-AUDIT/scratch/NILU/metagenomes/assembly_annotation_wf/dereplication_pilot/$PARTICIPANT
CA_CONTIGS=/maps/projects/hansen_ol-AUDIT/scratch/NILU/metagenomes/spa_coassembly_all/results_coassembly_all/assembly/megahit
SA_CONTIGS=/maps/projects/hansen_ol-AUDIT/scratch/NILU/metagenomes/spa_single_all/results_spa_single_all/assembly/spades
CA_FAA=/maps/projects/hansen_ol-AUDIT/scratch/NILU/metagenomes/assembly_annotation_wf/coassembly_production/gene_export
SA_FAA=/maps/projects/hansen_ol-AUDIT/scratch/NILU/metagenomes/assembly_annotation_wf/single_assembly_production/gene_export

source /opt/software/mamba/23.3.1/etc/profile.d/conda.sh

mkdir -p "$BASE/contigs/tmp" "$BASE/orfs/tmp"

echo "=== [$PARTICIPANT] Tier 1: pooling contigs ($COASSEMBLY + ${SINGLE_GROUPS[*]}) ==="
: > "$BASE/contigs/pooled.fa"
if [ -n "$COASSEMBLY" ]; then
  cat "$CA_CONTIGS/$COASSEMBLY/final.contigs.reformatted.fa" >> "$BASE/contigs/pooled.fa"
fi
for g in "${SINGLE_GROUPS[@]}"; do
  cat "$SA_CONTIGS/$g/final.contigs.reformatted.fa" >> "$BASE/contigs/pooled.fa"
done
grep -c "^>" "$BASE/contigs/pooled.fa"

conda activate seqtk
seqtk seq -L 1000 "$BASE/contigs/pooled.fa" > "$BASE/contigs/pooled.min1kb.fa"
conda deactivate
echo "contigs >= 1kb:"
grep -c "^>" "$BASE/contigs/pooled.min1kb.fa"

conda activate genomad
echo "=== [$PARTICIPANT] Tier 1: mmseqs easy-cluster (nucleotide, containment-aware) ==="
mmseqs easy-cluster "$BASE/contigs/pooled.min1kb.fa" "$BASE/contigs/clu" "$BASE/contigs/tmp" \
  --min-seq-id 0.95 -c 0.8 --cov-mode 2 --cluster-mode 2 --threads "$SLURM_CPUS_PER_TASK"

echo "=== [$PARTICIPANT] Tier 2: pooling ORFs ==="
: > "$BASE/orfs/pooled.faa"
if [ -n "$COASSEMBLY" ]; then
  cat "$CA_FAA/$COASSEMBLY.faa" >> "$BASE/orfs/pooled.faa"
fi
for g in "${SINGLE_GROUPS[@]}"; do
  cat "$SA_FAA/$g.faa" >> "$BASE/orfs/pooled.faa"
done
grep -c "^>" "$BASE/orfs/pooled.faa"

echo "=== [$PARTICIPANT] Tier 2: mmseqs easy-cluster (protein) ==="
mmseqs easy-cluster "$BASE/orfs/pooled.faa" "$BASE/orfs/clu" "$BASE/orfs/tmp" \
  --min-seq-id 0.95 -c 0.9 --cov-mode 1 --cluster-mode 2 --threads "$SLURM_CPUS_PER_TASK"
conda deactivate

echo "=== [$PARTICIPANT] clustering done, analyzing ==="
conda activate r-binner-compare
# NOTE: "GROUPS" is a bash special variable (the process's Unix group-ID list) --
# assigning a custom value to it silently collides with that instead of holding our
# own string, so $GROUPS below would NOT expand to what you'd expect. Confirmed the
# hard way (100% of cluster members came back "UNKNOWN" in a first real run this
# produced -- $GROUPS had actually expanded to a real numeric GID, not our group
# list). Use a differently-named variable for anything script-local.
PARTICIPANT_GROUPS="$COASSEMBLY ${SINGLE_GROUPS[*]}"
Rscript /maps/projects/hansen_ol-AUDIT/scratch/NILU/metagenomes/assembly_annotation_wf/dereplication_pilot/R/analyze_clusters.R \
  --cluster-tsv "$BASE/contigs/clu_cluster.tsv" --fasta "$BASE/contigs/pooled.min1kb.fa" \
  --groups $PARTICIPANT_GROUPS --tier contigs --participant "$PARTICIPANT" \
  --out "$BASE/contigs/cluster_membership.tsv" --summary "$BASE/contigs/summary.txt"

Rscript /maps/projects/hansen_ol-AUDIT/scratch/NILU/metagenomes/assembly_annotation_wf/dereplication_pilot/R/analyze_clusters.R \
  --cluster-tsv "$BASE/orfs/clu_cluster.tsv" --fasta "$BASE/orfs/pooled.faa" \
  --groups $PARTICIPANT_GROUPS --tier orfs --participant "$PARTICIPANT" \
  --out "$BASE/orfs/cluster_membership.tsv" --summary "$BASE/orfs/summary.txt"

echo "=== [$PARTICIPANT] DONE ==="
cat "$BASE/contigs/summary.txt" "$BASE/orfs/summary.txt"
