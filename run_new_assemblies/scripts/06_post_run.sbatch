#!/bin/bash
#SBATCH --account=cbmr
#SBATCH --partition=standardqueue
#SBATCH --job-name=post_run_checks
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=4:00:00
#SBATCH --output=post_run_checks_%j.log
# After funcscan + MAP finished: (1) MAP silent-failure check, (2) SignalP fill of each MAP combined report,
# (3) raw -> summary accounting per sample. Read-only on results except the new *_combined_report.signalp.tsv.
#   sbatch scripts/06_post_run.sbatch <run_dir>      (log: post_run_checks_<jobid>.log in the submit dir)
set -uo pipefail
RUN=$(readlink -f ${1:?run dir})
D=${REPO:-/maps/projects/hansen_ol-AUDIT/scratch/NILU/metagenomes/assembly_annotation_wf}   # repo root (sbatch runs a copy of this script)
M=$RUN/map/results; F=$RUN/funcscan/results; rc=0
echo "########## 1. MAP output check"; $D/scripts/check_map_outputs.sh $M $RUN/map/work || rc=1
for S in $(tail -n +2 $RUN/samples.tsv | cut -f1); do
  echo "########## 2. SignalP fill $S"
  python3 $D/scripts/add_signalp_to_report.py $M/$S/${S}_combined_report.tsv $M/$S/prediction/interproscan/$S.tsv.gz $M/$S/${S}_combined_report.signalp.tsv || rc=1
  echo "########## 3. completeness $S"
  python3 $D/scripts/check_annotation_completeness.py --sample $S --map-dir $M/$S --funcscan-dir $F --genes-gff $RUN/gene_export/$S.gff3 || rc=1
done
echo "=== DONE $(date), status $rc (0 = all OK)"; exit $rc
