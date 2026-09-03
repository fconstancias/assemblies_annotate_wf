#!/bin/bash
#SBATCH --job-name=funcscan_s3_sync_watcher
#SBATCH --output=/maps/projects/hansen_ol-AUDIT/scratch/NILU/metagenomes/assembly_annotation_wf/single_assembly_production/funcscan_s3_sync_watcher_%j.out
#SBATCH --cpus-per-task=2
#SBATCH --mem=2G
#SBATCH --time=7-00:00:00
#SBATCH --partition=standardqueue
#SBATCH --account=cbmr

# Companion to ../single_assembly_production/s3_sync_watcher.sh (that one covers MAP).
# funcscan's results/ layout is tool-first (results/<category>/<tool>/<sample>/...), not
# sample-first like MAP's results/<sample>/... -- confirmed across all 8 tool dirs
# (arg/{rgi,amrfinderplus,abricate,deeparg,fargene,hamronization}, amp/ampir,
# cazyme/dbcan/{cazyme_annotation,cgc,substrate}) in the coassembly run this is copied
# from. Slicing a granular per-sample upload across that many tool dirs is fragile and
# not worth it -- instead this does a periodic whole-tree `rclone copy --checksum`
# (idempotent, only transfers new/changed files, same as every other watcher this
# project uses) gated by counting how many samples have reached dbCAN's substrate step
# -- the last step of RUNDBCAN, funcscan's slowest tool chain here (see
# funcscan_run/nextflow.config's time/memory override comment) -- as the "most samples
# are basically done" proxy for when a sync pass is worth running.

module load rclone/1.65.1

GROUP_LIST=/maps/projects/hansen_ol-AUDIT/scratch/NILU/metagenomes/assembly_annotation_wf/single_assembly_production/group_list.txt
FUNCSCAN_RESULTS=/maps/projects/hansen_ol-AUDIT/scratch/NILU/metagenomes/assembly_annotation_wf/single_assembly_production/funcscan_run/results
FUNCSCAN_DST=s3bucket:recn-fac-fbm-dbc-slehtine-stool-sampling/assembly_annotation_wf/single_assembly_production/funcscan_run/results

STATE_DIR=/maps/projects/hansen_ol-AUDIT/scratch/NILU/metagenomes/assembly_annotation_wf/single_assembly_production/.s3_sync_state
mkdir -p "$STATE_DIR"

TOTAL=$(wc -l < "$GROUP_LIST")

while true; do
  done_count=0
  while read -r g; do
    marker="$FUNCSCAN_RESULTS/cazyme/dbcan/substrate/$g/${g}_substrate_prediction.tsv"
    [ -s "$marker" ] && done_count=$((done_count+1))
  done < "$GROUP_LIST"

  echo "$(date +%T) $done_count/$TOTAL samples reached dbCAN substrate step -- running whole-tree sync"
  rclone copy "$FUNCSCAN_RESULTS/" "$FUNCSCAN_DST/" --copy-links --transfers 8 --checkers 16 --checksum -v \
    >> "$STATE_DIR/funcscan_upload.log" 2>&1
  if [ $? -eq 0 ]; then
    echo "$(date +%T) sync cycle OK"
  else
    echo "$(date +%T) sync cycle FAILED (will retry next cycle)"
  fi

  if [ "$done_count" -eq "$TOTAL" ]; then
    echo "$(date +%T) all $TOTAL samples reached dbCAN substrate step and final sync ran -- exiting watcher"
    break
  fi

  sleep 900
done
