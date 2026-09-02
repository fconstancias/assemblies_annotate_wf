#!/bin/bash
#SBATCH --job-name=s3_sync_watcher
#SBATCH --output=/maps/projects/hansen_ol-AUDIT/scratch/NILU/metagenomes/assembly_annotation_wf/single_assembly_production/s3_sync_watcher_%j.out
#SBATCH --cpus-per-task=2
#SBATCH --mem=2G
#SBATCH --time=7-00:00:00
#SBATCH --partition=standardqueue
#SBATCH --account=cbmr

# Same pattern as ../coassembly_production/s3_sync_watcher.sh, scaled to the 277-group
# single-sample MAP run. Uses SANNTIS output as the "this group's MAP run is basically
# done" marker (BGC_ANNOTATION runs near the end of MAP's pipeline, after mobilome
# integration) -- same proxy the coassembly watcher already uses successfully.
# Runs as its own SLURM job (not a head-node background process) so it survives session
# disconnects for as long as the 277-group run takes (days, not hours).

module load rclone/1.65.1

GROUP_LIST=/maps/projects/hansen_ol-AUDIT/scratch/NILU/metagenomes/assembly_annotation_wf/single_assembly_production/group_list.txt
MAP_RESULTS=/maps/projects/hansen_ol-AUDIT/scratch/NILU/metagenomes/assembly_annotation_wf/single_assembly_production/map_run/results
MAP_DST=s3bucket:recn-fac-fbm-dbc-slehtine-stool-sampling/assembly_annotation_wf/single_assembly_production/map_run/results

STATE_DIR=/maps/projects/hansen_ol-AUDIT/scratch/NILU/metagenomes/assembly_annotation_wf/single_assembly_production/.s3_sync_state
mkdir -p "$STATE_DIR"

TOTAL=$(wc -l < "$GROUP_LIST")

while true; do
  LIVE_GROUPS=$(squeue -u "$USER" -h -o "%j" 2>/dev/null | grep -oP "\(spaS\d+\)" | tr -d '()' | sort -u)

  while read -r g; do
    sanntis="$MAP_RESULTS/$g/prediction/bgcs/sanntis/${g}_sanntis.gff.gz"
    marker="$STATE_DIR/${g}.done"
    if [ -s "$sanntis" ] && ! echo "$LIVE_GROUPS" | grep -qx "$g" && [ ! -f "$marker" ]; then
      echo "$(date +%T) $g newly complete -- uploading"
      rclone copy "$MAP_RESULTS/$g/" "$MAP_DST/$g/" --copy-links --transfers 8 --checkers 16 --checksum -v \
        >> "$STATE_DIR/${g}_upload.log" 2>&1
      if [ $? -eq 0 ]; then
        touch "$marker"
        echo "$(date +%T) UPLOADED: $g"
      else
        echo "$(date +%T) FAILED: $g upload (will retry next cycle)"
      fi
    fi
  done < "$GROUP_LIST"

  done_count=$(find "$STATE_DIR" -name "*.done" | wc -l)
  echo "$(date +%T) cycle done -- $done_count/$TOTAL groups uploaded so far"
  if [ "$done_count" -eq "$TOTAL" ]; then
    echo "$(date +%T) everything complete -- exiting watcher"
    break
  fi

  sleep 900
done
