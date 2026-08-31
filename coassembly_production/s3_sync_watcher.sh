#!/bin/bash
#SBATCH --job-name=s3_sync_watcher
#SBATCH --output=/maps/projects/hansen_ol-AUDIT/scratch/NILU/metagenomes/assembly_annotation_wf/coassembly_production/s3_sync_watcher_%j.out
#SBATCH --cpus-per-task=2
#SBATCH --mem=2G
#SBATCH --time=7-00:00:00
#SBATCH --partition=standardqueue
#SBATCH --account=cbmr

# Periodically checks which co-assembly groups have newly finished MAP's BGC/SANNTIS
# step or funcscan's CAZyme substrate step (the only pieces still held back from the
# earlier S3 uploads), and uploads exactly those groups' newly-complete output the
# moment they're done -- without re-uploading anything already on S3 (rclone copy
# with --checksum only transfers new/changed files, so re-running is always safe/cheap).
# Runs as its own SLURM job (not a head-node background process) since it needs to
# survive for as long as the production run does -- head-node processes get killed
# after ~15-17 min regardless of what they're doing (confirmed earlier this session).

module load rclone/1.65.1

ALL_GROUPS="mh_p110 mh_p125 mh_p292 mh_p341 mh_p398 mh_p408 mh_p470 mh_p523 mh_p550 mh_p584 mh_p601 mh_p644 mh_p722 mh_p728 mh_p789 mh_p813 mh_p860 mh_p894 mh_p97"
MAP_RESULTS=/maps/projects/hansen_ol-AUDIT/scratch/NILU/metagenomes/assembly_annotation_wf/coassembly_production/map_run/results
FUNCSCAN_RESULTS=/maps/projects/hansen_ol-AUDIT/scratch/NILU/metagenomes/assembly_annotation_wf/coassembly_production/funcscan_run/results
MAP_DST=s3bucket:recn-fac-fbm-dbc-slehtine-stool-sampling/assembly_annotation_wf/coassembly_production/map_run/results
FUNCSCAN_DST=s3bucket:recn-fac-fbm-dbc-slehtine-stool-sampling/assembly_annotation_wf/coassembly_production/funcscan_run/results

STATE_DIR=/maps/projects/hansen_ol-AUDIT/scratch/NILU/metagenomes/assembly_annotation_wf/coassembly_production/.s3_sync_state
mkdir -p "$STATE_DIR"

while true; do
  LIVE_GROUPS=$(squeue -u "$USER" -h -o "%j" 2>/dev/null | grep -oP "\(mh_p\d+\)" | tr -d '()' | sort -u)

  for g in $ALL_GROUPS; do
    # MAP BGC/SANNTIS: upload the whole group once sanntis output exists and no job
    # is currently touching this group (avoids picking up a mid-rewrite retry).
    sanntis="$MAP_RESULTS/$g/prediction/bgcs/sanntis/${g}_sanntis.gff.gz"
    marker="$STATE_DIR/map_bgc_${g}.done"
    if [ -s "$sanntis" ] && ! echo "$LIVE_GROUPS" | grep -qx "$g" && [ ! -f "$marker" ]; then
      echo "$(date +%T) MAP BGC newly complete for $g -- uploading"
      rclone copy "$MAP_RESULTS/$g/" "$MAP_DST/$g/" --copy-links --transfers 8 --checkers 16 --checksum -v \
        >> "$STATE_DIR/map_${g}_upload.log" 2>&1
      if [ $? -eq 0 ]; then
        touch "$marker"
        echo "$(date +%T) UPLOADED: MAP BGC for $g"
      else
        echo "$(date +%T) FAILED: MAP BGC upload for $g (will retry next cycle)"
      fi
    fi
  done

  # funcscan CAZyme substrate: only mh_p789 and mh_p813 were held back.
  for g in mh_p789 mh_p813; do
    marker="$STATE_DIR/funcscan_substrate_${g}.done"
    if ! echo "$LIVE_GROUPS" | grep -qx "$g" && [ ! -f "$marker" ]; then
      echo "$(date +%T) funcscan CAZyme substrate newly complete for $g -- uploading"
      rclone copy "$FUNCSCAN_RESULTS/cazyme/dbcan/substrate/$g/" "$FUNCSCAN_DST/cazyme/dbcan/substrate/$g/" \
        --copy-links --transfers 8 --checkers 16 --checksum -v \
        >> "$STATE_DIR/funcscan_${g}_upload.log" 2>&1
      if [ $? -eq 0 ]; then
        touch "$marker"
        echo "$(date +%T) UPLOADED: funcscan CAZyme substrate for $g"
      else
        echo "$(date +%T) FAILED: funcscan substrate upload for $g (will retry next cycle)"
      fi
    fi
  done

  done_count=$(find "$STATE_DIR" -name "map_bgc_*.done" | wc -l)
  echo "$(date +%T) cycle done -- $done_count/19 MAP BGC groups uploaded so far"
  if [ "$done_count" -eq 19 ] && [ -f "$STATE_DIR/funcscan_substrate_mh_p789.done" ] && [ -f "$STATE_DIR/funcscan_substrate_mh_p813.done" ]; then
    echo "$(date +%T) everything complete -- exiting watcher"
    break
  fi

  sleep 900
done
