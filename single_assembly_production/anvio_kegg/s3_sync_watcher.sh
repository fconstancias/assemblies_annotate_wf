#!/bin/bash
#SBATCH --job-name=s3_sync_watcher_kegg
#SBATCH --output=/maps/projects/hansen_ol-AUDIT/scratch/NILU/metagenomes/spa_single_all_anvio_kegg/s3_sync_watcher_%j.out
#SBATCH --cpus-per-task=2
#SBATCH --mem=2G
#SBATCH --time=7-00:00:00
#SBATCH --partition=standardqueue
#SBATCH --account=cbmr

# Periodically checks which single-assembly groups have finished BOTH
# anvi-run-kegg-kofams and anvi-estimate-metabolism (marker: {group}_metabolism_modules.txt
# exists, since anvi-estimate-metabolism only writes it on successful completion), and
# uploads that group's contigs-db + metabolism file to S3 the moment it's done. Runs as
# its own SLURM job (not a head-node process) since it needs to survive the full
# multi-day kofams array run.

module load rclone/1.65.1

GROUP_LIST=/maps/projects/hansen_ol-AUDIT/scratch/NILU/metagenomes/spa_single_all_anvio_kegg/group_list.txt
KEGG_DIR=/maps/projects/hansen_ol-AUDIT/scratch/NILU/metagenomes/spa_single_all_anvio_kegg
DST=s3bucket:recn-fac-fbm-dbc-slehtine-stool-sampling/spa_single_all_anvio_kegg
STATE_DIR=$KEGG_DIR/.s3_sync_state
mkdir -p "$STATE_DIR"

while true; do
  # Map currently-running array task IDs -> group names, so we never upload a group
  # whose kofams/estimate-metabolism job is still (re)writing it.
  LIVE_TASK_IDS=$(squeue -u "$USER" -h -t RUNNING -n kofams_single -o "%K" 2>/dev/null)
  LIVE_GROUPS=""
  for t in $LIVE_TASK_IDS; do
    g=$(sed -n "$((t + 1))p" "$GROUP_LIST")
    LIVE_GROUPS="$LIVE_GROUPS $g"
  done

  n_uploaded=0
  while read -r g; do
    marker="$STATE_DIR/${g}.done"
    modules_file="$KEGG_DIR/${g}_metabolism_modules.txt"
    [ -f "$marker" ] && { n_uploaded=$((n_uploaded + 1)); continue; }
    if [ -s "$modules_file" ] && ! echo " $LIVE_GROUPS " | grep -q " $g "; then
      echo "$(date +%T) $g newly complete -- uploading"
      rclone copy "$KEGG_DIR/contigs_dbs/${g}.db" "$DST/contigs_dbs/" --checksum -v >> "$STATE_DIR/${g}_upload.log" 2>&1
      rclone copy "$modules_file" "$DST/" --checksum -v >> "$STATE_DIR/${g}_upload.log" 2>&1
      if [ $? -eq 0 ]; then
        touch "$marker"
        n_uploaded=$((n_uploaded + 1))
        echo "$(date +%T) UPLOADED: $g"
      else
        echo "$(date +%T) FAILED: $g (will retry next cycle)"
      fi
    fi
  done < "$GROUP_LIST"

  echo "$(date +%T) cycle done -- $n_uploaded/277 groups uploaded so far"
  if [ "$n_uploaded" -eq 277 ]; then
    echo "$(date +%T) all groups complete -- exiting watcher"
    break
  fi

  sleep 900
done
