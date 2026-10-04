#!/usr/bin/env bash
# Launch MAP v5.0.0 in <run_dir>/map with the repo template config + DB paths, driver with nohup. Same command again
# = resume. Add --preview to only compile. `run_native_signalp` (template) is only kept if the checkout knows it
# (the esrum shared checkout carries an old native-SignalP experiment; a fresh install + patches/ does not, and MAP
# fails on unknown params).
set -euo pipefail
D=$(cd "$(dirname "$0")/.." && pwd); RUN=$(readlink -f ${1:?run dir}); EXTRA=${2:-}
CHK=$HOME/.nextflow/assets/EBI-Metagenomics/mobilome-annotation-pipeline
cd $RUN/map
if [ ! -s nextflow.config ]; then
  cp $D/../templates/map.nextflow.config nextflow.config
  grep -rq run_native_signalp $CHK/workflows $CHK/subworkflows 2>/dev/null || sed -i '/run_native_signalp/d' nextflow.config
fi
source $D/scripts/00_env.sh
nohup nextflow run EBI-Metagenomics/mobilome-annotation-pipeline -r v5.0.0 --input samplesheet.csv --outdir results \
  -c nextflow.config -c $D/config/map_db_paths.config -profile singularity -resume $EXTRA > nohup_map.log 2>&1 &
echo "MAP driver PID $! -> $RUN/map/nohup_map.log"
