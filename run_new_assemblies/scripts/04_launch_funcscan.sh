#!/usr/bin/env bash
# Launch nf-core/funcscan 4.0.0 (ARG + CAZyme + AMP) in <run_dir>/funcscan with the repo template config, driver with
# nohup on the login node (tasks go to SLURM). Same command again = resume. Add --preview to only compile.
set -euo pipefail
D=$(cd "$(dirname "$0")/.." && pwd); RUN=$(readlink -f ${1:?run dir}); EXTRA=${2:-}
cd $RUN/funcscan; [ -s nextflow.config ] || cp $D/../templates/funcscan.nextflow.config nextflow.config
source $D/scripts/00_env.sh
nohup nextflow run nf-core/funcscan -r 4.0.0 --input samplesheet.csv --outdir results \
  --run_arg_screening --run_cazyme_screening --run_amp_screening --amp_skip_amplify --amp_skip_macrel \
  -c nextflow.config -profile singularity -resume $EXTRA > nohup_funcscan.log 2>&1 &
echo "funcscan driver PID $! -> $RUN/funcscan/nohup_funcscan.log"
