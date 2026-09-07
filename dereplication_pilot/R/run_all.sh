#!/bin/bash
#SBATCH --job-name=derep_R_port
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=4:00:00
#SBATCH --partition=standardqueue
#SBATCH --account=cbmr
#SBATCH --output=/maps/projects/hansen_ol-AUDIT/scratch/NILU/metagenomes/assembly_annotation_wf/dereplication_pilot/logs/%x_%j.log

set -euo pipefail
source /opt/software/mamba/23.3.1/etc/profile.d/conda.sh
conda activate r-binner-compare
cd /maps/projects/hansen_ol-AUDIT/scratch/NILU/metagenomes/assembly_annotation_wf/dereplication_pilot
Rscript R/run_all.R
