# Source (don't run) in the SAME shell as every `nextflow run`:  source scripts/00_env.sh
source /opt/software/mamba/23.3.1/etc/profile.d/conda.sh
source /usr/share/Modules/init/bash
conda activate env_nf                  # Nextflow >= 25.10.4 (funcscan 4.0.0); `nextflow self-update` if older
module load singularity/3.8.7
export TMPDIR=/tmp                     # an inherited TMPDIR the containers can't see breaks AMRFinderPlus etc.
export NXF_OPTS="-Xms1g -Xmx4g"
