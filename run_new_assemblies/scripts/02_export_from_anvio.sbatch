#!/bin/bash
#SBATCH --account=cbmr
#SBATCH --partition=standardqueue
#SBATCH --job-name=anvio_gene_export
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=4:00:00
#SBATCH --output=anvio_gene_export_%A_%a.log
# Per sample (one array task per line of <run_dir>/samples.tsv: sample<TAB>contigs_db): contigs + gene calls from the
# anvi'o contigs DB, the input of funcscan and MAP. Protein / gene IDs <project>___<gene_callers_id>, so every
# annotation joins back to anvi'o. Read-only on the DB.
#   sbatch --array=1-$(tail -n +2 <run_dir>/samples.tsv | wc -l) scripts/02_export_from_anvio.sbatch <run_dir>
source /opt/software/mamba/23.3.1/etc/profile.d/conda.sh
conda activate anvio-9
set -euo pipefail
RUN=$(readlink -f ${1:?run dir}); i=${SLURM_ARRAY_TASK_ID:-1}
read -r S DB < <(tail -n +2 $RUN/samples.tsv | sed -n "${i}p")
O=$RUN/gene_export; mkdir -p $O; cd $O
anvi-export-contigs -c $DB -o $S.contigs.fa > /dev/null
anvi-get-sequences-for-gene-calls -c $DB --export-gff3 -o $S.gff3 > /dev/null
anvi-get-sequences-for-gene-calls -c $DB --get-aa-sequences --defline-format "{contigs_db_project_name}___{gene_caller_id}" -o $S.faa > /dev/null
anvi-get-sequences-for-gene-calls -c $DB --defline-format "{contigs_db_project_name}___{gene_caller_id}" -o $S.fna > /dev/null
n_gff=$(grep -c $'\tCDS\t' $S.gff3); n_faa=$(grep -c "^>" $S.faa); n_c=$(grep -c "^>" $S.contigs.fa)
[ "$n_gff" = "$n_faa" ] && echo "OK $S: $n_c contigs, $n_faa genes" || { echo "MISMATCH $S: gff3 $n_gff vs faa $n_faa"; exit 1; }
