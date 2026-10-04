#!/usr/bin/env bash
# funcscan + MAP samplesheets from <run_dir>/gene_export (after 02). gff_type = prodigal (anvi'o's minimal GFF3 works).
set -euo pipefail
RUN=$(readlink -f ${1:?run dir}); G=$RUN/gene_export; mkdir -p $RUN/funcscan $RUN/map
echo "sample,fasta,protein,gff,gff_type" > $RUN/funcscan/samplesheet.csv
echo "sample,assembly,proteins_gff,proteins_faa,virify_gff,interproscan_tsv" > $RUN/map/samplesheet.csv
for S in $(tail -n +2 $RUN/samples.tsv | cut -f1); do
  for f in contigs.fa gff3 faa; do [ -s $G/$S.$f ] || { echo "missing $G/$S.$f (run 02 first)"; exit 1; }; done
  echo "$S,$G/$S.contigs.fa,$G/$S.faa,$G/$S.gff3,prodigal" >> $RUN/funcscan/samplesheet.csv
  echo "$S,$G/$S.contigs.fa,$G/$S.gff3,$G/$S.faa,," >> $RUN/map/samplesheet.csv
done
echo "$(($(wc -l < $RUN/map/samplesheet.csv) - 1)) samples -> $RUN/funcscan/samplesheet.csv, $RUN/map/samplesheet.csv"
