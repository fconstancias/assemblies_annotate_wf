#!/bin/bash
#SBATCH --account=cbmr
#SBATCH --partition=standardqueue
#SBATCH --job-name=import_annotations_anvio
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=4:00:00
#SBATCH --output=import_annotations_anvio_%j.log
# TEMPLATE: funcscan + MAP results -> anvi'o functions, imported into a COPY of the contigs DB (the original is
# only read). Edit the variables below, then `sbatch import_annotations_into_anvio.sbatch`.
# Prerequisites: funcscan + MAP were run on THIS DB's gene export (protein IDs <project>___<gene_callers_id>,
# REPRODUCE.md §5); scripts/check_map_outputs.sh is clean (IntegronFinder dummies fixed, SignalP filled).
# Sources, tiers and the e_value convention (real e-value or 1 - score): header of scripts/build_anvio_functions.py.
# Checks: rows imported = rows exported per source; contigs_db_hash unchanged (profiles keep working).
# Used for matph (product test, dRep catalogue), 2026-10.
source /opt/software/mamba/23.3.1/etc/profile.d/conda.sh
conda activate anvio-9
set -euo pipefail
REPO=/maps/projects/hansen_ol-AUDIT/scratch/NILU/metagenomes/assembly_annotation_wf
DB="CHANGE_ME/project.db"                 # contigs DB (read only)
GFF="CHANGE_ME/gene_export/project.gff3"  # anvi'o gene export given to funcscan / MAP
SAMPLE="CHANGE_ME"                        # sample name in the MAP / funcscan samplesheets
MAP="CHANGE_ME/map/results/$SAMPLE"       # MAP --outdir/<sample>
FS="CHANGE_ME/funcscan/results"           # funcscan --outdir
VFDB_DMND=/maps/projects/hansen_ol-AUDIT/scratch/mobilome_db/virulence/VFDB_setB_pro.dmnd   # = MAP's --virulencefactors_db
OUT="CHANGE_ME/annotated"                 # writes $OUT/<project>_annotated.db + $OUT/functions/
case "$DB$GFF$SAMPLE$MAP$FS$OUT" in *CHANGE_ME*) echo "edit the CHANGE_ME variables first"; exit 1 ;; esac

mkdir -p $OUT; cd $OUT
VH=$OUT/vfdb_headers.txt                   # VF name + category for each VFG id
[ -s $VH ] || diamond getseq --db $VFDB_DMND | grep '^>' > $VH
python3 $REPO/scripts/build_anvio_functions.py --genes-gff $GFF --map-dir $MAP --sample $SAMPLE \
  --dbcan $FS/cazyme/dbcan/cazyme_annotation/$SAMPLE/${SAMPLE}_overview.tsv \
  --ampcombi $FS/reports/ampcombi2/$SAMPLE/${SAMPLE}_ampcombi.tsv \
  --vfdb-headers $VH --outdir functions
A=$OUT/$(basename $DB .db)_annotated.db
cp $DB $A
anvi-import-functions -c $A -i functions/functions_all.txt
anvi-export-functions -c $A -o functions/exported_all.txt
echo "== per source: in vs out"
join -1 2 -2 2 -a 2 <(awk -F'\t' 'NR>1{print $2}' functions/functions_all.txt | sort | uniq -c | sort -k2) \
                    <(awk -F'\t' 'NR>1{print $2}' functions/exported_all.txt | sort | uniq -c | sort -k2) \
  | awk '{printf "%-24s in=%-8s out=%s\n",$1,($3==""?"-":$2),($3==""?$2:$3)}'
echo "== contigs_db_hash original / annotated"
for f in $DB $A; do python3 -c "import sqlite3,sys; c=sqlite3.connect('file:'+sys.argv[1]+'?mode=ro',uri=True); print(c.execute(\"select value from self where key='contigs_db_hash'\").fetchone()[0], sys.argv[1])" $f; done
echo "=== DONE $(date)"
