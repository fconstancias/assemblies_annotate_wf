#!/usr/bin/env bash
# Post-run check of MAP (mobilome-annotation-pipeline v5.0.0) results for the SILENT failures seen so far:
# tasks that exit 0 but publish empty / dummy / incomplete output. Run after every MAP run, before using it.
#
#   scripts/check_map_outputs.sh <MAP --outdir> [<MAP work dir>]
#
# Per sample (= each <outdir>/<sample>/ with a prediction/ folder):
#   combined report      present and non-empty
#   signalP column       always "-" in MAP v5.0.0 -> fill with scripts/add_signalp_to_report.py
#   PathoFact2 gff       must have feature lines (patch pathofact2_integrator_200line_cap; see README Patches)
#   IntegronFinder       contig_dummy.* = MAP's wrapper bug (patch 5) unless IntegronFinder really found nothing;
#                        with the work dir given, looks for the real Results_Integron_Finder_<sample>_* there
#   geNomad / CheckV     summaries present
#   InterProScan         tsv present (needed for SignalP fill, SanntiS)
# Prints OK / WARN / FAIL lines; exit 1 if any FAIL.
set -uo pipefail
OUT=${1:?usage: check_map_outputs.sh <MAP outdir> [<MAP work dir>]}; WORK=${2:-}
fail=0
for d in "$OUT"/*/; do
  s=$(basename "$d"); [ -d "$d/prediction" ] || continue
  echo "== $s"
  r="$d/${s}_combined_report.tsv"
  if [ -s "$r" ]; then echo "OK   combined report: $(($(wc -l < "$r") - 1)) proteins"
  else echo "FAIL combined report missing/empty: $r"; fail=1; fi
  if [ -s "$r" ]; then
    c=$(head -1 "$r" | tr '\t' '\n' | grep -n '^signalP$' | cut -d: -f1)
    n=$(awk -F'\t' -v c="$c" 'NR>1 && $c!="-" && $c!=""' "$r" | wc -l)
    if [ -s "${r%.tsv}.signalp.tsv" ]; then echo "OK   SignalP-filled report present (${s}_combined_report.signalp.tsv)"
    elif [ "$n" -eq 0 ]; then echo "WARN signalP column empty -> python3 scripts/add_signalp_to_report.py $r $d/prediction/interproscan/$s.tsv.gz ${r%.tsv}.signalp.tsv"
    else echo "OK   signalP column: $n proteins"; fi
  fi
  g="$d/prediction/virulence/${s}_pathofact2.gff"
  if [ -s "$g" ] && grep -qv '^#' "$g"; then echo "OK   PathoFact2 gff: $(grep -vc '^#' "$g") features"
  else echo "FAIL PathoFact2 gff missing or without features (200-line-cap patch applied?): $g"; fail=1; fi
  i="$d/prediction/integronfinder"
  if [ -e "$i/contig_dummy.summary" ]; then
    msg="WARN IntegronFinder: dummy outputs (MAP wrapper bug, patch 5, or really no integron element)"
    if [ -n "$WORK" ]; then
      real=$(find "$WORK" -maxdepth 3 -type d -name "Results_Integron_Finder_${s}_100kb_contigs" 2>/dev/null | head -1)
      if [ -n "$real" ]; then
        ng=$(find "$real" -maxdepth 1 -name '*.gbk' | wc -l)
        [ "$ng" -gt 0 ] && msg="FAIL IntegronFinder: dummy published but $ng .gbk in $real -> copy *.gbk *.summary *.integrons back (README Patches 5)" && fail=1
      fi
    fi
    echo "$msg"
  elif compgen -G "$i/*.summary" > /dev/null; then
    echo "OK   IntegronFinder: $(ls "$i"/*.gbk 2>/dev/null | wc -l) gbk, $(grep -v '^#' "$i"/*.summary | awk -F'\t' 'NR>1{c+=$2;k+=$3;z+=$4} END{print "CALIN "c+0", complete "k+0", In0 "z+0}')"
  else echo "WARN IntegronFinder: no output in $i"; fi
  for f in "${s}_5kb_contigs_plasmid_summary.tsv" "${s}_5kb_contigs_virus_summary.tsv" checkv_quality_summary.tsv; do
    [ -s "$d/prediction/genomad/$f" ] && echo "OK   genomad/$f: $(($(wc -l < "$d/prediction/genomad/$f") - 1)) rows" || echo "WARN genomad/$f missing"
  done
  [ -s "$d/prediction/interproscan/$s.tsv.gz" ] && echo "OK   InterProScan tsv" || echo "WARN InterProScan tsv missing (SignalP fill / SanntiS need it)"
done
exit $fail
