#!/usr/bin/env python3
"""Fill the empty `signalP` column of a MAP combined report from MAP's own InterProScan output.

MAP v5.0.0 runs licensed SignalP 4.1 inside InterProScan (prediction/interproscan/<sample>.tsv.gz:
SignalP_GRAM_POSITIVE / SignalP_GRAM_NEGATIVE / SignalP_EUK), but its integrator leaves `signalP`
as "-" (also in the NILU production reports). This writes a copy of the report with `signalP` =
bacterial models with a predicted signal peptide and its span, e.g.
"GRAM_POSITIVE:1-23;GRAM_NEGATIVE:1-23" ("-" if none). SignalP_EUK is ignored (bacteria).
The original report is not modified.

Usage: add_signalp_to_report.py <combined_report.tsv> <interproscan.tsv.gz> <output.tsv>
"""
import csv
import gzip
import sys
from collections import defaultdict

report, ips, out = sys.argv[1:4]

hits = defaultdict(dict)  # protein -> {model: "start-end"}
with gzip.open(ips, "rt") as fh:
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if len(f) < 8 or f[3] not in ("SignalP_GRAM_POSITIVE", "SignalP_GRAM_NEGATIVE"):
            continue
        model = f[3].replace("SignalP_", "")
        hits[f[0]].setdefault(model, f"{f[6]}-{f[7]}")

n = n_sp = 0
with open(report) as fin, open(out, "w", newline="") as fout:
    r = csv.DictReader(fin, delimiter="\t")
    w = csv.DictWriter(fout, fieldnames=r.fieldnames, delimiter="\t", lineterminator="\n")
    w.writeheader()
    for row in r:
        n += 1
        h = hits.get(row["protein_id"])
        if h:
            n_sp += 1
            row["signalP"] = ";".join(f"{m}:{h[m]}" for m in ("GRAM_POSITIVE", "GRAM_NEGATIVE") if m in h)
        w.writerow(row)

print(f"{report}: {n} rows, {n_sp} with a bacterial signal peptide; "
      f"{len(hits)} proteins with a bacterial SignalP hit in InterProScan -> {out}")
