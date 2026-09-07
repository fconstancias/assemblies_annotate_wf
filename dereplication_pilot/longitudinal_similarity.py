#!/usr/bin/env python3
"""Gene-cluster similarity between every pair of a participant's own single-sample
timepoints, vs. real calendar days between them (from to_bin_sample_host_date.tsv).
Tests whether shared-gene-cluster fraction decays with time-distance -- a real
strain-turnover question, and informs whether "annotate once" stays valid over a
long study or needs periodic refresh.

Output: one row per (participant, group_i, group_j) pair, for the Rmd to plot.
"""
import csv
import glob
from collections import defaultdict
from datetime import datetime

METADATA = "/maps/projects/hansen_ol-AUDIT/scratch/NILU/metagenomes/to_bin_sample_host_date.tsv"
PILOT_DIR = "/maps/projects/hansen_ol-AUDIT/scratch/NILU/metagenomes/assembly_annotation_wf/dereplication_pilot"


def load_sample_dates():
    """spaS<n> -> (participant, datetime)"""
    dates = {}
    with open(METADATA) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            n = row["Sample"].replace("S_", "")
            group = f"spaS{n}"
            try:
                d = datetime.strptime(row["Time"], "%d/%m/%Y")
            except ValueError:
                continue
            dates[group] = (row["Subject"], d)
    return dates


def run_tier(tier, dates):
    out_rows = []

    for mfile in sorted(glob.glob(f"{PILOT_DIR}/p*/{tier}/cluster_membership.tsv")):
        participant = mfile.split("/")[-3]
        clusters_of_group = defaultdict(set)  # group -> set of cluster_ids it has a member in
        with open(mfile) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                clusters_of_group[row["source_group"]].add(row["cluster_id"])

        single_groups = sorted(g for g in clusters_of_group if not g.startswith("mh_p"))
        for i in range(len(single_groups)):
            for j in range(i + 1, len(single_groups)):
                gi, gj = single_groups[i], single_groups[j]
                if gi not in dates or gj not in dates:
                    continue
                _, di = dates[gi]
                _, dj = dates[gj]
                days = abs((dj - di).days)

                ci, cj = clusters_of_group[gi], clusters_of_group[gj]
                shared = len(ci & cj)
                union = len(ci | cj)
                jaccard = shared / union if union else 0.0

                out_rows.append((participant, gi, gj, days, shared, len(ci), len(cj), jaccard))

    out_path = f"{PILOT_DIR}/longitudinal_similarity_{tier}.tsv"
    with open(out_path, "w") as fh:
        fh.write("participant\tgroup_i\tgroup_j\tdays_between\tshared_clusters\tn_clusters_i\tn_clusters_j\tjaccard\n")
        for row in out_rows:
            fh.write("\t".join(str(x) for x in row) + "\n")

    print(f"[{tier}] wrote {len(out_rows)} pairs to {out_path}")
    if out_rows:
        days_list = [r[3] for r in out_rows]
        jac_list = [r[7] for r in out_rows]
        print(f"[{tier}] days_between range: {min(days_list)}-{max(days_list)}")
        print(f"[{tier}] jaccard range: {min(jac_list):.3f}-{max(jac_list):.3f}")


def main():
    dates = load_sample_dates()
    for tier in ("orfs", "contigs"):
        run_tier(tier, dates)


if __name__ == "__main__":
    main()
