#!/usr/bin/env python3
"""AMR-carrying plasmid contig clusters: within-participant vs. whole-cohort
(same two-level logic as amr_within_between.py, but for Tier-1 contig
clusters that are BOTH geNomad-plasmid-classified AND carry a real RGI AMR
hit -- the specific "AMR plasmid dynamics, within vs. between" question).
"""
import csv
from collections import defaultdict

PILOT_DIR = "/maps/projects/hansen_ol-AUDIT/scratch/NILU/metagenomes/assembly_annotation_wf/dereplication_pilot"


def main():
    # 1. First-level AMR-carrying plasmid clusters, from the already-computed join
    amr_plasmid_owner = {}  # first-level cluster_id -> participant
    with open(f"{PILOT_DIR}/plasmid_amr_join.tsv") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            if row["is_plasmid"] == "True" and row["has_amr"] == "True":
                amr_plasmid_owner[row["cluster_id"]] = row["participant"]

    within_total = len(amr_plasmid_owner)

    # 2. Second-level (cross-participant) contig clustering
    second_level_of = {}
    with open(f"{PILOT_DIR}/cross_participant/clu_cluster.tsv") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            rep_id, member_id = line.split("\t")
            second_level_of[member_id] = rep_id

    second_level_participants = defaultdict(set)
    unresolved = 0
    for cid, participant in amr_plasmid_owner.items():
        sl = second_level_of.get(cid)
        if sl is None:
            unresolved += 1
            continue
        second_level_participants[sl].add(participant)

    whole_cohort_distinct = len(second_level_participants)
    shared = sum(1 for s in second_level_participants.values() if len(s) > 1)
    specific = whole_cohort_distinct - shared

    lines = []
    lines.append("=== AMR-carrying plasmid contig clusters: within-participant vs. whole-cohort ===")
    lines.append(f"within-participant distinct AMR+plasmid clusters (summed across participants): {within_total}")
    lines.append(f"whole-cohort distinct AMR+plasmid clusters (cross-participant dedup too): {whole_cohort_distinct}")
    if within_total:
        lines.append(
            f"additional collapse from cross-participant sharing: "
            f"{within_total - whole_cohort_distinct} "
            f"({100.0*(within_total - whole_cohort_distinct)/within_total:.1f}%)"
        )
    lines.append("")
    if whole_cohort_distinct:
        lines.append(
            f"shared across >=2 participants: {shared} ({100.0*shared/whole_cohort_distinct:.1f}%)"
        )
        lines.append(
            f"participant-specific: {specific} ({100.0*specific/whole_cohort_distinct:.1f}%)"
        )
    if unresolved:
        lines.append(f"WARNING: {unresolved} first-level clusters had no second-level match")

    hist = defaultdict(int)
    for s in second_level_participants.values():
        hist[len(s)] += 1
    lines.append("")
    lines.append("histogram: n_participants_sharing_this_AMR+plasmid -> n_distinct_clusters")
    for k in sorted(hist):
        lines.append(f"  {k}: {hist[k]}")

    report = "\n".join(lines) + "\n"
    with open(f"{PILOT_DIR}/amr_plasmid_within_between.txt", "w") as fh:
        fh.write(report)
    print(report, end="")

    with open(f"{PILOT_DIR}/report_amr_plasmid_within_between.tsv", "w") as fh:
        fh.write("n_participants_sharing\tn_distinct_clusters\n")
        for k in sorted(hist):
            fh.write(f"{k}\t{hist[k]}\n")


if __name__ == "__main__":
    main()
