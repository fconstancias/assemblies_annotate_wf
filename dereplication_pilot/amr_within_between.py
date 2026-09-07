#!/usr/bin/env python3
"""AMR gene redundancy, properly quantified at BOTH levels:
  - within-participant (already have: report_amr_redundancy.tsv)
  - whole-cohort / between-participant: how many distinct AMR genes remain
    once cross-participant sharing (cross_participant_orfs/, the second-level
    clustering of each participant's own representatives) is also accounted
    for.

Mechanism: analyze_clusters.py sets cluster_id = the mmseqs representative ID
for that cluster, and that same representative sequence is exactly what got
pooled into cross_participant_orfs/ -- so a participant's first-level
AMR-carrying cluster_id is a valid member_id to look up in the second-level
clu_cluster.tsv directly, no extra ID translation needed.
"""
import csv
import glob
from collections import defaultdict

from validate_propagation import load_gene_calls_rgi

PILOT_DIR = "/maps/projects/hansen_ol-AUDIT/scratch/NILU/metagenomes/assembly_annotation_wf/dereplication_pilot"


def main():
    # 1. Per-participant: which first-level cluster_ids carry >=1 AMR hit?
    amr_cluster_owner = {}  # first-level cluster_id -> participant
    within_total = 0
    for mfile in sorted(glob.glob(f"{PILOT_DIR}/p*/orfs/cluster_membership.tsv")):
        participant = mfile.split("/")[-3]
        cluster_of = {}
        groups = set()
        with open(mfile) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                cluster_of[row["member_id"]] = row["cluster_id"]
                groups.add(row["source_group"])

        amr_clusters_here = set()
        for g in groups:
            for gene_id in load_gene_calls_rgi(g):
                cid = cluster_of.get(gene_id)
                if cid is not None:
                    amr_clusters_here.add(cid)

        for cid in amr_clusters_here:
            amr_cluster_owner[cid] = participant
        within_total += len(amr_clusters_here)

    # 2. Second-level clustering: first-level cluster_id (as a member) -> second-level cluster
    second_level_of = {}
    with open(f"{PILOT_DIR}/cross_participant_orfs/clu_cluster.tsv") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            rep_id, member_id = line.split("\t")
            second_level_of[member_id] = rep_id

    # 3. Whole-cohort distinct AMR genes = distinct second-level clusters touched
    second_level_participants = defaultdict(set)  # second-level cluster -> set of participants
    unresolved = 0
    for cid, participant in amr_cluster_owner.items():
        sl = second_level_of.get(cid)
        if sl is None:
            unresolved += 1
            continue
        second_level_participants[sl].add(participant)

    whole_cohort_distinct = len(second_level_participants)
    shared_across_participants = sum(1 for s in second_level_participants.values() if len(s) > 1)
    participant_specific = whole_cohort_distinct - shared_across_participants

    lines = []
    lines.append("=== AMR gene redundancy: within-participant vs. whole-cohort ===")
    lines.append(f"within-participant distinct AMR genes (summed across 25 participants): {within_total}")
    lines.append(f"whole-cohort distinct AMR genes (after cross-participant dedup too): {whole_cohort_distinct}")
    lines.append(
        f"additional collapse from cross-participant sharing: "
        f"{within_total - whole_cohort_distinct} "
        f"({100.0*(within_total - whole_cohort_distinct)/within_total:.1f}% of the within-participant count)"
    )
    lines.append("")
    lines.append(f"of the {whole_cohort_distinct} whole-cohort distinct AMR genes:")
    lines.append(
        f"  shared across >=2 participants: {shared_across_participants} "
        f"({100.0*shared_across_participants/whole_cohort_distinct:.1f}%)"
    )
    lines.append(
        f"  participant-specific (found in only 1 participant): {participant_specific} "
        f"({100.0*participant_specific/whole_cohort_distinct:.1f}%)"
    )
    if unresolved:
        lines.append(f"WARNING: {unresolved} first-level AMR clusters had no second-level match")

    # Histogram of how many participants share each whole-cohort AMR gene
    hist = defaultdict(int)
    for s in second_level_participants.values():
        hist[len(s)] += 1
    lines.append("")
    lines.append("histogram: n_participants_sharing_this_AMR_gene -> n_distinct_AMR_genes")
    for k in sorted(hist):
        lines.append(f"  {k}: {hist[k]}")

    report = "\n".join(lines) + "\n"
    with open(f"{PILOT_DIR}/amr_within_between.txt", "w") as fh:
        fh.write(report)
    print(report, end="")

    # Also write a clean TSV for the Rmd
    with open(f"{PILOT_DIR}/report_amr_within_between.tsv", "w") as fh:
        fh.write("n_participants_sharing\tn_distinct_amr_genes\n")
        for k in sorted(hist):
            fh.write(f"{k}\t{hist[k]}\n")


if __name__ == "__main__":
    main()
