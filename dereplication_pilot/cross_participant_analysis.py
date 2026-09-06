#!/usr/bin/env python3
"""Within- vs. between-participant sharing at the second clustering level
(cross_participant/ or cross_participant_orfs/, which cluster each
participant's own already-deduplicated representatives against everyone
else's). A cluster containing representatives from >1 participant means that
sequence is shared across different people, not just across one person's own
timepoints -- the real cross-participant redundancy question.

Also reports, for the ORF tier specifically, how many AMR-gene-carrying
representatives end up in a between-participant cluster -- i.e., is a given
resistance gene shared across different people, or person-specific?
"""
import argparse
import csv
import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(__file__))
from validate_propagation import load_gene_calls_rgi, is_coassembly

PILOT_DIR = "/maps/projects/hansen_ol-AUDIT/scratch/NILU/metagenomes/assembly_annotation_wf/dereplication_pilot"


def load_group_to_participant():
    g2p = {}
    for d in sorted(os.listdir(PILOT_DIR)):
        if not d.startswith("p") or not os.path.isdir(f"{PILOT_DIR}/{d}"):
            continue
        mfile = f"{PILOT_DIR}/{d}/contigs/cluster_membership.tsv"
        if not os.path.exists(mfile):
            continue
        participant = d
        with open(mfile) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                g2p[row["source_group"]] = participant
    return g2p


def group_of_member(member_id, groups_sorted_desc):
    for g in groups_sorted_desc:
        if member_id.startswith(g) and (len(member_id) == len(g) or member_id[len(g)] == "_"):
            return g
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tier", required=True, choices=["contigs", "orfs"])
    ap.add_argument("--cluster-tsv", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    g2p = load_group_to_participant()
    groups_sorted_desc = sorted(g2p, key=len, reverse=True)

    members_of = defaultdict(list)
    with open(args.cluster_tsv) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            rep_id, member_id = line.split("\t")
            members_of[rep_id].append(member_id)

    total_clusters = 0
    within_only = 0
    between = 0
    participant_counts_hist = defaultdict(int)
    amr_between_examples = []

    rgi_cache = {}

    for rep_id, members in members_of.items():
        total_clusters += 1
        participants = set()
        member_participants = []
        for m in members:
            g = group_of_member(m, groups_sorted_desc)
            p = g2p.get(g, "UNKNOWN")
            participants.add(p)
            member_participants.append((m, g, p))

        participant_counts_hist[len(participants)] += 1
        if len(participants) == 1:
            within_only += 1
        else:
            between += 1
            if args.tier == "orfs" and len(amr_between_examples) < 30:
                hit_members = []
                for m, g, p in member_participants:
                    if g not in rgi_cache:
                        rgi_cache[g] = load_gene_calls_rgi(g)
                    call = rgi_cache[g].get(m)
                    if call:
                        hit_members.append((p, g, m, call))
                if len(set(p for p, *_ in hit_members)) > 1:
                    amr_between_examples.append((rep_id, hit_members))

    lines = []
    lines.append(f"=== cross-participant {args.tier}: within- vs. between-participant sharing ===")
    lines.append(f"total second-level clusters: {total_clusters}")
    lines.append(
        f"within-participant only (single participant, even after cross-participant pass): "
        f"{within_only} ({100.0*within_only/total_clusters:.1f}%)"
    )
    lines.append(
        f"between-participant (shared across >=2 people): "
        f"{between} ({100.0*between/total_clusters:.1f}%)"
    )
    lines.append("")
    lines.append("histogram: n_participants_in_cluster -> n_clusters")
    for k in sorted(participant_counts_hist):
        lines.append(f"  {k}: {participant_counts_hist[k]}")

    if amr_between_examples:
        lines.append("")
        lines.append("AMR genes shared across DIFFERENT participants (real cross-person sharing, up to 30):")
        for rep_id, hit_members in amr_between_examples:
            lines.append(f"  cluster {rep_id}:")
            for p, g, m, call in hit_members:
                lines.append(f"    participant {p} ({g}:{m}) -> {call}")

    report = "\n".join(lines) + "\n"
    with open(args.out, "w") as fh:
        fh.write(report)
    print(report, end="")


if __name__ == "__main__":
    main()
