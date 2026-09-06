#!/usr/bin/env python3
"""Two more reports built on the same dereplication clustering, no new compute:

1. AMR redundancy: how many raw RGI hits exist across all 296 assemblies vs.
   how many *distinct* AMR genes that collapses to once cluster membership
   is accounted for (i.e., the real, non-redundant AMR gene count).
2. Co-assembly-unique contribution: for the 19 participants with a co-assembly,
   how many of its contigs/genes fall in a cluster with NO single-sample member
   at all -- i.e., only recoverable by pooling reads across timepoints, not by
   any individual single-sample assembly. Cross-referenced against AMR hits to
   flag any resistance gene the co-assembly caught that no single sample would
   have.
"""
import csv
import glob
from collections import defaultdict

from validate_propagation import (
    PILOT_DIR,
    is_coassembly,
    load_gene_calls_rgi,
)


def load_membership(path):
    cluster_of = {}
    members_of = defaultdict(list)
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            cluster_of[row["member_id"]] = row["cluster_id"]
            members_of[row["cluster_id"]].append((row["member_id"], row["source_group"]))
    return cluster_of, members_of


def amr_redundancy_report(out_path):
    lines = ["=== AMR redundancy (RGI Best_Hit_ARO), all participants ===", ""]
    total_raw = 0
    total_distinct_clusters = 0
    per_participant = []

    for mfile in sorted(glob.glob(f"{PILOT_DIR}/p*/orfs/cluster_membership.tsv")):
        participant = mfile.split("/")[-3]
        cluster_of, _ = load_membership(mfile)
        groups = set()
        with open(mfile) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                groups.add(row["source_group"])

        raw_hits = 0
        amr_clusters = set()
        for g in groups:
            calls = load_gene_calls_rgi(g)
            raw_hits += len(calls)
            for gene_id in calls:
                cid = cluster_of.get(gene_id)
                if cid is not None:
                    amr_clusters.add(cid)

        total_raw += raw_hits
        total_distinct_clusters += len(amr_clusters)
        per_participant.append((participant, raw_hits, len(amr_clusters)))

    lines.append(f"total raw RGI hits (across all groups, all participants): {total_raw}")
    lines.append(f"total distinct AMR gene clusters (deduplicated): {total_distinct_clusters}")
    if total_raw:
        lines.append(
            f"redundancy: {100.0 * (total_raw - total_distinct_clusters) / total_raw:.1f}% "
            f"of raw AMR hits are duplicates of another hit already counted"
        )
    lines.append("")
    lines.append("per-participant (raw hits, distinct AMR clusters):")
    for p, raw, distinct in per_participant:
        lines.append(f"  {p}: {raw}, {distinct}")

    report = "\n".join(lines) + "\n"
    with open(out_path, "w") as fh:
        fh.write(report)
    print(report, end="")


def coassembly_unique_report(out_path):
    lines = ["=== co-assembly-unique contigs/genes (19 participants with a co-assembly) ===", ""]

    for tier in ("contigs", "orfs"):
        lines.append(f"--- tier: {tier} ---")
        for mfile in sorted(glob.glob(f"{PILOT_DIR}/p*/{tier}/cluster_membership.tsv")):
            participant = mfile.split("/")[-3]
            cluster_of, members_of = load_membership(mfile)

            coassembly_group = None
            for cid, members in members_of.items():
                for member_id, group in members:
                    if is_coassembly(group):
                        coassembly_group = group
                        break
                if coassembly_group:
                    break
            if not coassembly_group:
                continue  # one of the 6 no-co-assembly participants

            ca_total = 0
            ca_unique = 0
            unique_amr_hits = []
            rgi_calls = load_gene_calls_rgi(coassembly_group) if tier == "orfs" else {}

            for cid, members in members_of.items():
                groups_in_cluster = set(g for _, g in members)
                ca_members = [m for m, g in members if g == coassembly_group]
                if not ca_members:
                    continue
                ca_total += len(ca_members)
                if groups_in_cluster == {coassembly_group}:
                    ca_unique += len(ca_members)
                    if tier == "orfs":
                        for m in ca_members:
                            if m in rgi_calls:
                                unique_amr_hits.append((m, rgi_calls[m]))

            pct = 100.0 * ca_unique / ca_total if ca_total else 0.0
            lines.append(
                f"  {participant} ({coassembly_group}): {ca_unique}/{ca_total} "
                f"({pct:.1f}%) co-assembly-only, no single-sample match"
            )
            if unique_amr_hits:
                for gene_id, call in unique_amr_hits:
                    lines.append(f"    AMR hit only found via co-assembly: {gene_id} -> {call}")
        lines.append("")

    report = "\n".join(lines) + "\n"
    with open(out_path, "w") as fh:
        fh.write(report)
    print(report, end="")


if __name__ == "__main__":
    import sys

    amr_redundancy_report(f"{PILOT_DIR}/amr_redundancy.txt")
    print()
    coassembly_unique_report(f"{PILOT_DIR}/coassembly_unique.txt")
