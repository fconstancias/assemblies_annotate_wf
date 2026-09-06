#!/usr/bin/env python3
"""Analyze an mmseqs easy-cluster _cluster.tsv output for cross-assembly redundancy.

Parses cluster membership, attributes each member back to its source assembly
group via its header prefix (headers in this project are already globally
unique and group-prefixed, e.g. mh_p110_000000000011, spaS30___0 -- no separate
provenance tagging needed), and reports how much of a participant's own
co-assembly + single-sample contigs/genes are redundant with each other.

See ../dereplication_brainstorm.md and /home/ljc444/.claude/plans/happy-moseying-clarke.md
for why this pilot exists.
"""
import argparse
import sys
from collections import defaultdict


def read_fasta_lengths(path):
    lengths = {}
    name = None
    length = 0
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    lengths[name] = length
                name = line[1:].split()[0].strip()
                length = 0
            else:
                length += len(line.strip())
    if name is not None:
        lengths[name] = length
    return lengths


def group_of(seq_id, groups_sorted_desc):
    for g in groups_sorted_desc:
        if seq_id.startswith(g) and (len(seq_id) == len(g) or seq_id[len(g)] == "_"):
            return g
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cluster-tsv", required=True)
    ap.add_argument("--fasta", required=True, help="pooled fasta used as clustering input (for lengths)")
    ap.add_argument("--groups", nargs="+", required=True)
    ap.add_argument("--tier", required=True, choices=["contigs", "orfs"])
    ap.add_argument("--participant", required=True)
    ap.add_argument("--out", required=True, help="cluster_membership.tsv output path")
    ap.add_argument("--summary", required=True, help="summary.txt output path")
    args = ap.parse_args()

    groups_sorted_desc = sorted(set(args.groups), key=len, reverse=True)
    lengths = read_fasta_lengths(args.fasta)

    clusters = defaultdict(list)  # rep_id -> [member_id, ...]
    unresolved = []
    with open(args.cluster_tsv) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            rep_id, member_id = line.split("\t")
            clusters[rep_id].append(member_id)

    rows = []
    per_group_total = defaultdict(int)
    per_group_total_bp = defaultdict(int)
    per_group_redundant = defaultdict(int)
    per_group_redundant_bp = defaultdict(int)

    n_singleton = 0
    n_multi = 0
    n_cross_group = 0
    total_seqs = 0
    total_bp = 0
    redundant_seqs = 0
    redundant_bp = 0

    for rep_id, members in clusters.items():
        member_groups = []
        for m in members:
            g = group_of(m, groups_sorted_desc)
            if g is None:
                unresolved.append(m)
                g = "UNKNOWN"
            member_groups.append(g)

        distinct_groups = set(member_groups)
        is_cross_group = len(distinct_groups) > 1

        if len(members) == 1:
            n_singleton += 1
        else:
            n_multi += 1
        if is_cross_group:
            n_cross_group += 1

        for m, g in zip(members, member_groups):
            length = lengths.get(m, 0)
            total_seqs += 1
            total_bp += length
            per_group_total[g] += 1
            per_group_total_bp[g] += length
            if is_cross_group:
                redundant_seqs += 1
                redundant_bp += length
                per_group_redundant[g] += 1
                per_group_redundant_bp[g] += length
            rows.append((rep_id, m, g, length, is_cross_group))

    with open(args.out, "w") as fh:
        fh.write("cluster_id\tmember_id\tsource_group\tlength\tcross_group_cluster\n")
        for rep_id, m, g, length, is_cross in rows:
            fh.write(f"{rep_id}\t{m}\t{g}\t{length}\t{int(is_cross)}\n")

    lines = []
    lines.append(f"=== {args.participant} / {args.tier} ===")
    lines.append(f"total sequences pooled: {total_seqs}  ({total_bp} bp/residues)")
    lines.append(f"total clusters: {len(clusters)}  (singleton: {n_singleton}, multi-member: {n_multi})")
    lines.append(f"cross-group clusters (redundant across >=2 assemblies): {n_cross_group}")
    if total_seqs:
        lines.append(
            f"redundant sequences: {redundant_seqs}/{total_seqs} "
            f"({100.0 * redundant_seqs / total_seqs:.1f}% by count, "
            f"{100.0 * redundant_bp / total_bp:.1f}% by bp)"
        )
    lines.append("per-group breakdown (redundant/total, count and bp):")
    for g in sorted(per_group_total):
        t, r = per_group_total[g], per_group_redundant[g]
        tb, rb = per_group_total_bp[g], per_group_redundant_bp[g]
        pct = 100.0 * r / t if t else 0.0
        pctb = 100.0 * rb / tb if tb else 0.0
        lines.append(f"  {g}: {r}/{t} ({pct:.1f}%) by count, {rb}/{tb} ({pctb:.1f}%) by bp")
    if unresolved:
        lines.append(f"WARNING: {len(unresolved)} member IDs did not match any known group prefix (first 5: {unresolved[:5]})")

    summary = "\n".join(lines) + "\n"
    with open(args.summary, "w") as fh:
        fh.write(summary)
    print(summary, end="")

    if unresolved:
        sys.exit(1)


if __name__ == "__main__":
    main()
