#!/usr/bin/env python3
"""Test #2 from the dereplication follow-up discussion: does the core
propagation assumption ("cluster members share the same real annotation call")
actually hold, checked against already-computed production output -- no new
annotation compute needed.

Two independent checks, using tools that only report positive hits (so the
testable set is naturally restricted to real hits, not every sequence):

- Tier 2 / genes: RGI's Best_Hit_ARO call, joined to our own gene_export IDs via
  (contig, start, stop) coordinates (RGI's own ORF_ID uses prodigal-style
  numbering, not our gene_export ID scheme, so coordinates are the only common
  key -- confirmed by inspecting both real output formats directly).
- Tier 1 / contigs: geNomad's plasmid/virus classification, joined back to our
  own contig IDs via MAP's own preprocessing/<group>_contigID.map (geNomad's
  input contigs get renamed to contig_<N> by MAP before running; the map file
  is what resolves that back).

For each cluster (from dereplication_pilot/<participant>/<tier>/cluster_membership.tsv)
with >=2 members that independently received a real hit from the tool, check
whether all those hits agree. Report the agreement rate -- this is the direct
test of whether propagating one member's annotation to the rest of its cluster
would actually be correct.
"""
import argparse
import csv
import glob
import os
import sys
from collections import defaultdict

REPO = "/maps/projects/hansen_ol-AUDIT/scratch/NILU/metagenomes/assembly_annotation_wf"
COASSEMBLY_GENE_EXPORT = f"{REPO}/coassembly_production/gene_export"
SINGLE_GENE_EXPORT = f"{REPO}/single_assembly_production/gene_export"
COASSEMBLY_FUNCSCAN = f"{REPO}/coassembly_production/funcscan_run/results"
SINGLE_FUNCSCAN = f"{REPO}/single_assembly_production/funcscan_run/results"
COASSEMBLY_MAP = f"{REPO}/coassembly_production/map_run/results"
SINGLE_MAP = f"{REPO}/single_assembly_production/map_run/results"
PILOT_DIR = f"{REPO}/dereplication_pilot"


def is_coassembly(group):
    return group.startswith("mh_p")


def gff3_path(group):
    base = COASSEMBLY_GENE_EXPORT if is_coassembly(group) else SINGLE_GENE_EXPORT
    return f"{base}/{group}.gff3"


def rgi_path(group):
    base = COASSEMBLY_FUNCSCAN if is_coassembly(group) else SINGLE_FUNCSCAN
    return f"{base}/arg/rgi/{group}/{group}.txt"


def contigid_map_path(group):
    base = COASSEMBLY_MAP if is_coassembly(group) else SINGLE_MAP
    return f"{base}/{group}/preprocessing/{group}_contigID.map"


def genomad_summary_paths(group):
    base = COASSEMBLY_MAP if is_coassembly(group) else SINGLE_MAP
    return (
        f"{base}/{group}/prediction/genomad/{group}_5kb_contigs_plasmid_summary.tsv",
        f"{base}/{group}/prediction/genomad/{group}_5kb_contigs_virus_summary.tsv",
    )


def load_gene_calls_rgi(group):
    """gene_id (our gene_export ID) -> Best_Hit_ARO, via (contig,start,stop) join."""
    gpath = gff3_path(group)
    rpath = rgi_path(group)
    if not (os.path.exists(gpath) and os.path.exists(rpath)):
        return {}

    coord_to_gene = {}
    with open(gpath) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            contig, start, stop = f[0], f[3], f[4]
            attrs = f[8]
            gene_id = None
            for kv in attrs.split(";"):
                if kv.startswith("ID="):
                    gene_id = kv[3:]
                    break
            if gene_id:
                coord_to_gene[(contig, start, stop)] = gene_id

    calls = {}
    with open(rpath) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            key = (row["Contig"], row["Start"], row["Stop"])
            gene_id = coord_to_gene.get(key)
            if gene_id:
                calls[gene_id] = row["Best_Hit_ARO"]
    return calls


def load_contig_calls_genomad(group):
    """our contig_id -> classification label ('plasmid'/'virus'), via contigID.map join."""
    mpath = contigid_map_path(group)
    if not os.path.exists(mpath):
        return {}
    genomad_to_ours = {}
    with open(mpath) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) != 2:
                continue
            genomad_id = f[0].lstrip(">")
            genomad_to_ours[genomad_id] = f[1]

    plasmid_tsv, virus_tsv = genomad_summary_paths(group)
    calls = {}
    if os.path.exists(plasmid_tsv):
        with open(plasmid_tsv) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                our_id = genomad_to_ours.get(row["seq_name"])
                if our_id:
                    calls[our_id] = "plasmid"
    if os.path.exists(virus_tsv):
        with open(virus_tsv) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                our_id = genomad_to_ours.get(row["seq_name"])
                if our_id:
                    # a contig flagged by both would be unusual; keep first (plasmid) call
                    # and note the double-hit rather than silently overwrite
                    if our_id in calls:
                        calls[our_id] = calls[our_id] + "+virus"
                    else:
                        calls[our_id] = "virus"
    return calls


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tier", required=True, choices=["orfs", "contigs"])
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    membership_files = sorted(glob.glob(f"{PILOT_DIR}/p*/{args.tier}/cluster_membership.tsv"))
    if not membership_files:
        print(f"No cluster_membership.tsv found for tier={args.tier}", file=sys.stderr)
        sys.exit(1)

    load_calls = load_gene_calls_rgi if args.tier == "orfs" else load_contig_calls_genomad
    tool_name = "RGI (Best_Hit_ARO)" if args.tier == "orfs" else "geNomad (plasmid/virus)"

    overall_testable = 0
    overall_consistent = 0
    overall_inconsistent = 0
    disagreement_examples = []
    per_participant_rows = []

    for mfile in membership_files:
        participant = mfile.split("/")[-3]
        groups_seen = set()
        cluster_members = defaultdict(list)  # cluster_id -> [(member_id, source_group)]
        with open(mfile) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                cluster_members[row["cluster_id"]].append((row["member_id"], row["source_group"]))
                groups_seen.add(row["source_group"])

        group_calls = {g: load_calls(g) for g in groups_seen}

        testable = 0
        consistent = 0
        inconsistent = 0
        for cluster_id, members in cluster_members.items():
            calls = []
            for member_id, group in members:
                c = group_calls.get(group, {}).get(member_id)
                if c is not None:
                    calls.append((member_id, group, c))
            if len(calls) >= 2:
                testable += 1
                distinct = set(c for _, _, c in calls)
                if len(distinct) == 1:
                    consistent += 1
                else:
                    inconsistent += 1
                    if len(disagreement_examples) < 25:
                        disagreement_examples.append((participant, cluster_id, calls))

        overall_testable += testable
        overall_consistent += consistent
        overall_inconsistent += inconsistent
        per_participant_rows.append((participant, testable, consistent, inconsistent))

    lines = []
    lines.append(f"=== propagation validation: {tool_name}, tier={args.tier} ===")
    lines.append(f"participants covered: {len(membership_files)}")
    lines.append(
        f"testable clusters (>=2 members independently hit by the tool): {overall_testable}"
    )
    if overall_testable:
        lines.append(
            f"consistent (all members agree): {overall_consistent} "
            f"({100.0 * overall_consistent / overall_testable:.1f}%)"
        )
        lines.append(
            f"inconsistent (members disagree): {overall_inconsistent} "
            f"({100.0 * overall_inconsistent / overall_testable:.1f}%)"
        )
    lines.append("")
    lines.append("per-participant breakdown (testable, consistent, inconsistent):")
    for participant, t, c, i in per_participant_rows:
        lines.append(f"  {participant}: {t}, {c}, {i}")

    if disagreement_examples:
        lines.append("")
        lines.append("disagreement examples (up to 25):")
        for participant, cluster_id, calls in disagreement_examples:
            lines.append(f"  [{participant}] cluster {cluster_id}:")
            for member_id, group, c in calls:
                lines.append(f"    {group}:{member_id} -> {c}")

    report = "\n".join(lines) + "\n"
    with open(args.out, "w") as fh:
        fh.write(report)
    print(report, end="")


if __name__ == "__main__":
    main()
