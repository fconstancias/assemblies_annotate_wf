#!/usr/bin/env python3
"""Presence/absence of each whole-cohort-distinct AMR gene (the same
second-level AMR gene clusters as amr_within_between.py) across every real
sample -- both single-sample timepoints and co-assemblies -- so within-
participant temporal patterns (does a gene appear/persist/disappear over
time?) and between-participant patterns can be visualized directly, not just
summarized as a redundancy percentage.

"Present in sample X" means: X has >=1 gene belonging to the first-level
cluster whose representative maps into that whole-cohort AMR cluster -- i.e.
X carries a copy of that gene at >=95% identity/90% coverage to the AMR
reference hit, whether or not RGI happened to also call a hit specifically
on X's own copy (a looser, presence/absence-appropriate definition than the
strict per-instance validation in section 7 of the report).
"""
import csv
import glob
from collections import defaultdict

from validate_propagation import load_gene_calls_rgi
from longitudinal_similarity import load_sample_dates

PILOT_DIR = "/maps/projects/hansen_ol-AUDIT/scratch/NILU/metagenomes/assembly_annotation_wf/dereplication_pilot"


def main():
    dates = load_sample_dates()  # group -> (participant, datetime), single-sample only

    # 1. Per participant: first-level AMR clusters -> (representative ARO call, all member groups)
    amr_cluster_info = {}  # first-level cluster_id -> (participant, aro_call, set(groups))
    for mfile in sorted(glob.glob(f"{PILOT_DIR}/p*/orfs/cluster_membership.tsv")):
        participant = mfile.split("/")[-3]
        cluster_of = {}
        groups_of_cluster = defaultdict(set)
        groups = set()
        with open(mfile) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                cluster_of[row["member_id"]] = row["cluster_id"]
                groups_of_cluster[row["cluster_id"]].add(row["source_group"])
                groups.add(row["source_group"])

        amr_calls_by_cluster = defaultdict(list)
        for g in groups:
            for gene_id, call in load_gene_calls_rgi(g).items():
                cid = cluster_of.get(gene_id)
                if cid is not None:
                    amr_calls_by_cluster[cid].append(call)

        for cid, calls in amr_calls_by_cluster.items():
            amr_cluster_info[cid] = (participant, calls[0], groups_of_cluster[cid])

    # 2. Second-level (cross-participant) clustering
    second_level_of = {}
    with open(f"{PILOT_DIR}/cross_participant_orfs/clu_cluster.tsv") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            rep_id, member_id = line.split("\t")
            second_level_of[member_id] = rep_id

    # 3. Whole-cohort AMR gene -> (aro_call, set of (group, participant))
    amr_gene_presence = defaultdict(set)
    amr_gene_call = {}
    for cid, (participant, call, groups) in amr_cluster_info.items():
        sl = second_level_of.get(cid)
        if sl is None:
            continue
        for g in groups:
            amr_gene_presence[sl].add((g, participant))
        if sl not in amr_gene_call:
            amr_gene_call[sl] = call

    # 4. Write long-format presence/absence table
    out_rows = []
    for amr_gene, group_participants in amr_gene_presence.items():
        call = amr_gene_call[amr_gene]
        for g, participant in group_participants:
            is_coassembly = g.startswith("mh_p")
            date = ""
            if not is_coassembly and g in dates:
                date = dates[g][1].strftime("%Y-%m-%d")
            out_rows.append((amr_gene, call, g, participant, is_coassembly, date))

    out_path = f"{PILOT_DIR}/report_amr_presence_absence.tsv"
    with open(out_path, "w") as fh:
        fh.write("amr_gene\taro_call\tgroup\tparticipant\tis_coassembly\tdate\n")
        for row in out_rows:
            fh.write("\t".join(str(x) for x in row) + "\n")

    # 5. Prevalence summary
    prevalence = {g: len(gp) for g, gp in amr_gene_presence.items()}
    with open(f"{PILOT_DIR}/report_amr_prevalence.tsv", "w") as fh:
        fh.write("amr_gene\taro_call\tn_samples\n")
        for g, n in sorted(prevalence.items(), key=lambda x: -x[1]):
            fh.write(f"{g}\t{amr_gene_call[g]}\t{n}\n")

    print(f"wrote {len(out_rows)} presence rows, {len(amr_gene_presence)} distinct whole-cohort AMR genes")
    print(f"prevalence range: {min(prevalence.values())}-{max(prevalence.values())} samples")


if __name__ == "__main__":
    main()
