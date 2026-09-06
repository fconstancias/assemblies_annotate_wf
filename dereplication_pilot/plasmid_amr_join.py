#!/usr/bin/env python3
"""Contig cluster x geNomad plasmid/virus classification x AMR gene content.

For each contig CLUSTER (not raw contig -- avoids re-counting the same real
plasmid once per timepoint it happens to reassemble in), determine:
  - was any member ever classified as plasmid/virus by geNomad?
  - does any member carry >=1 real AMR gene hit (RGI), via each gene's parent
    contig (from the gff3)?
Reports the AMR-carriage rate for plasmid-classified clusters vs. everything
else, to check for real enrichment (plasmids are a classic AMR-gene
reservoir/vector) rather than assuming it.
"""
import csv
import glob
from collections import defaultdict

from validate_propagation import (
    PILOT_DIR,
    is_coassembly,
    load_contig_calls_genomad,
    rgi_path,
    gff3_path,
)


def load_contig_has_amr(group):
    """contig_id -> True if >=1 gene on it has an RGI hit."""
    gpath = gff3_path(group)
    rpath = rgi_path(group)
    if not (__import__("os").path.exists(gpath) and __import__("os").path.exists(rpath)):
        return set()

    coord_to_contig = {}
    with open(gpath) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            contig, start, stop = f[0], f[3], f[4]
            coord_to_contig[(contig, start, stop)] = contig

    contigs_with_amr = set()
    with open(rpath) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            key = (row["Contig"], row["Start"], row["Stop"])
            contig = coord_to_contig.get(key)
            if contig:
                contigs_with_amr.add(contig)
    return contigs_with_amr


def main():
    out_rows = []  # participant, cluster_id, is_plasmid, is_virus, has_amr

    for mfile in sorted(glob.glob(f"{PILOT_DIR}/p*/contigs/cluster_membership.tsv")):
        participant = mfile.split("/")[-3]
        members_of = defaultdict(list)
        groups_seen = set()
        with open(mfile) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                members_of[row["cluster_id"]].append((row["member_id"], row["source_group"]))
                groups_seen.add(row["source_group"])

        genomad_calls = {g: load_contig_calls_genomad(g) for g in groups_seen}
        amr_contigs = {g: load_contig_has_amr(g) for g in groups_seen}

        for cluster_id, members in members_of.items():
            is_plasmid = False
            is_virus = False
            has_amr = False
            for member_id, group in members:
                call = genomad_calls.get(group, {}).get(member_id)
                if call and "plasmid" in call:
                    is_plasmid = True
                if call and "virus" in call:
                    is_virus = True
                if member_id in amr_contigs.get(group, set()):
                    has_amr = True
            out_rows.append((participant, cluster_id, is_plasmid, is_virus, has_amr, len(members)))

    out_path = f"{PILOT_DIR}/plasmid_amr_join.tsv"
    with open(out_path, "w") as fh:
        fh.write("participant\tcluster_id\tis_plasmid\tis_virus\thas_amr\tn_members\n")
        for row in out_rows:
            fh.write("\t".join(str(x) for x in row) + "\n")

    total = len(out_rows)
    plasmid_clusters = [r for r in out_rows if r[2]]
    virus_clusters = [r for r in out_rows if r[3]]
    other_clusters = [r for r in out_rows if not r[2] and not r[3]]

    def amr_rate(rows):
        if not rows:
            return 0.0, 0
        n_amr = sum(1 for r in rows if r[4])
        return 100.0 * n_amr / len(rows), len(rows)

    print(f"total contig clusters (all participants): {total}")
    r, n = amr_rate(plasmid_clusters)
    print(f"plasmid-classified clusters: {n}, AMR-carrying: {r:.2f}%")
    r, n = amr_rate(virus_clusters)
    print(f"virus-classified clusters: {n}, AMR-carrying: {r:.2f}%")
    r, n = amr_rate(other_clusters)
    print(f"unclassified/other clusters: {n}, AMR-carrying: {r:.2f}%")
    print(f"wrote {out_path}")


if __name__ == "__main__":
    main()
