#!/usr/bin/env python3
"""Consolidate every per-participant summary.txt + the various analysis
outputs into clean, machine-readable TSVs for the Rmd report to read directly
-- avoids parsing free-text prose in R.
"""
import csv
import glob
import re
from collections import defaultdict

PILOT_DIR = "/maps/projects/hansen_ol-AUDIT/scratch/NILU/metagenomes/assembly_annotation_wf/dereplication_pilot"


def parse_summary(path):
    """Returns (total_seqs, total_bp, total_clusters, singleton, multi,
    cross_group_clusters, redundant_count, redundant_count_pct,
    redundant_bp_pct, per_group rows)."""
    text = open(path).read()
    total_seqs, total_bp = map(
        int, re.search(r"total sequences pooled: (\d+)\s+\((\d+) bp", text).groups()
    )
    total_clusters, singleton, multi = map(
        int,
        re.search(
            r"total clusters: (\d+)\s+\(singleton: (\d+), multi-member: (\d+)\)", text
        ).groups(),
    )
    cross_group = int(re.search(r"cross-group clusters.*?: (\d+)", text).group(1))
    m = re.search(
        r"redundant sequences: (\d+)/(\d+) \(([\d.]+)% by count, ([\d.]+)% by bp\)", text
    )
    redundant_count, _, redundant_pct, redundant_bp_pct = m.groups()

    per_group = []
    for line in text.splitlines():
        gm = re.match(
            r"\s+(\S+): (\d+)/(\d+) \(([\d.]+)%\) by count, (\d+)/(\d+) \(([\d.]+)%\) by bp",
            line,
        )
        if gm:
            g, r_c, t_c, pct_c, r_b, t_b, pct_b = gm.groups()
            per_group.append((g, int(r_c), int(t_c), float(pct_c), int(r_b), int(t_b), float(pct_b)))

    return {
        "total_seqs": total_seqs,
        "total_bp": total_bp,
        "total_clusters": total_clusters,
        "singleton": singleton,
        "multi": multi,
        "cross_group_clusters": cross_group,
        "redundant_count": int(redundant_count),
        "redundant_pct": float(redundant_pct),
        "redundant_bp_pct": float(redundant_bp_pct),
        "per_group": per_group,
    }


def main():
    # 1. Per-participant, per-tier overview
    overview_rows = []
    per_group_rows = []
    for tier in ("contigs", "orfs"):
        for f in sorted(glob.glob(f"{PILOT_DIR}/p*/{tier}/summary.txt")):
            participant = f.split("/")[-3]
            d = parse_summary(f)
            overview_rows.append(
                (
                    participant,
                    tier,
                    d["total_seqs"],
                    d["total_bp"],
                    d["total_clusters"],
                    d["singleton"],
                    d["multi"],
                    d["redundant_count"],
                    d["redundant_pct"],
                    d["redundant_bp_pct"],
                )
            )
            for g, r_c, t_c, pct_c, r_b, t_b, pct_b in d["per_group"]:
                is_coassembly = g.startswith("mh_p")
                per_group_rows.append(
                    (participant, tier, g, is_coassembly, r_c, t_c, pct_c, r_b, t_b, pct_b)
                )

    with open(f"{PILOT_DIR}/report_overview.tsv", "w") as fh:
        fh.write(
            "participant\ttier\ttotal_seqs\ttotal_bp\ttotal_clusters\tsingleton\tmulti\t"
            "redundant_count\tredundant_pct\tredundant_bp_pct\n"
        )
        for row in overview_rows:
            fh.write("\t".join(str(x) for x in row) + "\n")

    with open(f"{PILOT_DIR}/report_per_group.tsv", "w") as fh:
        fh.write(
            "participant\ttier\tgroup\tis_coassembly\tredundant_count\ttotal_count\t"
            "redundant_pct\tredundant_bp\ttotal_bp\tredundant_bp_pct\n"
        )
        for row in per_group_rows:
            fh.write("\t".join(str(x) for x in row) + "\n")

    # 2. Sensitivity sweep (contigs cov-mode, orfs min-seq-id) -- both original + sweep summaries
    sweep_rows = []
    for participant in ("p110", "p97"):
        d0 = parse_summary(f"{PILOT_DIR}/{participant}/contigs/summary.txt")
        sweep_rows.append((participant, "contigs", "cov-mode 2 (original)", d0["redundant_pct"], d0["redundant_bp_pct"]))
        d1 = parse_summary(f"{PILOT_DIR}/{participant}/contigs/sweep_covmode0/summary.txt")
        sweep_rows.append((participant, "contigs", "cov-mode 0", d1["redundant_pct"], d1["redundant_bp_pct"]))
        d2 = parse_summary(f"{PILOT_DIR}/{participant}/orfs/summary.txt")
        sweep_rows.append((participant, "orfs", "min-seq-id 0.95 (original)", d2["redundant_pct"], d2["redundant_bp_pct"]))
        d3 = parse_summary(f"{PILOT_DIR}/{participant}/orfs/sweep_minid099/summary.txt")
        sweep_rows.append((participant, "orfs", "min-seq-id 0.99", d3["redundant_pct"], d3["redundant_bp_pct"]))

    with open(f"{PILOT_DIR}/report_sweep.tsv", "w") as fh:
        fh.write("participant\ttier\tsetting\tredundant_pct\tredundant_bp_pct\n")
        for row in sweep_rows:
            fh.write("\t".join(str(x) for x in row) + "\n")

    # 3. Propagation validation (RGI + geNomad) per-participant
    def parse_validation(path):
        text = open(path).read()
        rows = []
        in_block = False
        for line in text.splitlines():
            if line.startswith("per-participant breakdown"):
                in_block = True
                continue
            if in_block:
                m = re.match(r"\s+(\S+): (\d+), (\d+), (\d+)", line)
                if m:
                    p, t, c, i = m.groups()
                    rows.append((p, int(t), int(c), int(i)))
                elif line.strip() == "" or line.startswith("disagreement"):
                    break
        return rows

    with open(f"{PILOT_DIR}/report_validation.tsv", "w") as fh:
        fh.write("tool\tparticipant\ttestable\tconsistent\tinconsistent\n")
        for p, t, c, i in parse_validation(f"{PILOT_DIR}/validate_orfs_rgi.txt"):
            fh.write(f"RGI\t{p}\t{t}\t{c}\t{i}\n")
        for p, t, c, i in parse_validation(f"{PILOT_DIR}/validate_contigs_genomad.txt"):
            fh.write(f"geNomad\t{p}\t{t}\t{c}\t{i}\n")

    # 4. AMR redundancy per-participant
    text = open(f"{PILOT_DIR}/amr_redundancy.txt").read()
    with open(f"{PILOT_DIR}/report_amr_redundancy.tsv", "w") as fh:
        fh.write("participant\traw_hits\tdistinct_clusters\n")
        in_block = False
        for line in text.splitlines():
            if line.startswith("per-participant"):
                in_block = True
                continue
            if in_block:
                m = re.match(r"\s+(\S+): (\d+), (\d+)", line)
                if m:
                    fh.write(f"{m.group(1)}\t{m.group(2)}\t{m.group(3)}\n")

    # 5. Co-assembly unique, per-participant per-tier + AMR examples
    text = open(f"{PILOT_DIR}/coassembly_unique.txt").read()
    coassembly_rows = []
    amr_unique_rows = []
    current_tier = None
    current_participant = None
    for line in text.splitlines():
        tm = re.match(r"--- tier: (\w+) ---", line)
        if tm:
            current_tier = tm.group(1)
            continue
        pm = re.match(r"\s+(\S+) \((\S+)\): (\d+)/(\d+) \(([\d.]+)%\)", line)
        if pm:
            p, ca_group, uniq, total, pct = pm.groups()
            current_participant = p
            coassembly_rows.append((p, ca_group, current_tier, int(uniq), int(total), float(pct)))
            continue
        am = re.match(r"\s+AMR hit only found via co-assembly: (\S+) -> (.+)", line)
        if am:
            gene_id, call = am.groups()
            amr_unique_rows.append((current_participant, gene_id, call))

    with open(f"{PILOT_DIR}/report_coassembly_unique.tsv", "w") as fh:
        fh.write("participant\tcoassembly_group\ttier\tunique_count\ttotal_count\tunique_pct\n")
        for row in coassembly_rows:
            fh.write("\t".join(str(x) for x in row) + "\n")

    with open(f"{PILOT_DIR}/report_amr_coassembly_unique.tsv", "w") as fh:
        fh.write("participant\tgene_id\taro_call\n")
        for row in amr_unique_rows:
            fh.write("\t".join(str(x) for x in row) + "\n")

    # 6. Cross-participant within/between (contigs + orfs) histograms
    def parse_cross(path):
        text = open(path).read()
        total = int(re.search(r"total second-level clusters: (\d+)", text).group(1))
        within = int(re.search(r"within-participant only.*?: (\d+)", text).group(1))
        between = int(re.search(r"between-participant.*?: (\d+)", text).group(1))
        hist = {}
        in_hist = False
        for line in text.splitlines():
            if line.startswith("histogram"):
                in_hist = True
                continue
            if in_hist:
                m = re.match(r"\s+(\d+): (\d+)", line)
                if m:
                    hist[int(m.group(1))] = int(m.group(2))
                elif line.strip() == "":
                    break
        return total, within, between, hist

    with open(f"{PILOT_DIR}/report_cross_participant.tsv", "w") as fh:
        fh.write("tier\tn_participants_in_cluster\tn_clusters\n")
        for tier, path in (
            ("contigs", f"{PILOT_DIR}/cross_participant_contigs_analysis.txt"),
            ("orfs", f"{PILOT_DIR}/cross_participant_orfs_analysis.txt"),
        ):
            total, within, between, hist = parse_cross(path)
            for k, v in sorted(hist.items()):
                fh.write(f"{tier}\t{k}\t{v}\n")

    # 7. Plasmid x AMR -- already a clean TSV (plasmid_amr_join.tsv), just summarize
    plasmid_amr_summary = []
    rows = list(csv.DictReader(open(f"{PILOT_DIR}/plasmid_amr_join.tsv"), delimiter="\t"))
    for category, pred in (
        ("plasmid", lambda r: r["is_plasmid"] == "True"),
        ("virus", lambda r: r["is_virus"] == "True" and r["is_plasmid"] != "True"),
        ("other", lambda r: r["is_plasmid"] != "True" and r["is_virus"] != "True"),
    ):
        sub = [r for r in rows if pred(r)]
        n_amr = sum(1 for r in sub if r["has_amr"] == "True")
        plasmid_amr_summary.append((category, len(sub), n_amr))

    with open(f"{PILOT_DIR}/report_plasmid_amr.tsv", "w") as fh:
        fh.write("category\tn_clusters\tn_amr_carrying\n")
        for row in plasmid_amr_summary:
            fh.write("\t".join(str(x) for x in row) + "\n")

    print("Wrote all report_*.tsv files to", PILOT_DIR)


if __name__ == "__main__":
    main()
