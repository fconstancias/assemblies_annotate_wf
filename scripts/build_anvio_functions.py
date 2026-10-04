#!/usr/bin/env python3
"""funcscan + MAP outputs -> anvi'o functions-txt, one source per annotation / confidence tier.
Usage: see scripts/import_annotations_into_anvio.sbatch (and REPRODUCE.md §9). Input: MAP --outdir/<sample>,
funcscan dbCAN overview + AMPcombi table, VFDB headers (`diamond getseq --db VFDB_setB_pro.dmnd | grep '>'`).
Works for any anvi'o gene export (protein IDs <group>___<gene_callers_id>): product (spa_lib710216)
and dRep catalogue (spa_drep_reps).

Rules
- One row per gene and source; several hits are joined with '!!!' (anvi'o convention), e_value = best.
- function = clean label (gene symbol, family, VF name...) so enrichment groups identical functions;
  the evidence (tools, identity, probabilities) goes to the evidence TSVs, not into the label.
- Tiers are separate sources (e.g. AMR_all / AMR_high); thresholds are in TIERS below and in the README.
- e_value (lower = better): real e-value for Pfam / NCBIfam / InterPro / VFDB*; otherwise 1 - score:
  geNomad* = 1 - geNomad score, CheckV_quality = 1 - completeness/100, toxins = 1 - PathoFact2 toxin prob,
  AMP* = 1 - ampir prob, AMR* = 1 - best identity/100, dbCAN* = (3 - n tools)/3, BGC = 1 - SanntiS support;
  0 for SignalP and MGE_*. e.g. e_value <= 0.1 <=> score >= 0.9.

Outputs (in --outdir): functions_all.txt (for anvi-import-functions), evidence_amr.tsv,
evidence_vf_toxin.tsv, sources_summary.tsv.
"""
import argparse, csv, gzip, os, re, sys
from collections import defaultdict

TIERS = dict(
    amr_high_min_ident_multi=80.0,   # >= 2 tools AND best identity >= this
    vf_high_max_evalue=1e-10,        # VFDB hit e-value
    vf_high_min_prob=0.99,           # PathoFact2 VF probability
    tox_min_prob=0.9,                # PathoFact2 toxin probability
    amp_high_min_prob=0.9,           # ampir probability (or any DRAMP hit)
    mge_min_cds_frac=0.9,            # gene inside MGE feature (same as MAP's report)
)

ap = argparse.ArgumentParser()
ap.add_argument("--genes-gff", required=True, help="anvi'o gene export GFF3")
ap.add_argument("--map-dir", required=True, help="MAP results/<sample> dir")
ap.add_argument("--sample", required=True, help="MAP / funcscan sample name")
ap.add_argument("--dbcan", required=True, help="funcscan dbCAN *_overview.tsv")
ap.add_argument("--ampcombi", required=True, help="funcscan AMPcombi <sample>_ampcombi.tsv")
ap.add_argument("--vfdb-headers", required=True, help="VFDB fasta headers (diamond getseq | grep '>')")
ap.add_argument("--outdir", required=True)
a = ap.parse_args()
os.makedirs(a.outdir, exist_ok=True)
M, S = a.map_dir, a.sample

def gid(protein_id):
    return int(protein_id.rsplit("___", 1)[1])

def opn(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)

def clean(s):
    return re.sub(r"[\t\n\r]+", " ", str(s)).strip() or "-"

def fnum(x, default=None):
    try:
        return float(x)
    except (TypeError, ValueError):
        return default

# rows[source][gene] = list of (accession, function, e_value)
rows = defaultdict(lambda: defaultdict(list))
def add(source, g, acc, func, ev=0.0):
    rows[source][g].append((clean(acc), clean(func), 0.0 if ev is None else ev))

# ---------- genes (coordinates, for MGE overlap) ----------
genes_by_contig = defaultdict(list)
for line in open(a.genes_gff):
    if line.startswith("#"): continue
    f = line.rstrip("\n").split("\t")
    if len(f) < 9 or f[2] != "CDS": continue
    m = re.search(r"ID=([^;]+)", f[8])
    genes_by_contig[f[0]].append((int(f[3]), int(f[4]), gid(m.group(1))))
n_genes = sum(len(v) for v in genes_by_contig.values())

# ---------- CAZymes (dbCAN) ----------
for r in csv.DictReader(open(a.dbcan), delimiter="\t"):
    g, fam = gid(r["Gene ID"]), r["Recommend Results"]
    ev_t = (3 - int(r["#ofTools"])) / 3
    add("dbCAN", g, fam, fam, ev_t)
    subs = [s for s in dict.fromkeys(r["Substrate"].split(";")) if s and s != "-"]
    if subs:
        add("dbCAN_substrate", g, fam, ";".join(subs), ev_t)

# ---------- SignalP (InterProScan, all proteins) + Pfam / NCBIfam / InterPro ----------
signalp = defaultdict(set)
for line in opn(f"{M}/prediction/interproscan/{S}.tsv.gz"):
    f = line.rstrip("\n").split("\t")
    g, an = gid(f[0]), f[3]
    ev = fnum(f[8], 0.0)
    if an in ("SignalP_GRAM_POSITIVE", "SignalP_GRAM_NEGATIVE"):  # bacteria: SignalP_EUK ignored (as add_signalp_to_report.py)
        signalp[g].add(an.replace("SignalP_", ""))
    elif an in ("Pfam", "NCBIfam"):
        add(an, g, f[4], f[5] if f[5] != "-" else f[4], ev)
    if len(f) > 12 and f[11] not in ("-", ""):
        add("InterPro", g, f[11], f[12], ev)
for g, modes in signalp.items():
    add("SignalP", g, ";".join(sorted(modes)), ";".join(sorted(modes)))

# ---------- AMR (MAP: AMRFinderPlus, RGI, DeepARG; protein mode) ----------
amr = defaultdict(lambda: {"tools": set(), "ident": []})
for r in csv.DictReader(open(f"{M}/prediction/amr_genes/amrfinderplus/{S}.tsv"), delimiter="\t"):
    e = amr[gid(r["Protein id"])]; e["tools"].add("amrfinderplus")
    e["afp_symbol"], e["afp_class"], e["afp_subclass"] = r["Element symbol"], r["Class"], r["Subclass"]
    e["afp_ident"], e["afp_cov"] = r["% Identity to reference"], r["% Coverage of reference"]
    e["ident"].append(fnum(r["% Identity to reference"], 0))
for r in csv.DictReader(open(f"{M}/prediction/amr_genes/rgi/{S}.txt"), delimiter="\t"):
    e = amr[gid(r["ORF_ID"].split()[0])]; e["tools"].add("rgi")
    e["rgi_cutoff"], e["rgi_aro"], e["rgi_ident"] = r["Cut_Off"], r["Best_Hit_ARO"], r["Best_Identities"]
    e["rgi_class"], e["rgi_family"] = r["Drug Class"], r["AMR Gene Family"]
    e["ident"].append(fnum(r["Best_Identities"], 0))
for r in csv.DictReader(open(f"{M}/prediction/amr_genes/deeparg/{S}.mapping.ARG"), delimiter="\t"):
    e = amr[gid(r["read_id"])]; e["tools"].add("deeparg")
    e["deeparg_gene"], e["deeparg_class"] = r["best-hit"].split("|")[-1], r["predicted_ARG-class"]
    e["deeparg_ident"], e["deeparg_prob"] = r["identity"], r["probability"]
    e["ident"].append(fnum(r["identity"], 0))
amr_cols = ["gene_callers_id", "label", "drug_class", "n_tools", "tools", "best_identity", "AMR_high",
            "afp_symbol", "afp_class", "afp_subclass", "afp_ident", "afp_cov",
            "rgi_aro", "rgi_cutoff", "rgi_ident", "rgi_class", "rgi_family",
            "deeparg_gene", "deeparg_class", "deeparg_ident", "deeparg_prob"]
with open(f"{a.outdir}/evidence_amr.tsv", "w") as o:
    o.write("\t".join(amr_cols) + "\n")
    for g in sorted(amr):
        e = amr[g]
        label = e.get("afp_symbol") or e.get("rgi_aro") or e.get("deeparg_gene")
        dclass = (e.get("afp_class") or e.get("deeparg_class") or e.get("rgi_class", "").split(";")[0]).lower()
        best = max(e["ident"]) if e["ident"] else 0
        high = ("amrfinderplus" in e["tools"] or e.get("rgi_cutoff") == "Perfect"
                or (len(e["tools"]) >= 2 and best >= TIERS["amr_high_min_ident_multi"]))
        e.update(label=label, drug_class=dclass, n_tools=len(e["tools"]), best_identity=best, AMR_high=int(high))
        ev_a = round(1 - best / 100, 4)
        add("AMR_all", g, dclass, label, ev_a)
        if high:
            add("AMR_high", g, dclass, label, ev_a)
            add("AMR_class_high", g, dclass, dclass, ev_a)
        o.write("\t".join(clean(",".join(sorted(e["tools"])) if c == "tools" else e.get(c, "-") if c != "gene_callers_id" else g)
                          for c in amr_cols) + "\n")

# ---------- Virulence / toxins (VFDB + PathoFact2, MAP combined report) ----------
vfdb = {}  # VFG id -> (gene, vf_name, category)
for line in open(a.vfdb_headers):
    h = line[1:].strip()
    vid = h.split("(")[0]
    m = re.search(r"\[([^\[\]]+?) \((VF\d+)\) - ([^\[\]]+?) \((VFC\d+)\)\]", h)
    gm = re.match(r"\S+ \(([^)]*)\)", h)
    vfdb[vid] = (gm.group(1) if gm else "-", m.group(1) if m else "-", m.group(3) if m else "-")
vt_cols = ["gene_callers_id", "vfdb_hit", "vfdb_gene", "vf_name", "vf_category", "vfdb_evalue", "pathofact2_vf_prob",
           "pathofact2_tox_prob", "signalP", "cdd_annotation", "VF_high", "Toxin_high", "Toxin_secreted"]
with open(f"{M}/{S}_combined_report.tsv") as fh, open(f"{a.outdir}/evidence_vf_toxin.tsv", "w") as o:
    o.write("\t".join(vt_cols) + "\n")
    for r in csv.DictReader(fh, delimiter="\t"):
        g = gid(r["protein_id"])
        hit, ev = r["vfdb_hit"], fnum(r["vfdb_blastp_eval"])
        vprob, tprob = fnum(r["pathofact2_vf_prob"]), fnum(r["pathofact2_tox_prob"])
        vgene, vname, vcat = vfdb.get(hit, ("-", "-", "-")) if hit != "-" else ("-", "-", "-")
        sp = ";".join(sorted(signalp.get(g, []))) or "-"
        vf_high = hit != "-" and ev is not None and ev <= TIERS["vf_high_max_evalue"] and (vprob or 0) >= TIERS["vf_high_min_prob"]
        tox = (tprob or 0) >= TIERS["tox_min_prob"]
        tox_high, tox_sec = tox and hit != "-", tox and sp != "-"
        if hit != "-":
            add("VFDB", g, hit, vname if vname != "-" else hit, ev)
            add("VFDB_category", g, hit, vcat, ev)
        if vf_high:
            add("VF_high", g, hit, vname if vname != "-" else hit, ev)
        if tox:
            tlabel = vname if vname != "-" else (r["cdd_annotation"].split(",")[0] if r["cdd_annotation"] != "-" else "predicted toxin")
            ev_x = round(1 - tprob, 6)
            add("PathoFact2_toxin", g, hit, tlabel, ev_x)
            if tox_high: add("Toxin_high", g, hit, tlabel, ev_x)
            if tox_sec: add("Toxin_secreted", g, hit, tlabel, ev_x)
        if hit != "-" or tox or vprob is not None:
            vals = dict(gene_callers_id=g, vfdb_hit=hit, vfdb_gene=vgene, vf_name=vname, vf_category=vcat,
                        vfdb_evalue=r["vfdb_blastp_eval"], pathofact2_vf_prob=r["pathofact2_vf_prob"],
                        pathofact2_tox_prob=r["pathofact2_tox_prob"], signalP=sp, cdd_annotation=r["cdd_annotation"],
                        VF_high=int(vf_high), Toxin_high=int(tox_high), Toxin_secreted=int(tox_sec))
            o.write("\t".join(clean(vals[c]) for c in vt_cols) + "\n")

# ---------- AMPs (AMPcombi) ----------
for r in csv.DictReader(open(a.ampcombi), delimiter="\t"):
    g, p = gid(r["CDS_id"]), fnum(r["prob_ampir"], 0)
    dr = r.get("DRAMP_ID", "") not in ("", "NA", "-")
    label = clean(r["Name"]) if dr and r.get("Name") not in ("", "NA") else "AMP (ampir)"
    ev = round(1 - p, 6)
    add("AMP_all", g, r["DRAMP_ID"] if dr else "-", label, ev)
    if dr or p >= TIERS["amp_high_min_prob"]:
        add("AMP_high", g, r["DRAMP_ID"] if dr else "-", label, ev)

# ---------- BGCs (MAP: SanntiS only in this setup; antiSMASH / GECCO skipped) ----------
for line in open(f"{M}/prediction/bgcs/{S}_bgcs.gff"):
    f = line.rstrip("\n").split("\t")
    if line.startswith("#") or len(f) < 9 or f[2] != "CDS": continue
    at = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
    add("BGC_SanntiS", gid(at["ID"]), at.get("nearest_MiBIG", "-"), at.get("nearest_MiBIG_class", "-"),
        round(1 - fnum(at.get("bgc_support"), 0.0), 4))

# ---------- MGE context (gene >= 90% inside a MAP mobilome feature) ----------
skip = {"CDS", "direct_repeat_element", "inverted_repeat_element"}
for line in opn(f"{M}/gff/{S}_mobilome.gff.gz"):
    f = line.rstrip("\n").split("\t")
    if line.startswith("#") or len(f) < 9 or f[2] in skip: continue
    c, s, e = f[0], int(f[3]), int(f[4])
    met = re.search(r"mobile_element_type=([^;]+)", f[8])
    label = f[2]
    if f[2] == "insertion_sequence" and met:
        label = "IS:" + met.group(1).split("_")[0]
    for gs, ge, g in genes_by_contig.get(c, []):
        ov = min(ge, e) - max(gs, s) + 1
        if ov > 0 and ov / (ge - gs + 1) >= TIERS["mge_min_cds_frac"]:
            add("MGE_context", g, f[2], label)
            if f[2] != "compositional_outlier":
                add("MGE_strong", g, f[2], label)

# ---------- geNomad + CheckV (contigs >= 5 kb; MAP names contig_N -> anvi'o via contigID map) ----------
cmap = {}
for line in open(f"{M}/preprocessing/{S}_contigID.map"):
    f = line.lstrip(">").split()
    if len(f) >= 2: cmap[f[0]] = f[1]
checkv = {r["contig_id"]: r for r in csv.DictReader(open(f"{M}/prediction/genomad/checkv_quality_summary.tsv"), delimiter="\t")}
gcols = ["contig", "map_contig", "genomad_class", "region_start", "region_end", "length", "topology", "score",
         "n_hallmarks", "marker_enrichment", "conjugation_genes", "amr_genes", "virus_taxonomy",
         "checkv_quality", "checkv_completeness", "checkv_contamination", "miuvig_quality", "genomad_high"]
def tag_genes(contig, s, e, cls, high, extra, ev_g):
    for gs, ge, g in genes_by_contig.get(contig, []):
        if s is None or (gs >= s and ge <= e):
            add("geNomad", g, cls, cls, ev_g)
            if high: add("geNomad_high", g, cls, cls, ev_g)
            for src, val, ev_x in extra:
                if val not in ("", "-", "NA"): add(src, g, cls, val, ev_x)
with open(f"{a.outdir}/contig_genomad.tsv", "w") as o:
    o.write("\t".join(gcols) + "\n")
    for cls, fn, sc in (("plasmid", "plasmid", "plasmid_score"), ("virus", "virus", "virus_score")):
        for r in csv.DictReader(open(f"{M}/prediction/genomad/{S}_5kb_contigs_{fn}_summary.tsv"), delimiter="\t"):
            name = r["seq_name"]; mc = name.split("|")[0]; c = cmap[mc]
            s = e = None; kind = cls
            if "|provirus" in name:
                kind = "provirus"; s, e = (int(x) for x in r["coordinates"].split("-"))
            score, nh = float(r[sc]), int(r["n_hallmarks"])
            high = score >= 0.9 and (cls == "plasmid" or nh >= 1)
            tax = [t for t in r.get("taxonomy", "").split(";") if t]
            cv = checkv.get(name, {})
            vals = dict(contig=c, map_contig=mc, genomad_class=kind, region_start=s or "-", region_end=e or "-",
                        length=r["length"], topology=r["topology"], score=r[sc], n_hallmarks=nh,
                        marker_enrichment=r["marker_enrichment"], conjugation_genes=r.get("conjugation_genes", "-"),
                        amr_genes=r.get("amr_genes", "-"), virus_taxonomy=";".join(tax) or "-",
                        checkv_quality=cv.get("checkv_quality", "-"), checkv_completeness=cv.get("completeness", "-"),
                        checkv_contamination=cv.get("contamination", "-"), miuvig_quality=cv.get("miuvig_quality", "-"),
                        genomad_high=int(high))
            o.write("\t".join(clean(vals[k]) for k in gcols) + "\n")
            ev_g = round(1 - score, 6)
            ev_cv = round(1 - fnum(cv.get("completeness"), 0.0) / 100, 4)
            tag_genes(c, s, e, kind, high, [("geNomad_virus_taxonomy", tax[-1] if (cls == "virus" and tax) else "-", ev_g),
                                            ("CheckV_quality", cv.get("checkv_quality", "-") if cls == "virus" else "-", ev_cv)], ev_g)

# ---------- write ----------
with open(f"{a.outdir}/functions_all.txt", "w") as o, open(f"{a.outdir}/sources_summary.tsv", "w") as sm:
    o.write("gene_callers_id\tsource\taccession\tfunction\te_value\n")
    sm.write("source\tn_genes\n")
    for src in sorted(rows):
        for g in sorted(rows[src]):
            pairs = list(dict.fromkeys((h[0], h[1]) for h in rows[src][g]))  # accession/function stay paired
            acc = "!!!".join(p[0] for p in pairs)
            fun = "!!!".join(p[1] for p in pairs)
            ev = min(h[2] for h in rows[src][g])
            o.write(f"{g}\t{src}\t{acc}\t{fun}\t{ev:g}\n")
        sm.write(f"{src}\t{len(rows[src])}\n")
print(f"{n_genes} genes in GFF; sources:")
print(open(f"{a.outdir}/sources_summary.tsv").read())
