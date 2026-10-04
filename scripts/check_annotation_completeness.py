#!/usr/bin/env python3
"""Exhaustive raw -> summary accounting for one sample's funcscan + MAP outputs (and optionally the anvi'o import).
Every raw hit must be in the summary, or be removed by a documented rule; nothing may be in a summary without a
raw source. Complements scripts/check_map_outputs.sh (which looks for empty / dummy outputs).

  check_annotation_completeness.py --sample S --map-dir <MAP outdir>/S --funcscan-dir <funcscan outdir> \
      --genes-gff <anvi'o gene export .gff3> [--functions <build_anvio_functions outdir>/functions_all.txt]

Rules encoded (from the pipelines' code, v5.0.0 / funcscan 4.0.0; update here if versions change):
  funcscan  hAMRonization keeps every row of abricate / amrfinderplus / rgi / deeparg (+ potential) / fargene,
            except AMRFinderPlus hits with Method INTERNAL_STOP (dropped; also from argNorm) -> WARN listing them;
            rows filtered by input_file_name, as batch runs share one summary;
            argNorm normalises abricate / amrfinderplus / deeparg row by row;
            AMPcombi keeps ampir predictions with prob >= amp_ampcombi_parsetables_cutoff and length <= ..._aalength
            (read from pipeline_info/params_*.json);
            dbCAN overview = every gene hit by dbCAN_hmm, dbCAN-sub or DIAMOND, #ofTools = number of those.
  MAP       mobilome: geNomad plasmid / virus (provirus -> prophage) with score > 0.8; ISEScan complete ('c') only;
            IntegronFinder complete integrons only; compositional outliers (2 BED lines each); minus
            <sample>_discarded_mge.txt (mge<500bp, no_cds, CO_overlap_with_MGE, RNAs_in_window).
            combined report rows = PathoFact2 gff proteins + AMR proteins (seeds); amr_tool = tools that hit it;
            bgc_type for seeds in the BGC gff; mge_type = mobilome features covering >= 90% of the CDS.
Prints OK / WARN / FAIL per check with counts and examples; exit 1 if any FAIL.
"""
import argparse, csv, glob, gzip, json, re, sys
from collections import Counter, defaultdict

ap = argparse.ArgumentParser()
ap.add_argument("--sample", required=True); ap.add_argument("--map-dir", required=True)
ap.add_argument("--funcscan-dir", required=True); ap.add_argument("--genes-gff", required=True)
ap.add_argument("--functions", help="functions_all.txt from build_anvio_functions.py (optional)")
a = ap.parse_args()
S, M, F = a.sample, a.map_dir, a.funcscan_dir
csv.field_size_limit(sys.maxsize)
res = []
def report(sec, name, ok, detail, level="FAIL"):
    st = "OK" if ok else level
    res.append(st); print(f"{st:4} [{sec}] {name}: {detail}")
def cmp_sets(sec, name, raw, summ, level="FAIL", note=""):
    miss, extra = raw - summ, summ - raw
    ok = not miss and not extra
    d = f"raw {len(raw)}, summary {len(summ)}"
    if miss: d += f"; MISSING from summary {len(miss)} e.g. {sorted(miss)[:3]}"
    if extra: d += f"; EXTRA in summary {len(extra)} e.g. {sorted(extra)[:3]}"
    report(sec, name, ok, d + (f" ({note})" if note and not ok else ""), level)
def opn(p): return gzip.open(p, "rt") if p.endswith(".gz") else open(p)
def tsv(p, comment=None):
    with opn(p) as fh:
        lines = [l for l in fh if l.strip() and not (comment and l.startswith(comment))]
    return list(csv.DictReader(lines, delimiter="\t"))
def gid(x): return x.split()[0]
def ids_in(gff, ftype="CDS"):
    out = set()
    for l in opn(gff):
        f = l.rstrip("\n").split("\t")
        if len(f) >= 9 and f[2] == ftype:
            m = re.search(r"ID=([^;]+)", f[8]); out.add(m.group(1))
    return out

# ======================= funcscan =======================
FA = f"{F}/arg"
ham = tsv(f"{F}/reports/hamronization_summarize/hamronization_combined_report.tsv")
hcount = Counter(r["analysis_software_name"] for r in ham)
hrows = defaultdict(Counter)
ham = [r for r in ham if r["input_file_name"].split(".")[0] == S]  # batch runs share one summary
hcount = Counter(r["analysis_software_name"] for r in ham)
for r in ham: hrows[r["analysis_software_name"]][r["input_sequence_id"]] += 1
def rows(p, comment=None): return tsv(p, comment) if glob.glob(p) else []
raw = {
    "abricate": rows(f"{FA}/abricate/{S}/{S}.txt"),
    "amrfinderplus": rows(f"{FA}/amrfinderplus/{S}/{S}.tsv"),
    "rgi": rows(f"{FA}/rgi/{S}/{S}.txt"),
    "deeparg": rows(f"{FA}/deeparg/{S}/{S}.mapping.ARG") + rows(f"{FA}/deeparg/{S}/{S}.mapping.potential.ARG"),
}
# hAMRonization (and argNorm, which reads its output) do not output AMRFinderPlus hits with Method INTERNAL_STOP
# (likely pseudogene / sequencing error). Reported as WARN with the genes, excluded from the row comparisons.
istop = [r for r in raw["amrfinderplus"] if r.get("Method") == "INTERNAL_STOP"]
raw["amrfinderplus"] = [r for r in raw["amrfinderplus"] if r.get("Method") != "INTERNAL_STOP"]
report("funcscan", "AMRFinderPlus INTERNAL_STOP hits (not in hAMRonization / argNorm)", not istop,
       f"{len(istop)}" + (": " + ", ".join(f"{r['Element symbol']} on {r['Contig id']}" for r in istop) if istop else ""), "WARN")
idcol = {"abricate": ["SEQUENCE"], "amrfinderplus": ["Contig id", "Protein id", "Protein identifier"], "rgi": ["Contig", "ORF_ID"], "deeparg": ["read_id"]}
for t, rr in raw.items():
    sc = Counter({k: v for k, v in hrows.get(t, Counter()).items()})
    rc = Counter()
    for c in idcol[t]:  # hAMRonization's input_sequence_id = whichever raw column matches
        if rr and c in rr[0]:
            rc = Counter(gid(r[c]) for r in rr)
            if rc == sc: break
    ok = rc == sc
    report("funcscan", f"hAMRonization <- {t}", ok,
           f"raw {sum(rc.values())} rows / summary {sum(sc.values())} rows" +
           ("" if ok else f"; ids differing: {sorted(set(rc) ^ set(sc))[:3]}"))
fg = 0
for f in glob.glob(f"{FA}/fargene/{S}/*/results_summary.txt"):
    for l in open(f):
        m = re.search(r"predicted genes.*?:\s*(\d+)", l, re.I)
        if m and "full" in l.lower(): fg += int(m.group(1))
report("funcscan", "hAMRonization <- fargene", hcount.get("fargene", 0) == fg,
       f"fARGene full-length predicted genes {fg}, summary rows {hcount.get('fargene', 0)}")
for t, p in (("abricate", f"{FA}/argnorm/abricate/{S}.normalized.tsv"),
             ("amrfinderplus", f"{FA}/argnorm/amrfinderplus/{S}.normalized.tsv"),
             ("deeparg", f"{FA}/argnorm/deeparg/{S}.ARG.normalized.tsv")):
    n_norm = len(rows(p)); n_raw = len(rows(f"{FA}/deeparg/{S}/{S}.mapping.ARG")) if t == "deeparg" else len(raw[t])
    report("funcscan", f"argNorm <- {t}", n_norm == n_raw, f"raw {n_raw} rows / argNorm {n_norm} rows")
par = {}
for p in glob.glob(f"{F}/pipeline_info/params_*.json"): par.update(json.load(open(p)))
cut, maxlen = float(par.get("amp_ampcombi_parsetables_cutoff", 0.6)), int(par.get("amp_ampcombi_parsetables_aalength", 120))
amp_raw = {r["seq_name"] for r in tsv(f"{F}/amp/ampir/{S}/{S}.ampir.tsv") if float(r["prob_AMP"]) >= cut and len(r["seq_aa"]) <= maxlen}
amp_sum = {r["CDS_id"] for r in tsv(f"{F}/reports/ampcombi2/{S}/{S}_ampcombi.tsv")}
cmp_sets("funcscan", f"AMPcombi <- ampir (prob >= {cut}, length <= {maxlen})", amp_raw, amp_sum)
D = f"{F}/cazyme/dbcan/cazyme_annotation/{S}/{S}"
hmm = {r["Target Name"] for r in tsv(f"{D}_dbCAN_hmm_results.tsv")}
sub = {r["Target Name"] for r in tsv(f"{D}_dbCANsub_hmm_results.tsv")}
dia = {r["Gene ID"] for r in tsv(f"{D}_diamond.out")}
ov = tsv(f"{D}_overview.tsv")
cmp_sets("funcscan", "dbCAN overview <- hmm | sub | diamond", hmm | sub | dia, {r["Gene ID"] for r in ov})
bad = [r["Gene ID"] for r in ov if int(r["#ofTools"]) != (r["Gene ID"] in hmm) + (r["Gene ID"] in sub) + (r["Gene ID"] in dia)]
report("funcscan", "dbCAN overview #ofTools", not bad, f"{len(ov) - len(bad)}/{len(ov)} consistent" + (f"; e.g. {bad[:3]}" if bad else ""))

# ======================= MAP =======================
P = f"{M}/prediction"
rep = tsv(f"{M}/{S}_combined_report.tsv")
rp = {r["protein_id"]: r for r in rep}
amr = {
    "amrfinderplus": {gid(r["Protein id"]) for r in rows(f"{P}/amr_genes/amrfinderplus/{S}.tsv")},
    "rgi": {gid(r["ORF_ID"]) for r in rows(f"{P}/amr_genes/rgi/{S}.txt")},
    "deeparg": {gid(r["read_id"]) for r in rows(f"{P}/amr_genes/deeparg/{S}.mapping.ARG")},
}
for t, s in amr.items():
    cmp_sets("MAP", f"report amr_tool <- {t}", s, {p for p, r in rp.items() if t in r["amr_tool"].split(",")})
amr_all = set().union(*amr.values())
cmp_sets("MAP", "integrated AMR gff <- all AMR tools", amr_all, ids_in(f"{P}/amr_genes/integrated_{S}.gff"))
pf = ids_in(f"{P}/virulence/{S}_pathofact2.gff")
cmp_sets("MAP", "report rows <- PathoFact2 gff + AMR (seeds)", pf | amr_all, set(rp))
pfa = {}
for l in open(f"{P}/virulence/{S}_pathofact2.gff"):
    f = l.rstrip("\n").split("\t")
    if len(f) >= 9 and f[2] == "CDS":
        at = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv); pfa[at["ID"]] = at
mism = [p for p, at in pfa.items() if p in rp and
        (("pathofact2_vf_prob" in at) != (rp[p]["pathofact2_vf_prob"] != "-") or ("pathofact2_tox_prob" in at) != (rp[p]["pathofact2_tox_prob"] != "-"))]
report("MAP", "report vf/tox probabilities <- PathoFact2 gff", not mism, f"{len(pfa) - len(mism)}/{len(pfa)} consistent" + (f"; e.g. {mism[:3]}" if mism else ""))
bg = ids_in(f"{P}/bgcs/{S}_bgcs.gff")
cmp_sets("MAP", "report bgc_type <- BGC gff (seed proteins only)", bg & set(rp), {p for p, r in rp.items() if r["bgc_type"] != "-"})
reg_b = len(ids_in(f"{P}/bgcs/{S}_bgcs.gff", "bgc_region"))
report("MAP", "BGC regions (SanntiS only in this setup)", reg_b > 0 or not bg, f"{reg_b} regions, {len(bg)} CDS in BGC gff", "WARN")

# mobilome accounting
mob = defaultdict(int); feats = defaultdict(list)
for l in opn(f"{M}/gff/{S}_mobilome.gff.gz"):
    f = l.rstrip("\n").split("\t")
    if len(f) < 9 or l.startswith("#"): continue
    mob[f[2]] += 1
    if f[2] not in ("CDS", "direct_repeat_element", "inverted_repeat_element"):
        lab = re.search(r"mobile_element_type=([^;]+)", f[8]).group(1)
        feats[f[0]].append((int(f[3]), int(f[4]), lab))
disc = Counter()
for l in open(f"{M}/{S}_discarded_mge.txt"):
    f = l.rstrip("\n").split("\t")
    if len(f) >= 3:
        k = "co" if f[0] == "CO" else re.sub(r"_\d+$", "", f[0])  # 'CO <id> <reason>' or '<mge_id> ... <reason>'
        disc[k] += 1
G = f"{P}/genomad"
pl = [r for r in rows(f"{G}/{S}_5kb_contigs_plasmid_summary.tsv") if float(r["plasmid_score"]) > 0.8]
vi = [r for r in rows(f"{G}/{S}_5kb_contigs_virus_summary.tsv") if float(r["virus_score"]) > 0.8]
nprov = sum("provirus" in r["seq_name"] for r in vi)
exp = {"plasmid": len(pl) - disc["plas"], "viral_sequence+prophage": len(vi) - disc["vir1"]}
report("MAP", "mobilome plasmid <- geNomad score > 0.8 - discarded", mob["plasmid"] == exp["plasmid"],
       f"geNomad {len(rows(f'{G}/{S}_5kb_contigs_plasmid_summary.tsv'))}, > 0.8 {len(pl)}, discarded {disc['plas']}, mobilome {mob['plasmid']}")
report("MAP", "mobilome viral + prophage <- geNomad score > 0.8 - discarded", mob["viral_sequence"] + mob["prophage"] == exp["viral_sequence+prophage"],
       f"geNomad {len(rows(f'{G}/{S}_5kb_contigs_virus_summary.tsv'))}, > 0.8 {len(vi)} ({nprov} provirus), discarded {disc['vir1']}, mobilome viral {mob['viral_sequence']} + prophage {mob['prophage']}")
iss = rows(f"{P}/isescan/{S}_1kb_contigs.fasta.tsv")
isc = sum(r["type"] == "c" for r in iss)
report("MAP", "mobilome insertion_sequence <- ISEScan complete - discarded", mob["insertion_sequence"] == isc - disc["iss"],
       f"ISEScan {len(iss)}, complete {isc}, discarded {disc['iss']}, mobilome {mob['insertion_sequence']}")
nco = sum(1 for b in glob.glob(f"{P}/compositional_outliers/*.bed") for l in open(b) if l.strip() and not l.startswith("#")) // 2
report("MAP", "mobilome compositional_outlier <- BED - discarded", mob["compositional_outlier"] == nco - disc["co"],
       f"outliers {nco}, discarded {disc['co']}, mobilome {mob['compositional_outlier']}")
isum = glob.glob(f"{P}/integronfinder/*.summary"); comp = 0
for f in isum:
    for l in open(f):
        x = l.rstrip("\n").split("\t")
        if len(x) > 2 and x[2].isdigit(): comp += int(x[2])
dummy = any("contig_dummy" in f for f in isum)
report("MAP", "mobilome integron <- IntegronFinder complete", mob.get("integron", 0) + disc.get("int", 0) == comp and not dummy,
       f"complete integrons {comp}, mobilome {mob.get('integron', 0)}" + ("; DUMMY IntegronFinder output (patch 5 / check_map_outputs.sh)" if dummy else ""),
       "WARN")
others = {k: v for k, v in mob.items() if k not in ("CDS", "direct_repeat_element", "inverted_repeat_element", "plasmid", "viral_sequence",
                                                    "prophage", "insertion_sequence", "compositional_outlier", "integron")}
report("MAP", "mobilome other feature types", not others, f"{others or 'none'}", "WARN")

# report mge_type <- recomputed CDS overlap (>= 90 %)
cds = {}
for l in open(a.genes_gff):
    f = l.rstrip("\n").split("\t")
    if len(f) >= 9 and f[2] == "CDS":
        cds[re.search(r"ID=([^;]+)", f[8]).group(1)] = (f[0], int(f[3]), int(f[4]))
def norm(lab): return "IS" if re.match(r"IS|ISL|ISNCY|Tn", lab) and lab not in ("plasmid",) else lab
bad = []
for p, r in rp.items():
    c, s, e = cds[p]
    want = {norm(lab) for fs, fe, lab in feats.get(c, []) if (min(e, fe) - max(s, fs) + 1) / (e - s + 1) >= 0.9}
    got = {norm(x) for x in r["mge_type"].split(",") if x != "-"}
    if want != got: bad.append((p, sorted(want), sorted(got)))
report("MAP", "report mge_type <- mobilome features >= 90% of CDS", not bad,
       f"{len(rp) - len(bad)}/{len(rp)} consistent" + (f"; e.g. {bad[:3]}" if bad else ""))
sp_ips = set()
for l in gzip.open(f"{P}/interproscan/{S}.tsv.gz", "rt"):
    f = l.split("\t")
    if len(f) > 3 and f[3] in ("SignalP_GRAM_POSITIVE", "SignalP_GRAM_NEGATIVE"): sp_ips.add(f[0])
spf = glob.glob(f"{M}/{S}_combined_report.signalp.tsv")
if spf:
    cmp_sets("MAP", "SignalP-filled report <- InterProScan (bacterial models), report proteins", sp_ips & set(rp),
             {r["protein_id"] for r in tsv(spf[0]) if r["signalP"] != "-"})
else:
    report("MAP", "SignalP-filled report", False, "missing -> scripts/add_signalp_to_report.py", "WARN")

# ======================= anvi'o import (optional) =======================
if a.functions:
    fx = defaultdict(set)
    for l in open(a.functions):
        f = l.split("\t", 3)
        if f[0] != "gene_callers_id": fx[f[1]].add(f[0])
    G2 = lambda s: {x.rsplit("___", 1)[1] for x in s}
    cmp_sets("anvio", "dbCAN <- overview", G2({r["Gene ID"] for r in ov}), fx["dbCAN"])
    cmp_sets("anvio", "AMR_all <- MAP AMR tools", G2(amr_all), fx["AMR_all"])
    cmp_sets("anvio", "VFDB <- report vfdb_hit", G2({p for p, r in rp.items() if r["vfdb_hit"] != "-"}), fx["VFDB"])
    cmp_sets("anvio", "PathoFact2_toxin <- report tox prob >= 0.9", G2({p for p, r in rp.items() if r["pathofact2_tox_prob"] != "-" and float(r["pathofact2_tox_prob"]) >= 0.9}), fx["PathoFact2_toxin"])
    cmp_sets("anvio", "AMP_all <- AMPcombi", G2(amp_sum), fx["AMP_all"])
    cmp_sets("anvio", "BGC_SanntiS <- BGC gff", G2(bg), fx["BGC_SanntiS"])
    cmp_sets("anvio", "SignalP <- InterProScan (bacterial)", G2(sp_ips), fx["SignalP"])

n = Counter(res); print(f"== {S}: {n['OK']} OK, {n['WARN']} WARN, {n['FAIL']} FAIL")
sys.exit(1 if n["FAIL"] else 0)
