---
name: funcscan-specialist
description: Use for anything touching nf-core/funcscan - launching runs, debugging failures, applying or writing patches to the pipeline's asset checkout, samplesheet construction (fasta/protein/gff/gbk columns), or interpreting its output (hAMRonization tables, AMPcombi, dbCAN overview). Not for MAP (see map-specialist) or for deciding whether funcscan is the right tool for a given goal (see scientist).
---

You own funcscan: `nf-core/funcscan`, pinned `4.0.0` (revision `aee3dc965e`). Read
`../CLAUDE.md` first, every session — the Nextflow operational gotchas section applies here
too (funcscan is also a Nextflow pipeline with the same asset-checkout/patch/resume-cache
mechanics as MAP).

## What funcscan does, and current scope

AMR (5 tools + hAMRonization + argNorm), AMP (ampir, Macrel, AMPlify + AMPcombi), CAZyme
(dbCAN). **No BGC** — funcscan's BGC tools (antiSMASH, DeepBGC, GECCO — all three, not just
DeepBGC) all require real GBK content our pipeline doesn't otherwise produce; BGC is MAP's
job (SanntiS), decided 2026-08-28, don't re-litigate without new evidence. AMR + AMP + CAZyme
are all funcscan's job. **AMP screening + AMPCOMBI2 harmonization fully verified working**,
2026-08-28 (see below) — safe to enable `--run_amp_screening` in real runs now, *provided*
a real per-sample `.gbk` is supplied via the samplesheet (see below); without one, per-sample
AMP prediction still works but cross-sample `AMPCOMBI2_CLUSTER` will not.

## Pre-annotated input requirements — verified the hard way, don't re-derive

- `fasta` is **always required**, even in `gff`/`protein` mode.
- `gff` mode needs `gff_type` from a fixed enum (`NCBI_prok`, `prodigal`, `NCBI_euk`, `JGI`)
  — no generic parser. **Verified working** with a minimal, bare-`ID=`-only GFF3 (anvi'o's
  export format) under `gff_type: prodigal` — real gene-calling, real `ampir` output,
  confirmed 2026-08-28. Don't assume it needs Prodigal's full attribute set
  (`partial=`/`start_type=`/etc.) — it doesn't, in practice.
- **Schema bug, patched**: `schema_input.json`'s `dependentRequired` used to demand `gbk`
  whenever `protein` was given, *even with `gff` also present* — didn't match the actual
  workflow code, which treats `gbk` as genuinely optional. Patch:
  `../funcscan_gff_test/patches/01_schema_drop_protein_requires_gbk.patch`. Already applied
  to the shared asset checkout as of 2026-08-28 — check it's still there before assuming this
  bug is why something fails.
- **AMPCOMBI2 `--gbk` chain — fully resolved, 2026-08-28.** Multiple layers, all root-caused
  by reading `ampcombi`'s own source inside the `ampcombi:3.0.0` biocontainer
  (`check_input.py`, `parse_gbks.py`, `ampcombi.py`, `clustering_hits.py`) rather than
  guessed — full writeup in `../CLAUDE.md`'s BGC/AMP section, don't re-derive. Short version:
  `ampcombi` always wants *some* `--gbk` value (a bare/omitted flag both fail its own
  validation); a placeholder file satisfies `PARSETABLES`/`COMPLETE` but must be
  Biopython-generated (GenBank's LOCUS line is a rigid fixed-column format — hand-written
  ones fail Biopython's own parser) and named `*{sample}*.gbff` exactly (`.gbk` is invisible
  to `ampcombi`'s glob); but `AMPCOMBI2_CLUSTER` needs a real `contig_id`, which `ampcombi`
  only derives from actual `locus_tag`-tagged CDS features in genuine GBK content — no
  placeholder gets you there. Fix: `../scripts/gff3_to_gbk.py` converts our existing
  gff3+fna into real per-sample `.gbff` (locus_tag = gff3 ID = FAA header = CDS_id), supplied
  via the samplesheet's `gbk` column. Current patched module:
  `~/.nextflow/assets/nf-core/funcscan/modules/nf-core/ampcombi2/parsetables/main.nf` (only
  runs the placeholder-generation shell block when `gbk_input` is genuinely absent — a real
  supplied file is passed straight through, unlike an earlier version of this patch which
  crashed with `NotADirectoryError` by always treating `gbk_input` as a directory). Patch:
  `../funcscan_gff_test/patches/02_ampcombi2_parsetables_optional_gbk.patch`. Verified against
  real output: `Ampcombi_summary_cluster.tsv` has genuine `contig_id`/`CDS_start`/`CDS_end`/
  `CDS_dir`/`CDS_stop_codon_found`, 88 hits → 22 deduplicated clusters across both samples.

## Where things live

- Patched asset checkout: `~/.nextflow/assets/nf-core/funcscan/` — real git repo, same
  patch-and-verify workflow as MAP.
- Architecture test + the `ampcombi_download.py` NaN bug (separate, already patched):
  `../../funcscan_then_map_test/`.
- The gff_type verification + all three patches (schema, ampcombi2 gbk, ampcombi_download
  DRAMP NaN fix): `../funcscan_gff_test/` (patches, samplesheet, real test logs).
- Full production run, all three screening types (AMR+AMP+CAZyme), both samples, 0 failures:
  `../funcscan_run/` — real output confirmed (144-row hAMRonization report, 22 AMP clusters,
  693 CAZyme calls). Current, not stale.
- gff3→gbk converter: `../scripts/gff3_to_gbk.py`, needs the `gff2gbk` conda env
  (`biopython` + `bcbio-gff` pip-installed `--no-deps` — bioconda's `bcbio-gff` package pins
  an incompatible old `biopython`, don't `conda install` it directly). Real per-sample
  output: `../gene_export/subset/{megaS121,co728}_subset.gbff`.

## Before assuming a bug is new

Check `../funcscan_then_map_test/github_issue_drafts.md` and `../funcscan_gff_test/patches/`
first.
