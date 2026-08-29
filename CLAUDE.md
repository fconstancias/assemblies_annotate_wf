# CLAUDE.md — Metagenome functional annotation & longitudinal quantification

Read this first, in every session, regardless of which agent you are. It's the distilled
operational knowledge from the work that preceded this directory (`../assembly_to_MGE/`,
`../funcscan_then_map_test/`, and `../assembly_to_MGE/native_signalp_test/`) — not a
transcript, just what you need to not re-derive things the hard way a second time.

## The big picture

For every metagenome assembly (single-sample and coassembly), determine per-gene: is it AMR?
AMP? a BGC? a CAZyme? a virulence factor? does it sit on a plasmid/MGE-classified contig?
which MAG does it belong to? — then quantify it: bowtie2 mapping information, plus
potentially the anvi'o databases already built for these assemblies, to get contig-level and
gene-level coverage (breadth, depth), tracked across the longitudinal sample series per
participant.
KEGG/metabolic annotation is deferred ("maybe later"), not blocking anything here.

**Strategy**: combine existing, maintained community pipelines (MAP + nf-core/funcscan)
rather than hand-build tool wiring — no single one covers everything needed, and both are
actively maintained upstream rather than something we'd have to keep patching forever
ourselves (we patch what we must, minimally, and document every patch — see below).

## Upstream: the assembly/binning pipeline

Repo: `fconstancias/shotgun_metagenomics_binning` (`dev` branch), local checkout at
`../binning_smk/`. Produces contigs (single-sample + coassembly), bins/MAGs (concoct,
metabat2, vamb, binette, dRep-dereplicated), GTDB-Tk taxonomy — and, as of 2026-08-26, gene
calls too: GFF3 + protein FASTA per sample/group, exported from each group's existing anvi'o
contigs-db via `anvi-get-sequences-for-gene-calls` (see
`../spa_single_all/results_spa_single_all/export_gff3_faa.sh`). This means **neither MAP nor
funcscan needs to do its own gene-calling anymore** — a real architecture change from the
earlier funcscan-first/MAP-consumes design (that existed only to share ONE tool's
gene-calling with the other; now a third, upstream source does it for both).

**Known quirk, fix verified end-to-end 2026-08-28**: the GFF3 and FAA are two *separate*
`anvi-get-sequences-for-gene-calls` calls (`--export-gff3` and `--get-aa-sequences` are
mutually exclusive in one invocation; `--export-gff3` also rejects `--defline-format`
outright — "not compatible with the GFF3 output mode", confirmed by trying). By default
anvi'o's FASTA defline uses only `{gene_caller_id}` while `--export-gff3`'s `ID=` attribute
uses `{contigs_db_project_name}___{gene_caller_id}` — **not `{contig_name}`, a wrong guess
corrected after actually checking**: on `megaS121.db`, the GFF's `ID=megaS121___0` uses the
anvi'o project name ("megaS121"), while the contig itself is named
`megaS121_000000000007` — visibly different strings, easy to guess wrong from the pattern
alone. **Neither downstream pipeline's schema validates content** (both only check file
extensions), so this would silently produce wrong/empty results if not fixed, not an error.

**Verified real fix**: run `--get-aa-sequences` with `--defline-format
"{contigs_db_project_name}___{gene_caller_id}"` — confirmed on `megaS121.db` directly:
GFF's `ID=megaS121___0/1/2` now matches FAA's `>megaS121___0/1/2` exactly. Also confirmed
clean: GFF's CDS feature count (94,863) exactly equals FAA's sequence count (94,863) — the
"50 gene calls have empty AA sequences and skipped" warning anvi'o prints is `rRNA` features
(non-coding, correctly excluded from FAA, still present in GFF3 as non-`CDS` features
downstream tools should already ignore) — not a gap.

Two ways to apply this: (a) patch `export_gff3_faa.sh` in `../binning_smk` (`dev` branch) and
re-export there, or (b) run `anvi-get-sequences-for-gene-calls` directly ourselves against
the already-existing `BINNING/01_ANVIO_DBs/{megaS121,co728}/{sample}.db` (same underlying
assemblies as everything else in this doc; no binning-pipeline rebuild needed, no dev-branch
change, lower risk) — done this way for our own test pair, see
`gene_export/megaS121_test.{gff3,faa}` (verified) and `gene_export/co728.{gff3,faa}` (real
coassembly, larger, in progress as of this writing).

## Downstream: the two annotation pipelines

- **MAP** (`EBI-Metagenomics/mobilome-annotation-pipeline`, pinned `v5.0.0`): MGE/plasmid/
  virus (geNomad, ICEfinder2, IntegronFinder, ISEScan), virulence/toxin (PathoFact2), AMR
  (AMRFinderPlus, RGI, DeepARG), BGC (SanntiS, GECCO, antiSMASH).
- **funcscan** (`nf-core/funcscan`, pinned `4.0.0`): AMR (5 tools + hAMRonization + argNorm),
  AMP (ampir, Macrel, AMPlify + AMPcombi), CAZyme (dbCAN). **No BGC** — see BGC-ownership
  finding below, `--run_bgc_screening` stays off; MAP/SanntiS owns BGC entirely.
- **BGC ownership — DECIDED, 2026-08-28: MAP/SanntiS owns BGC screening. funcscan's BGC
  tools (DeepBGC, antiSMASH, GECCO — all three, not just DeepBGC) are not usable at all with
  our current pre-annotated gff3+faa export**, confirmed by directly reading
  `subworkflows/local/bgc.nf`: every one of the three takes `gbks` (a real, content-parsed
  GBK file) as its input channel, unconditionally — there is no gff+faa code path for BGC in
  funcscan. Confirmed by running it for real (`bgc_comparison/`, `--run_bgc_screening
  --bgc_skip_antismash --bgc_skip_gecco`, DeepBGC only): `DEEPBGC_DOWNLOAD` ran, but
  `DEEPBGC_PIPELINE` itself was scheduled **zero times** for either sample and the pipeline
  finished with `WARN: No hits found by BGC tools` for both — not a crash, just silently
  nothing to do, because our samplesheet's `gbk` resolves to `[]` (traced in
  `workflows/funcscan.nf` line ~93: `gbk_found != null ? gbk_found : []`). And unlike the
  AMPCOMBI2 `--gbk` bug (below), this isn't a placeholder-file trick away from working —
  DeepBGC/antiSMASH/GECCO need *real* CDS content in the GBK to make calls at all, and anvi'o
  has no GBK export capability (confirmed directly: only `anvi-script-process-genbank*`
  *import* tools exist, nothing exports). Building a real GBK from our gff3+fna+faa would be
  a genuine engineering task, not a quick fix — not worth it when SanntiS already works
  today with zero extra plumbing and produces real calls (9 BGC/sample on the curated
  subset, `map_run/results/*/prediction/bgcs/sanntis/*_sanntis.gff.gz` — real
  `nearest_MiBIG`/`nearest_MiBIG_class`/`score` fields) and gets IPS reuse for free (see
  "SignalP / export prediction" below). **funcscan's `--run_bgc_screening` should stay off
  for this project.**
- **Run independently, not sequentially.** Now that gene-calling comes from the binning
  pipeline, neither tool depends on the other's output — no efficiency reason left to
  sequence them (the old funcscan-first ordering is obsolete).
- Both need the underlying assembly FASTA too, not just gff+faa: MAP's `assembly` column,
  funcscan's `fasta` column, both required regardless of pre-annotated input.
- **funcscan `gff_type: prodigal` — verified working, 2026-08-28.** Real test: ran `ampir`
  AMP screening against our anvi'o-exported minimal GFF3 (bare `ID=` only, none of real
  Prodigal's `partial=`/`start_type=`/etc.) + matching FAA, curated subset, real assembly
  fasta. Ran clean through gene-calling/GFF-parsing and produced real, correct output
  (real gene IDs, real sequences, real AMP probability scores) — the minimal attribute set
  is not a problem in practice.
- **Second schema bug found and patched along the way**: funcscan's `schema_input.json` had
  `"protein": ["gbk"]` in `dependentRequired` — i.e. supplying `protein` always demanded a
  `gbk` file too, *even when `gff` was also given*. But the actual workflow code
  (`workflows/funcscan.nf`) treats `gbk` as genuinely optional (`gbk_found != null ?
  gbk_found : []`; `preannotated: gff != [] || gbk != []` — either alone qualifies) — a real
  schema/implementation mismatch, not intentional. Patched: dropped `"protein": ["gbk"]`
  from the schema (kept `"gbk": ["protein"]` and `"gff": ["protein"]`, both still correct).
  See `funcscan_gff_test/patches/01_schema_drop_protein_requires_gbk.patch`.
- **Third bug, AMPCOMBI2/`ampcombi` `--gbk` chain — RESOLVED, 2026-08-28, real end-to-end
  verified output.** With the schema fix applied and `gff`+`protein` supplied, the pipeline
  ran correctly through `ampir` itself, then hit a sequence of `ampcombi`-side issues, each
  root-caused by reading `ampcombi`'s own source (`check_input.py`, `parse_gbks.py`,
  `ampcombi.py`, `clustering_hits.py` inside the `ampcombi:3.0.0` biocontainer) rather than
  guessed:
  1. `AMPCOMBI2_PARSETABLES` unconditionally passes `--gbk ${gbk}`; empty/absent renders as
     a bare `--gbk` with no argument, which `ampcombi`'s argparse rejects.
  2. Omitting `--gbk` entirely doesn't work either — `ampcombi`'s own Python validation
     (`check_gbk_path()`) still demands *some* value, file or directory.
  3. A placeholder directory alone isn't enough — `check_gbk_path()` globs strictly for
     `*{sample}*.gbff` (case-sensitive; a `.gbk`-named file is invisible to it, despite an
     inline comment implying broader matching).
  4. Once correctly named, the file's *content* is genuinely parsed too
     (`Bio.SeqIO.parse(..., "genbank")`), and GenBank's LOCUS line is a rigid, exact-column
     format a hand-written placeholder failed. Fixed by generating the placeholder with
     Biopython itself (guarantees round-trip validity) instead of hand-matching the spec.
  5. With that placeholder, `AMPCOMBI2_PARSETABLES`/`_COMPLETE` ran clean, but
     `AMPCOMBI2_CLUSTER` then crashed (`KeyError: 'contig_id'`) — `contig_id` is only
     populated by matching real `locus_tag` values found in actual GBK CDS features
     (`ampcombi.py:355`, `parse_gbks.py:212-216`); an empty placeholder has none. This is the
     same root architectural fact as the BGC finding above: some `ampcombi`/funcscan stages
     need *real* annotated GBK content, not just a syntactically-valid file — but unlike BGC
     tools (which need real functional/product content we can't produce), AMPCOMBI2 only
     needs real coordinates + locus_tag-per-CDS, and we already have exactly that in our own
     gff3. Built `scripts/gff3_to_gbk.py` (uses `BCBio.GFF` + Biopython — installed in a
     dedicated `gff2gbk` conda env; note bioconda's `bcbio-gff` package has a stale
     dependency pin incompatible with its own code's Biopython API requirement, worked around
     by pip-installing `bcbio-gff --no-deps` over a modern conda `biopython`) to convert our
     existing gff3+fna into real, multi-record `.gbff` files — one CDS feature per gene, each
     `locus_tag` set to the GFF3 `ID` (== FAA header == `ampcombi`'s `CDS_id`), so the
     locus_tag match actually succeeds. Verified: round-trips cleanly through the real
     `ampcombi` container's Biopython (1.80) before ever touching Nextflow; CDS counts match
     the known curated-subset gene counts exactly (megaS121: 100 contigs/5,075 CDS; co728: 99
     contigs/6,022 CDS).
  6. Found and fixed one more, self-inflicted bug while wiring the real gbk in: the
     `parsetables/main.nf` patch's placeholder-generation logic unconditionally treated
     `gbk_input` as a directory to `mkdir -p`/`ls` into. A *real* supplied gbk is a plain
     *file* (which `check_gbk_path()` already accepts directly) — crashed with
     `NotADirectoryError`. Fixed: the placeholder-generation shell block now only runs when
     `gbk_input` is genuinely absent; a real supplied file/dir is passed straight through.
  Generated real `.gbff` for both curated-subset samples
  (`gene_export/subset/{megaS121,co728}_subset.gbff`), added as the `gbk` column in
  `funcscan_gff_test/samplesheet.csv`, full rerun succeeded end to end (`PARSETABLES` →
  `COMPLETE` → `CLUSTER` → `MULTIQC`, 0 failures). Real output confirmed: 38 AMP hits
  (megaS121) + 50 (co728) → 88-row merged summary → 22 deduplicated cross-sample clusters,
  with genuine `contig_id`/`CDS_start`/`CDS_end`/`CDS_dir`/`CDS_stop_codon_found` columns
  populated from the real GBK content (a real value-add beyond just unblocking the crash).
  Patches: `funcscan_gff_test/patches/02_ampcombi2_parsetables_optional_gbk.patch` (final
  version), `03_ampcombi_download_dramp_nan_fix.patch` (unrelated DRAMP-DB-download NaN
  crash, patched along the way — pandas reads an empty `Sequence` cell as `NaN`, not `""`;
  `re.match()` on a float crashes `download_ref_db()`).
  **Next**: rerun `funcscan_run/` (production-scope curated-subset test) with
  `--run_amp_screening` now included — it was excluded from that run pending this fix.

## SignalP / export prediction: use the MAP/IPS route, not the custom build

A native SignalP6-in-PathoFact2 module was built and validated (`../assembly_to_MGE/
native_signalp_test/` — full patch set, working, documented) but turned out **not** to be
meaningfully faster than just running real InterProScan (which SanntiS needs as a
prerequisite anyway) once real candidate-set sizes were measured — the "small VF/TOX
candidate subset" assumption was wrong in practice (a classifier-calibration issue, not an
extraction-script one; see that directory's README for the full investigation). Current
plan (per 2026-08-28 decision): **drop the custom SignalP6 module**, run MAP with
`--skip_sanntis false` (SanntiS needs IPS; `interpro_licensed_software=true` already gets you
real SignalP4.1 inside that IPS run for free) and skip antiSMASH/GECCO (`--skip_gecco true
--skip_antismash true`, redundant/expensive, funcscan/MAP split covers BGC without them both).

**Issue 10/4 still applies and isn't optional**: MAP never routes its own internally-computed
IPS output into `COMBINEREPORTER` — the `signalP` column stays `-` unless you *also* supply
`interproscan_tsv` in the samplesheet pointing back at the IPS output the same run just
produced. Two-pass: run once without it, then add it and `-resume` once IPS has completed.
Side effect, accepted: this also flips `PATHOFACT2_INTEGRATOR`'s annotation source from `cdd`
to `ips` (not a bug, just a real, documented tradeoff — see
`../assembly_to_MGE/README.md`'s troubleshooting section).

**`--goterms --pathways` for InterProScan**: feasible, wanted. Needs `--iprlookup` added too
(GO/pathway mappings hang off the InterPro cross-reference lookup, not the raw member-DB
hits). Requires a small patch to MAP's `conf/modules.config` `INTERPROSCAN` block (currently
`--applications CDD,TIGRFAM,GENE3D,PRINTS,PROSITEPATTERNS,PFAM,SIGNALP`). Expected cheap
(annotation lookup on hits already found, not a new search) but not yet measured — verify
before promising a number.

## Testing discipline — the pattern that's already paid off twice

**Build a curated, categorized subset from already-known signal — never a naive/random
subset, never wait on a full-scale run to validate a pipeline change.** Concretely: pick
~100 contigs per sample spanning plasmid-flagged (geNomad), virus-flagged (geNomad),
virulence-positive (strongest VFDB hits), AMR-positive, random-long (crosses MAP's 100kb
RENAME threshold — otherwise IntegronFinder/compositional-outlier-detection never get
exercised), and random-tiny (crosses the 1kb/5kb thresholds) — all grounded in a prior real
run's actual output (`combined_report.tsv`, geNomad's own summaries), never guessed. See
`../assembly_to_MGE/native_signalp_test/tiny_assembly_subset/build_subset.py` for the working
implementation. This caught real things a naive subset wouldn't have (e.g. would never have
noticed the GFF/FAA ID mismatch above without deliberately building a representative test).

Payoff, measured: full pipeline validation in ~25-30 minutes instead of the ~1.5 days the
real full-scale assemblies need for the expensive steps (was InterProScan/SignalP; will be
whatever the current expensive step is for new work).

**Scope rule for this directory: every test runs on exactly one coassembly + one
single-sample assembly, subsampled** — the same pair used throughout all prior work
(`megaS121` single + `co728` coassembly, participant 728), for direct comparability. Not
full-size, not more than one of each — one pair already exercises both assembly types
(contig-count scale differs ~4x between them). Move to full-size / more participants only
once a change is validated on this pair.

## Nextflow operational gotchas (all confirmed the hard way — don't rediscover)

- **`nextflow`/`conda`/`singularity` aren't on PATH in a fresh non-interactive shell.**
  `~/.bashrc` has an explicit `case $- in *i*) ;; *) return ;; esac` guard partway through —
  everything below it (conda init included) only runs for *interactive* shells, which a
  fresh non-interactive Bash invocation isn't. Every session needs, explicitly, every time:
  `source /opt/software/mamba/23.3.1/etc/profile.d/conda.sh && conda activate env_nf` (the
  dedicated nextflow env — found via `find ~/.conda/envs -iname nextflow`, not obvious from
  the env name list) plus `source /usr/share/Modules/init/bash && module load
  singularity/3.8.7`. No shell state carries between separate tool calls, so this has to be
  in the same command chain as the actual `nextflow run`.
- **`env_nf`'s conda-pinned nextflow (25.10.2) is older than funcscan 4.0.0 needs
  (>=25.10.4)** — fails fast with a clear version-mismatch message, not a cryptic one, but
  easy to lose time on. Fix is one-time and persists in the env: `nextflow self-update`
  (it's a self-updating launcher, not a normal conda package — jumped straight to 26.04.6).
- **`$TMPDIR` leaking into containers breaks tools that actually write real temp files.**
  This session's shell has `TMPDIR=/scratch/tmp/...` set (matches the scratchpad
  convention), which `nextflow`/singularity pass straight through to the container process —
  but singularity only auto-binds `/tmp`, `$HOME`, and cwd by default, *not* arbitrary custom
  paths, so a container trying to actually create a directory there fails outright
  (`AMRFinderPlus`: "Error creating a temporary directory in /scratch/tmp/" — cost a full,
  otherwise-successful production run). Fix: `export TMPDIR=/tmp` in the same command chain
  as the `nextflow run`, every time — same one-shell-per-call caveat as the PATH gotcha
  above.
- **`bin/`-staged script edits don't invalidate cache.** Nextflow doesn't hash the content of
  a module's `resources/usr/bin/*` or nf-core `bin/`-staged helper scripts (they're on PATH,
  not declared inputs). Editing one and expecting `-resume` to notice: it won't, unless the
  *calling process's* own declared input/signature also changed (structurally forces a new
  task hash). Otherwise, manually delete the specific cached task work dirs (found via
  `grep "PROCESS_NAME" .nextflow.log | grep -o "workDir: [^ ]*"`) before resuming.
- **Cross-directory resume needs the cache index copied, not just a shared `workDir`.** The
  task-hash → output index (`db` + `index.*` files) lives in the *launch directory's own*
  `.nextflow/cache/<session-uuid>/`, not inside `workDir`. A fresh directory pointed at an
  old `workDir` still says "never run this project before." Fix: `cp -r
  <old_dir>/.nextflow/. <new_dir>/.nextflow/` (contents, not the directory itself — nesting
  gotcha), then `-resume <exact-session-uuid>` (from the old dir's `.nextflow/history`), not
  a run name (names are local-history lookups only).
- **Remote-form invocation only** (`nextflow run ORG/repo -r TAG`), never a local clone path —
  switching invalidates the *entire* resume cache, not just the changed part. Patch the real
  `~/.nextflow/assets/ORG/repo/` checkout directly (it's a git repo — `git diff`/`git stash`/
  `git apply --check` all work normally for building and verifying patches).
- **`unset TMPDIR` before every nextflow invocation** — an inherited TMPDIR pointing somewhere
  that doesn't exist inside Singularity containers breaks container startup silently/oddly.
- **A killed Nextflow controller doesn't cancel its already-submitted Slurm jobs.** They keep
  running independently and their output is still usable on the next `-resume` — don't
  reflexively `scancel` them, check what they actually are first (`scontrol show job`).
- **Licensed/no-public-container tools** (SignalP4.1 via IPS's bind-mount, SignalP6 for the
  now-deprecated custom module): install at task runtime into an isolated `pip install
  --target=` dir inside a plain public base image, rather than trying to bake a custom
  image (no `singularity build --fakeroot` available on this cluster, no authenticated
  remote builder). `PYTHONNOUSERSITE=1` is required or pip silently falls back to `$HOME/
  .local` if anything there already satisfies a dependency. `./` prefix required on any
  locally-staged package directory passed to `pip install` (a bare Nextflow-staged basename
  reads as a PyPI package name, not a path).
- **`executor.queueSize` copied from a tiny test config silently throttles a real production
  run** — 50 (fine for a 2-sample curated-subset test) capped both funcscan and MAP to 50
  concurrent submitted tasks each during the 19-participant production run, while the
  cluster sat with hundreds of CPUs idle (confirmed via `squeue`/`sinfo`). This is an
  executor-level submission-throttle setting, not a per-task directive — safe to bump and
  `-resume` freely, doesn't invalidate the cache (contrast with the next point). Bumping
  50→200 took actual concurrent jobs from 125 (330 CPUs idle) to 265 (partition at ~full
  utilization). Check real idle capacity (`sinfo -p <queue> -o "%D %C"`) before assuming a
  copied-over queueSize is still appropriate at a new scale.
- **Editing a `process{}`-level directive (e.g. adding/changing a `memory` override)
  invalidates the resume cache for *every* process, not just the one touched** — cost real,
  already-completed work (a full fARGene re-run) via "Unable to resume cached task" across
  the board after one `memory` edit. Executor-level settings (previous point) don't have
  this problem; only `process{}`-level directives (`memory`, `cpus`, `errorStrategy`, etc.)
  feed into the cache hash. Know which kind of edit you're making before touching a config
  mid-run.
- **A `memory = { task.memory * task.attempt }` retry-scaling closure recurses infinitely**
  and crashes the run (`java.lang.StackOverflowError`) — Nextflow resolves `task.memory` by
  *calling this same closure*, so referencing it inside its own definition is circular, not
  "whatever the process already had." Use a fixed literal base instead:
  `{ [16.GB * task.attempt, 128.GB].min() }`, never `task.memory`, in a directive that's
  scaling that same directive.
- **Some tools' own wrappers swallow a SIGKILL (OOM) and re-raise as a plain exit 1**, hiding
  the real signal from Nextflow's exit-status-based retry logic — MAP's `GENOMAD` (wraps
  `mmseqs2` via its own Python CLI) did this; the real cause only showed up in
  `.command.err` ("died with `<Signals.SIGKILL: 9>`"), while Nextflow only saw exit 1 and
  (with a retry policy scoped to 130-145/104/255 only) sent it straight to `finish` with no
  retry — silently dropping the affected samples' results for that tool. Broaden
  `errorStrategy` to also retry on exit 1 if this class of run needs to be resilient; accepts
  wasting some retries on genuinely-broken exit-1 failures as the cost of not silently
  losing samples.
- **Heavy concurrent Singularity use can exhaust the system's loop devices**: `failed to
  find loop device... resource temporarily unavailable`. Transient infra contention, not a
  real task failure — covered by the same retry-with-backoff policy as the exit-1 case
  above (`errorStrategy`/`maxRetries`, both configs in `coassembly_production/`).
- **Real production-scale data (hundreds of thousands of genes/sample) can exceed a
  pipeline's default resource labels even when the exact same pipeline+flags worked fine on
  a curated subset** — hit a genuine OOM (`AMRFinderPlus`'s `tblastn`, exit 137) only at
  19-participant/230K-696K-gene scale, never during any subset-scale validation. Don't
  assume subset-validated resource defaults carry over to real scale without a retry/scaling
  safety net (previous points).
- **`anvi-run-kegg-kofams`/`anvi-estimate-metabolism` mutate the contigs-db they're given
  in place** (write KOfam/KEGG_Module/KEGG_BRITE functions directly into it) — pointing a
  production job at a shared upstream contigs-db (as opposed to our own copy) breaks the
  "never mutate upstream source data" convention followed everywhere else in this project.
  Caught before any actual damage in the 19-participant production run (confirmed via both
  log content and file mtimes that none had written yet) but cost ~3h45m of already-elapsed
  HMM search, redone against copies under `coassembly_production/anvio_kegg/contigs_dbs/`
  instead. Check *what a tool actually writes to, and where*, before pointing any anvi'o
  annotation step at a path outside your own workflow directory — not just funcscan/MAP
  which only ever read gff3/faa/fasta inputs and write to their own `results/`.

## Agent team for this directory

Defined and written 2026-08-28: `.claude/agents/*.md` (hidden dir, won't show in an IDE file
tree unless it's configured to show dotfiles — the files are real, check with `ls -la`).
Each one reads `CLAUDE.md` (this file) itself as the shared source of truth and is scoped to
not duplicate it — the summaries below are just an index.

- **`scientist`** — keeps the overall project arc in view: the upstream binning pipeline, the
  real research question, and whether a technical decision actually serves it. Consult
  before design decisions that aren't pure pipeline mechanics — not for routine execution of
  already-agreed work. Encodes the actual stated priority order: **plasmid, AMR, virulence
  (and specifically export-status via the MAP/SanntiS/IPS route, not a separate SignalP
  build) come first; BGC is nice-to-have, not yet decided whether funcscan or MAP should own
  it; redundant AMR between MAP and funcscan is explicitly fine, since the point is
  funcscan's harmonized hAMRonization table, not novel coverage.** Also holds the "always
  test on the subsampled pair first" discipline as a standing check, alongside
  `pipeline-validator` owning the *mechanics* of that.
- **`map-specialist`** — owns MAP's patch set (asset checkout is a real git repo), DB config
  (`../assembly_to_MGE/my_paths.config`), and its specific bugs/workarounds (10 documented,
  `github_issue_drafts.md`). Knows the `run_native_signalp` undefined-param gotcha from the
  still-applied (but now unused) native-SignalP6 patches.
- **`funcscan-specialist`** — same, for funcscan. Knows the exact state of all 3 real bugs
  found in this directory's own testing: the schema `dependentRequired` bug (patched), the
  `AMPCOMBI2_PARSETABLES` `--gbk` bug (patched one layer), and the deeper `ampcombi` 3.0.0
  bug underneath that (not yet fixed, parked — see "Downstream" section above for the full
  writeup).
- **`pipeline-validator`** — owns the curated-subset methodology (`build_subset.py` pattern)
  and the one-coassembly-one-single-sample scope rule. The one to reach for before scaling
  any change up to full-size assemblies, or when something needs testing fast rather than
  waiting on a multi-hour/day run.

## Where prior work lives (reference, not duplicated here)

- `../assembly_to_MGE/` — original MAP test, DB setup scripts, 10 documented upstream bugs
  (`github_issue_drafts.md`), full patch history.
- `../funcscan_then_map_test/` — funcscan architecture test, full ARG/AMP/BGC/CAZyme run
  completed for both samples, its own bug (`ampcombi_download.py` NaN crash) + patch.
- `../assembly_to_MGE/native_signalp_test/` — the native SignalP6 build (superseded, see
  above, but the patches/README are the reference for the "small subset assumption was
  wrong" investigation and the tiny-assembly testing methodology).

## Coverage / quantification — the actual end goal, not yet built

Everything above is annotation. This is what it's *for*: once genes are annotated
(AMR/AMP/BGC/CAZyme/virulence/MGE), quantify each one via coverage, tracked across the
longitudinal sample series per participant. Leaning anvi'o-native, since contigs-dbs already
exist for these assemblies:
- [`anvi-profile-blitz`](https://anvio.org/help/9/programs/anvi-profile-blitz/) — fast
  contig-level coverage/detection, no full profile-db needed.
- [`anvi-export-gene-coverage-and-detection`](https://anvio.org/help/8/programs/anvi-export-gene-coverage-and-detection/)
  — per-gene, if a profile-db exists.
- [`anvi-script-get-coverage-from-bam`](https://anvio.org/help/main/programs/anvi-script-get-coverage-from-bam/)
  — BAM-only, no contigs-db/profile-db at all.
- [`anvi-export-splits-and-coverages`](https://anvio.org/help/main/programs/anvi-export-splits-and-coverages/)
  — per-split.

**Open, not yet scoped**: how to bring binette's bin assignments into an anvi'o
**collection**, so coverage can be reported per-bin/MAG, not just per-contig/per-gene. Check
`../assembly_to_MGE/CLAUDE.md`'s older "Coverage / quantification" section first (CoverM vs.
`anvi-script-get-coverage-from-bam` comparison, an output-layout sketch already exists there)
before re-deriving from scratch.

## Production run: all 19 real co-assemblies — started 2026-08-29, under a hard deadline

**Context that overrides normal pace**: cluster access ends 2026-08-31 (Monday) — everything
in this section is being done under real time pressure, prioritizing throughput and
maximum background parallelism over the careful one-sample-at-a-time verification style used
everywhere else in this doc. See `coassembly_production/README.md` for full detail; summary
here.

**Scope correction, important**: earlier in this same day, real assembler conventions were
initially guessed wrong (assumed megahit was the project's uniform default, based only on
which two samples happened to get used for testing all session — `megaS121`/`co728`). The
actual, established project design (`../spa_single_all/SINGLE_ASSEMBLY.md`,
`../spa_coassembly_all/CO_ASSEMBLY.md`) is: **single assemblies use SPAdes** (277 assembly
groups across all 25 participants, `spa_single_all/`), **co-assemblies use MEGAHIT**
(per-participant, `spa_coassembly_all/`, reusing `B01_megahit_co_assembly/sel_all/
comega{participant}/`). `co728`/`co894`/`cospa728`/`cospa894` under `BINNING/01_ANVIO_DBs/`
are explicitly an older/different run, superseded — not used going forward. Don't re-guess
this; read the two README.md files above if it comes up again.

**What's running now** (co-assemblies only so far — single assemblies, 277 groups, not
started, likely out of scope given the deadline unless explicitly revisited):
19 of 25 participants (6 excluded per `CO_ASSEMBLY.md`: pending more samples). Real contigs-
dbs already exist, gene calls already done —
`../spa_coassembly_all/results_coassembly_all/concoct/mh_p{ID}/mh_p{ID}.db` (named
`concoct/` for historical reasons in that pipeline, confirmed via its `anvi_gen_contigs_db`
Snakemake rule to be the genuine full anvi'o contigs-db, not a CONCOCT-internal artifact).
Real scale: 230,757-695,576 genes/participant, ≈6.6M total — 15-55x the curated subset used
for every prior test in this project. Running: gene export (gff3+faa, same defline-format
fix), `anvi-run-kegg-kofams`+`anvi-estimate-metabolism` (SLURM array job, NOT direct-on-node
like the subset tests — genuinely heavy compute at this scale), funcscan (AMR+AMP+CAZyme)
and MAP (full scope) as 19-row multi-sample Nextflow runs. AMPCOMBI2 uses the placeholder-gbk
patch (real per-sample AMP hits, no cross-sample clustering) — building real per-participant
`.gbk` was deprioritized given the deadline. Some of the largest participants may not finish
every tool before access ends; that's accepted, not a failure.

## Open items, as of 2026-08-28

1. ✅ **Done — GFF3/FAA ID mismatch, real fix verified.** Turned out simpler than planned:
   both `megaS121.db` and `co728.db` (same assemblies as everything else in this doc)
   already existed in `BINNING/01_ANVIO_DBs/`, un-exported — no need to touch `../binning_smk`
   or rebuild anything. Ran `anvi-get-sequences-for-gene-calls` directly against them with
   `--defline-format "{contigs_db_project_name}___{gene_caller_id}"` on the FAA side; IDs
   match exactly, CDS-count = FAA-sequence-count exactly for both samples (megaS121: 94,863;
   co728: 333,487). Filtered both down to the existing curated 100-contig subset selection
   (`selections.pkl`, reused as-is). Files: `gene_export/{megaS121,co728}.{gff3,faa}` (full),
   `gene_export/subset/{megaS121,co728}_subset.{gff3,faa}` (curated subset, 5,075 / 6,022
   genes), plus `subset/{megaS121,co728}_subset.gbff` (real GenBank, converted from the
   subset gff3+fna via `scripts/gff3_to_gbk.py` — see the AMPCOMBI2 writeup above; same CDS
   counts, `locus_tag` = gff3 `ID` = FAA header). **Fixed upstream too, 2026-08-28**: turned out `export_gff3_faa.sh` was never a
   checked-in pipeline script — the real source of truth is `../binning_smk/README.md`'s
   documented export commands (added 2026-08-26, `33bd354`); that ad-hoc script was just a
   one-off run directly from them. Added the same `--defline-format` fix there, with a
   "required, not optional" note explaining why (same silent-failure risk — neither MAP nor
   funcscan's schema validates ID correspondence, only file existence/format). Committed and
   pushed to `dev`: `7426a0f`. Every future export from this pipeline's own documented
   commands now gets correct IDs automatically.
2. Coassembly GFF3/FAA — not a separate task anymore, folded into item 1 (co728's `.db`
   already had it done the same way).
3. ✅ **Done — funcscan `gff_type: prodigal` verified working**, plus schema
   `dependentRequired` bug, patched. See "Downstream" section above for the full writeup.
4. ✅ **Done — AMPCOMBI2 `--gbk` chain, fully resolved, real gff3→gbk converter built.** See
   "Downstream" section's full writeup and `scripts/gff3_to_gbk.py`. Real end-to-end output
   verified (88 AMP hits → 22 deduplicated clusters, both samples, genuine `contig_id`/
   coordinate columns). AMP screening is now safe to enable in real runs, provided the
   samplesheet's `gbk` column points at a real per-sample `.gbff`.
5. ✅ **Done — full MAP + funcscan run, both pipelines, both samples, full scope including
   AMP, against the real curated subset.** funcscan (`funcscan_run/`): AMR + CAZyme + AMP, all
   three screening types together, 35 completed + 31 cached, 0 failed. Real output confirmed
   across all three: 144-row `hamronization_combined_report.tsv` (AMR, 5 tools harmonized),
   22-cluster `Ampcombi_summary_cluster.tsv` (AMP, same real result as the isolated test),
   693 combined CAZyme calls (`{sample}_overview.tsv`, dbCAN). Hit one real infra gotcha along
   the way (fixed, see gotchas section): this session's `$TMPDIR=/scratch/tmp/...` leaking
   into the container broke `AMRFinderPlus` outright ("Error creating a temporary directory")
   — `export TMPDIR=/tmp` before the `nextflow run` fixed it. MAP (`map_run/`): full scope
   (`--skip_sanntis false --skip_gecco true --skip_antismash true`, real InterProScan +
   SanntiS + PathoFact2 + AMR), 44m53s, 92/92 tasks succeeded — **this also answers whether
   MAP's own GFF parser accepts the anvi'o-exported format**: yes, no errors, real output
   (real gene IDs matching the fixed `{sample}___{gene_caller_id}` scheme, real VFDB hits,
   real PathoFact2 probabilities, real MGE-type calls). One config gotcha hit: the
   now-unused native-SignalP6 patches from `../assembly_to_MGE/native_signalp_test/` are
   still applied to the shared MAP asset checkout, so `params.run_native_signalp` must be
   explicitly set (`false`) or the run errors on an undefined param — see `map_run/
   nextflow.config`. **This is the first full, real, working confirmation of the entire new
   architecture end to end — both pipelines, every planned screening type, on the real
   curated subset, zero failures.**
6. Decide on `--goterms --pathways` for InterProScan (feasible, not yet built).
7. **Research, not yet scoped**: does [DRAM-IT](https://dramit.readthedocs.io/en/latest/params_doc.html)
   replace or complement anything in the MAP/funcscan split? Unknown until actually checked
   — don't assume either way going in.
8. Bringing binette's bins into an anvi'o collection for per-MAG coverage — see the Coverage
   section above, this is the main open piece of the actual end goal.
9. ✅ **Done — BGC ownership decided.** MAP/SanntiS owns BGC entirely; funcscan's
   `--run_bgc_screening` stays off (all three of its BGC tools need real GBK content our
   pipeline doesn't otherwise produce — a genuinely different, heavier requirement than
   AMPCOMBI2's, which item 4's converter already satisfies). See "Downstream" section above.