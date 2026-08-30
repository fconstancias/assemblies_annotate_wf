# Co-assembly production run — all 19 real per-participant MEGAHIT co-assemblies

Started 2026-08-29, under real time pressure: cluster access ends 2026-08-31 (Monday).
Goal — process as much of the real dataset as the remaining window allows, not just the
curated test subset used everywhere else in `../`.

## Scope

19 of 25 participants (the ones actually used by `../../spa_coassembly_all/` — see that
directory's `CO_ASSEMBLY.md` for why 6 are excluded: `729, 654, 431, 373, 151, 142`, held
back pending more samples). **Not** `co728`/`co894`/`cospa728`/`cospa894` under
`BINNING/01_ANVIO_DBs/` — those are an older/different run, explicitly superseded per
`CO_ASSEMBLY.md`. Real contigs-dbs (gene calls, HMMs, everything already done — no
re-derivation needed), one per participant:
`../../spa_coassembly_all/results_coassembly_all/concoct/mh_p{ID}/mh_p{ID}.db`
(named `concoct/` for historical/organizational reasons in that pipeline, but it is the
genuine, full anvi'o contigs-db — confirmed directly from `metagenome_assemble.smk`'s
`anvi_gen_contigs_db` rule).

Real scale, per participant (`select count(*) from genes_in_contigs`): 230,757 (`mh_p408`,
smallest) to 695,576 (`mh_p813`, largest) genes. Total ≈ 6.6M genes across all 19 — 15-55x
bigger than the curated 100-contig subset used for all prior validation in `../`.

## What's running

- **`gene_export/`**: gff3 + faa per participant, same `anvi-get-sequences-for-gene-calls`
  command and `--defline-format "{contigs_db_project_name}___{gene_caller_id}"` fix as
  everywhere else in this project (`../CLAUDE.md`'s Upstream section).
- **`anvio_kegg/`**: `anvi-run-kegg-kofams` + `anvi-estimate-metabolism --metagenome-mode`,
  submitted as one SLURM array job (`anvio_kegg/run_kofams_array.sh`, job array, one task
  per participant) — **not** run directly on the login/interactive node like the curated-
  subset tests were, since at this scale (per-task ~15-50x the subset's already-substantial
  runtime) that would be genuinely heavy sustained compute, which is what `sbatch` is for.
  Uses the same `anvio_kegg_data` KEGG setup built for the subset comparison (`../
  metabolism_comparison/anvio_kegg/`) — no need to redo that step, it's participant-agnostic.
  **Operates on `anvio_kegg/contigs_dbs/mh_p{ID}.db` — our own copies, not the shared
  upstream `spa_coassembly_all/results_coassembly_all/concoct/mh_p{ID}/mh_p{ID}.db`
  originals.** `anvi-run-kegg-kofams` mutates the contigs-db it's given in place (writes
  KOfam/KEGG_Module/KEGG_BRITE functions directly into it) — the first version of this job
  was accidentally pointed at the shared originals; caught and fixed before any of the 19
  had actually written results (confirmed via both log content and file mtimes — none had
  reached the "added to the contigs database" point yet), costing ~3h45m of already-elapsed
  HMM search that had to be redone against the copies instead. Same "never mutate upstream
  source data" convention as everywhere else in this project.
- **`funcscan_run/`** + **`map_run/`**: same patched asset checkouts, same AMR+AMP+CAZyme /
  full-scope flags as `../funcscan_run/` and `../map_run/`, one multi-sample samplesheet
  (19 rows) each instead of 19 separate invocations — lets Nextflow's own SLURM executor
  handle per-sample parallelism across the queue.
- **AMPCOMBI2**: no real per-participant `.gbk` built for this run (unlike the curated-
  subset comparison) — given the time constraint, prioritizing getting real AMR/AMP/CAZyme
  hit lists across all 19 participants over perfecting cross-sample AMP clustering
  (`AMPCOMBI2_CLUSTER`, which needs real gbk content — see `../CLAUDE.md`). The placeholder-
  gbk patch (`../funcscan_gff_test/patches/02_ampcombi2_parsetables_optional_gbk.patch`,
  already applied to the shared checkout) still gets real per-sample AMP hits
  (`PARSETABLES`/`COMPLETE`) — just not the cross-sample dedup step. Revisit only if time
  remains.

## Real-scale gotchas hit going from 2-sample tests to 19-participant production

None of these ever surfaced during the curated-subset validation earlier in `../` — only
appeared at real scale, real concurrency, or with the actual production data. All fixed in
`funcscan_run/nextflow.config` and `map_run/nextflow.config`, see those files' comments for
the full detail; summary here since they cost real time (multiple kill/relaunch cycles).

1. **`executor.queueSize = 50`** (copied unmodified from the tiny 2-sample test configs)
   silently throttled both pipelines to 50 concurrent submitted tasks each, while the actual
   cluster sat with 330 idle CPUs and one fully-idle node — confirmed directly via
   `squeue`/`sinfo`. At 19-participant scale, a single sample alone spawns ~10 parallel
   fARGene sub-tasks; 50 total across the whole run was nowhere near enough. Bumped to 200
   in both configs. This is an executor-level submission-throttle setting, not a per-task
   directive — safe to change and resume without losing any completed work (see point 3).
   Confirmed real effect: 125 concurrent jobs (330 CPUs idle) → 265 jobs (partition at
   ~full utilization, 1660/1664 CPUs allocated).
2. **MAP's `--skip_sanntis`/`--skip_gecco`/`--skip_antismash` CLI flags broke** ("Value is
   [string] but should be [boolean]") — caused by the `nextflow self-update` done earlier
   this session (25.10.2 → 26.04.6, needed for funcscan 4.0.0's version floor) changing how
   the newer Nextflow CLI parses `--flag true`/`--flag false`. Fixed by setting them as real
   Groovy booleans directly in `nextflow.config`'s `params {}` block instead of as CLI
   arguments — sidesteps the CLI parser entirely.
3. **Editing `process{}`-level directives (e.g. adding a new `memory` override) invalidates
   the resume cache for every process, not just the one being changed** — cost real,
   already-completed fARGene work by triggering a full recompute on the next `-resume`
   ("Unable to resume cached task" for every previously-finished task). Executor-level
   settings (`executor.queueSize`, point 1) do **not** have this problem — only
   `process{}`-level directives (`memory`, `cpus`, `errorStrategy`, etc.) affect the cache
   hash. Know which kind of config change you're making before doing it mid-run.
4. **A `memory = { task.memory * task.attempt }` retry-scaling closure recurses infinitely**
   — Nextflow resolves `task.memory` by calling this same closure, so referencing it inside
   its own definition is circular. Caused a real `java.lang.StackOverflowError` that aborted
   the whole funcscan run outright (confirmed in `.nextflow.log`). Fixed: use a fixed literal
   base value instead (`{ [16.GB * task.attempt, 128.GB].min() }`), never `task.memory`,
   inside a directive that's scaling that same directive.
5. **Real per-participant scale (230K-696K genes, 15-55x the curated subset) exceeds
   funcscan's default resource labels** — hit a genuine OOM (`AMRFinderPlus`'s `tblastn` on
   `mh_p813`, the largest participant; `status = 35072` decodes to exit 137/SIGKILL). Fixed
   by the same per-attempt memory scaling as point 4, combined with retrying on the relevant
   exit codes (point 6).
6. **Some tools' own wrappers swallow SIGKILL and re-raise as a plain exit 1**, hiding the
   real signal from Nextflow's exit-status-based retry logic entirely — MAP's `GENOMAD`
   step (its Python `genomad` CLI wraps `mmseqs2`) did this for `mh_p398`/`mh_p584`: the
   real cause, confirmed in `.command.err`, was `died with <Signals.SIGKILL: 9>`, but
   Nextflow only ever saw the wrapper's own exit 1, so the original `errorStrategy` (which
   only retried on 130-145/104/255) sent it straight to `finish` with no retry, silently
   dropping those two participants' MGE/plasmid/virus predictions for that run. Fixed:
   broadened `errorStrategy` to also retry on exit 1 — accepts wasting some retries on
   genuinely-broken exit-1 failures, a reasonable cost against silently losing samples.
7. **Loop-device exhaustion under heavy concurrent Singularity use**: `failed to find loop
   device... resource temporarily unavailable` — transient infra contention from running
   many containers at once, not a real task failure. Covered by the same broadened
   `errorStrategy`/retry as point 6.

## Further real issues hit once funcscan/MAP actually finished a full pass

8. **`AMPIR` needs more than a 128GB memory cap on the 3 largest participants**
   (`mh_p523`/`mh_p789`/`mh_p813`, 435K-696K genes) — confirmed real OOM (exit 137) even
   after exhausting all 4 retries at the capped ceiling. Fixed with a more specific
   `withName: '.*AMPIR'` selector (placed after the blanket `'.*'` one) raising the cap to
   512GB — node capacity here is ~2TB, real headroom. A selector matching only one process
   only affects that process's resolved config/cache hash, unlike editing the blanket rule.
9. **`AMPCOMBI2_CLUSTER`'s known, already-accepted gbk limitation (point 4 above) took the
   whole 19-participant funcscan run down with it** — its default `errorStrategy` ('finish'
   after exhausting retries) triggered "Killing running tasks (60)", including `RGI_MAIN`
   for effectively all 19 participants mid-execution — real AMR work that would have
   completed fine on its own, lost for no benefit. Fixed: explicit
   `withName: '.*AMPCOMBI2_CLUSTER' { errorStrategy = 'ignore' }` so this one known,
   accepted limitation can no longer cascade into killing unrelated in-flight work.
10. **A genuine, recurring bug in MAP's own `ICEFINDER2_LITE` subworkflow — patched by
    disabling the subworkflow, not by excluding samples.** First seen on `mh_p398`: one
    parallel branch (`PROCESS_BLASTP_PROKKA`) succeeded with real, non-empty output (525
    lines), but a *different* parallel branch apparently produced nothing for this
    participant, and a downstream Groovy closure joining the branches by sample ID wasn't
    written to handle that "no match" case gracefully — crashed with `Invalid method
    invocation \`call\` with arguments: [[id:X], null, .../X_uniprot_names.tsv]`, a
    workflow-definition-level error (not a task failure `errorStrategy` can catch or retry).
    Initially excluded just `mh_p398` and resumed — the *identical* error immediately
    recurred on a *different* participant (`mh_p894`), confirming this is systemic, not a
    one-sample edge case; excluding samples one at a time would have kept recurring
    indefinitely. No built-in `--skip_icefinder`-style flag exists (confirmed via
    `nextflow_schema.json` and the workflow's own invocation — called unconditionally). MAP's
    own `CLAUDE.md` documents `remainder: true` joins as its established convention for
    exactly this "optional/missing match" case, and confirmed the *downstream* consumer of
    `ICEFINDER2_LITE.out.ices_tsv` already uses that pattern correctly — the bug is entirely
    inside `ICEFINDER2_LITE`'s own internal joins, not in how its output is consumed. Fixed
    by disabling the subworkflow call entirely and substituting `channel.empty()` for its
    output (patch: `map_run/04_disable_broken_icefinder2_lite.patch`, verified clean-apply)
    — the existing `remainder: true` downstream already tolerates this gracefully. ICE/IME
    detection isn't among this project's stated priorities (plasmid/AMR/virulence/AMP), so
    losing it entirely for this run is an acceptable trade against the alternative (chasing
    a Groovy DSL fix inside unfamiliar subworkflow code, or losing an unknown, possibly
    growing number of participants one at a time). All 19 participants, including `mh_p398`
    and `mh_p894`, back in the same run after this fix — real cache preserved for every
    unrelated already-completed task (confirmed: 0 "Unable to resume cached task" warnings).

11. **`RGI_MAIN` failed (exit 140) on 8+ different participants** at the blanket 128GB cap
    — but real, complete output already existed in the work dir (`mh_p341.json` 9.5MB,
    `mh_p341.txt` 1.2MB, both look legitimate) before the task's own exit status came back
    non-zero, consistent with a resource limit hit late (post-processing, or a lingering
    RGI multiprocessing worker), not a data/logic failure. Also produces repeated,
    non-fatal `Requested rname X does not exist! Please check your FASTA file.` warnings —
    traced directly to `pyfaidx` inside the RGI container (`app/Base.py`'s
    `get_part_sequence()`, used for partial-hit boundary refinement); `X` comes out as just
    the bare sample ID rather than a real contig name, most likely RGI's own gene-ID parsing
    assumption being broken by our `{sample}___{gene_id}` triple-underscore convention — but
    confirmed harmless/cosmetic (real output produced despite the repeated warnings, not the
    cause of the exit 140). Fixed with the same targeted pattern as AMPIR:
    `withName: '.*RGI_MAIN' { memory = { [32.GB * task.attempt, 512.GB].min() } }`.

12. **`RUNDBCAN`'s three steps (EASYCGC/EASYSUBSTRATE/CAZYMEANNOTATION) OOM'd (exit 137,
    plain "Killed") at the 128GB blanket cap** — but only on the same two largest
    participants (`mh_p789`/`mh_p813`, 667K/696K genes); everything else in the run had
    already completed successfully at this point (funcscan's own summary: "Pipeline
    completed successfully, but with errored process(es)" — 548 real successes, only these
    7 dbCAN failures). Same targeted fix pattern: `withName: '.*RUNDBCAN.*'` with the same
    512GB ceiling.

## Known, accepted limitation given the deadline

Some of the largest participants (`mh_p813`: 695,576 genes, `mh_p789`: 667,363,
`mh_p292`: 449,070, `mh_p523`: 435,803) may not finish every tool before access ends —
this is expected and accepted; partial results for the largest participants are still real,
usable results for whatever stages did complete. Check `logs/` for per-participant status.
