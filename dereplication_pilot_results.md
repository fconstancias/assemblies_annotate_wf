# Dereplication pilot — results (2 participants, mmseqs2)

Companion to `dereplication_brainstorm.md` (the design) and
`dereplication_pilot/` (the code/output). Ran 2026-09-06, per
`/home/ljc444/.claude/plans/happy-moseying-clarke.md`.

## Setup

For each participant, pooled that participant's own co-assembly + all of their own
single-sample assemblies (contigs and, separately, ORFs), then ran
`mmseqs easy-cluster`:

- **Tier 1 (contigs, nucleotide)**: `≥1kb` length filter, then
  `--min-seq-id 0.95 -c 0.8 --cov-mode 2 --cluster-mode 2` (containment-aware:
  query-coverage mode lets a short fragment cluster into a longer representative
  as long as most of *the fragment itself* aligns).
- **Tier 2 (ORFs, protein)**: `--min-seq-id 0.95 -c 0.9 --cov-mode 1 --cluster-mode 2`.

| Participant | co-assembly | single-sample groups |
|---|---|---|
| 110 | `mh_p110` | `spaS30, spaS37, spaS222, spaS223` (4) |
| 97 | `mh_p97` | `spaS25, spaS29, spaS35, spaS43, spaS250, spaS251, spaS252, spaS253` (8) |

**Real bug hit and fixed**: `run_pilot.sh` originally stored the group-name list
in a shell variable named `GROUPS` — that collides with bash's own special
`GROUPS` variable (the process's Unix group-ID list), so the intended string got
silently clobbered and the first real run came back with 100% of cluster members
unattributable to any source group. Fixed by renaming to `PARTICIPANT_GROUPS`;
re-ran just the (cheap) analysis step against the already-computed (expensive)
mmseqs clustering output — no need to re-cluster.

## Results

### Tier 1 — contigs

| Participant | total contigs (≥1kb) | redundant (count) | redundant (bp) |
|---|---|---|---|
| 110 | 174,051 | 43.9% | 71.4% |
| 97 | 191,613 | 42.0% | 73.3% |

Per-group, the co-assembly itself is consistently the *most* redundant slice
(56.3% / 54.5% of its own contigs land in a cluster shared with ≥1 single-sample
assembly) — expected, since the co-assembly is literally built by pooling those
same single-sample reads. Redundancy is bp-weighted much higher than count-weighted
(~71-73% vs ~42-44%) — longer, well-covered contigs are far more likely to reassemble
consistently across samples than short/fragmentary ones.

### Tier 2 — ORFs

| Participant | total genes | redundant (count) | redundant (bp) |
|---|---|---|---|
| 110 | 1,022,879 | 93.2% | 96.0% |
| 97 | 1,286,297 | 96.5% | 98.0% |

Far more extreme than Tier 1 — the overwhelming majority of predicted genes across
a participant's own co-assembly + single-sample assemblies are near-identical
copies of each other, not just similar. Every single-sample group individually
shows 92-99% redundancy against the rest of that participant's own data.

## Parameter sensitivity sweep (2026-09-06)

Reran clustering against the already-pooled fasta (no re-pooling needed) at
alternate settings, per `dereplication_pilot/run_sweep.sh`, to check whether the
above numbers depend heavily on the specific thresholds chosen:

**Tier 1 — `--cov-mode 0` (bidirectional, mmseqs's own default) vs. the original
`--cov-mode 2` (query-coverage/containment-aware)**:

| Participant | cov-mode 2 (original) | cov-mode 0 (bidirectional) |
|---|---|---|
| 110 | 43.9% count / 71.4% bp | 20.8% count / 34.4% bp |
| 97 | 42.0% count / 73.3% bp | 25.1% count / 47.8% bp |

Roughly halves. Expected, not a red flag: bidirectional coverage requires *both*
sequences to be well-covered, so a short single-sample contig fully contained
within a much longer co-assembly contig of the same region fails that test even
though it's a genuine match. Since "the same region assembles to different
lengths at different depths" is exactly the real scenario this pilot targets,
`--cov-mode 2` is the deliberately-correct choice, not an inflated one — but the
Tier 1 number is conditional on treating containment as redundancy; cov-mode 0's
result is the more conservative "near-full-length reciprocal match only"
baseline.

**Tier 2 — `--min-seq-id 0.99` (stricter) vs. the original `0.95`**:

| Participant | 0.95 (original) | 0.99 (stricter) |
|---|---|---|
| 110 | 93.2% count / 96.0% bp | 90.6% count / 94.0% bp |
| 97 | 96.5% count / 98.0% bp | 94.4% count / 96.4% bp |

Only a 2-3 point drop under a much stricter identity cutoff — **the gene-level
redundancy signal is robust**, not an artifact of a loose threshold. Genes really
are near-identical (>99%) copies across a participant's own assemblies, not just
family-level similar.

**Conclusion: keep the original settings.** Tier 2's finding holds regardless of
identity threshold; Tier 1's containment-aware setting matches the actual
motivating use case rather than being an arbitrary inflation. Sweep output lives
alongside the main results at `dereplication_pilot/{p110,p97}/{contigs,orfs}/sweep_*/`.

## Interpretation / recommendation

**Real, large redundancy confirmed at both tiers — strong case for building the
full dereplicate-then-annotate pipeline.** Tier 2 in particular implies that, per
participant, the vast majority of funcscan's AMR/AMP/CAZyme gene-level annotation
calls made separately per assembly are duplicated effort: annotating cluster
representatives only and propagating back could plausibly cut that compute by
~10-20x per participant (93-98% redundant → only ~2-7% of genes need independent
annotation).

Tier 1's lower but still substantial redundancy (42-44% by count) means real
savings for geNomad/MAP's contig-level tools too, though less dramatic than Tier 2 —
consistent with the brainstorm's expectation that contig-level and gene-level
redundancy are different, independent quantities.

**This justifies moving past the pilot.** Suggested next steps (not yet started):

1. Stress-test at the larger end (participant 813, 35 single-sample assemblies) to
   confirm the redundancy fraction holds up at bigger scale before committing to
   all 25 participants — cost/benefit could shift if very large participants behave
   differently (e.g. lower redundancy from more strain diversity accumulating over
   more timepoints).
2. Design the actual representative-annotation + propagation pipeline: pick a
   cluster-representative-selection rule (likely just mmseqs's own greedy-mode
   choice, already length-favoring), decide which of MAP/funcscan's tools are safe
   to propagate for (per the brainstorm's still-open "per-tool propagation audit"
   question — nothing checked yet), and design the propagation-back table format
   building on this pilot's `cluster_membership.tsv` shape.
3. Decide retroactive vs. forward-only scope (brainstorm's open question #3) — all
   296 assemblies are already fully annotated, so this could be purely a future-work
   optimization rather than a redo of anything already done.

## Where the real output lives

`dereplication_pilot/{p110,p97}/{contigs,orfs}/`:
- `summary.txt` — the numbers above, per participant/tier.
- `cluster_membership.tsv` — `cluster_id, member_id, source_group, length, cross_group_cluster`,
  the real connective-tissue artifact the brainstorm flagged as not existing yet.
- Large intermediates (pooled input fasta, mmseqs' own `clu_all_seqs.fasta`/
  `clu_rep_seq.fasta`/`clu_cluster.tsv`, `tmp/`) are `.gitignore`d — regeneratable
  via `dereplication_pilot/run_pilot.sh <participant_tag> <coassembly_group>
  <single_group1> [...]` against the real upstream assemblies/gene_export, same
  logic as every other regeneratable-output exclusion in this repo.
