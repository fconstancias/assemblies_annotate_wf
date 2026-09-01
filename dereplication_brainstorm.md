# Dereplication + annotation — brainstorm, not a decided plan

Started 2026-09-01, prompted by real pain this session: annotating every contig/gene in
every one of 296 assemblies (19 co-assemblies + 277 single-sample groups) means re-running
the same heavy tools (RGI, SANNTIS, IPS, dbCAN, geNomad, ICEfinder...) on what is very
likely, for a longitudinal series from the same 25 participants, a lot of near-identical
sequence. Two independent dereplication axes, not one — contigs and ORFs answer different
questions and are annotated by different tool sets.

**Also the likely right foundation for the deferred coverage/quantification goal** (item 4,
explicitly "not yet" this session) — a dereplicated catalog + a cluster-membership table
mapping back to (sample, contig/gene) *is* the standard gene-catalog/abundance-profiling
architecture (same shape as IGC/GMGC/UHGP-style studies). Worth treating this brainstorm and
that eventual work as the same design effort, not two separate ones.

## Tier 1 — contig-level dereplication + contig-level annotation (geNomad, etc.)

**What it saves**: geNomad (plasmid/virus classification), and by the same logic
ICEfinder2, IntegronFinder, ISEScan, PathoFact2's virulence context — every MAP tool that
classifies a *whole contig*, not a gene. Gene-level clustering (Tier 2) does nothing for
this redundancy — a shared plasmid gets reclassified by geNomad once per sample regardless
of how deduplicated its genes are.

**Method sketch**:
1. Pool contigs across all 296 assemblies (probably only ≥1kb or ≥5kb — matching MAP's own
   `RENAME` length tiers already used everywhere else in this project; tiny contigs are
   both less likely to be structurally interesting and expensive to cluster at scale for
   little benefit).
2. `mmseqs easy-cluster` on the nucleotide sequences directly. Real technical point already
   worked out in conversation: use a **containment-aware `--cov-mode`** (target- or
   query-coverage, not the symmetric default), because the same underlying plasmid/region
   often assembles to *different lengths* in different samples depending on coverage depth
   — naive bidirectional coverage would fail to cluster a partial fragment with a fuller
   assembly of the same thing elsewhere. Identity threshold should be tight (~95-99%) since
   this is the same participant population over time, not a cross-study comparison — real
   near-identity is plausible and expected, not just homology.
3. Pick a representative per cluster — favor the longest/most-complete member (more
   informative input for geNomad's own classification than a fragment would be).
4. Run geNomad (+ the other contig-level MAP tools) on representatives only.
5. Propagate each cluster's calls back to every member contig via the cluster-membership
   table.

**Open risk, not yet checked**: does every contig-level MAP output actually depend on
sequence alone? If any tool's call is influenced by sample-specific signal (e.g. anything
coverage/depth-based, if `COMPOSITIONAL_OUTLIER_DETECTION` or similar turns out to use more
than raw sequence) propagating a single annotation to all cluster members would be wrong
for that specific output — needs an audit before trusting propagation universally, tool by
tool, not assumed.

## Tier 2 — ORF/gene-level dereplication + gene-level annotation

**What it saves**: RGI, AMRFinderPlus, DeepARG, AMPIR/AMPCOMBI2, dbCAN, and IPS's domain
work (SANNTIS's own prerequisite) — every tool that scores an individual protein.

**Method sketch**:
1. Pool the `.faa` exports already produced for every sample (this catalog essentially
   already exists as a byproduct of this session's gene_export/ work, just not yet merged
   across samples).
2. `mmseqs easy-cluster` in protein space — the well-trodden path here (this is exactly how
   IGC/GMGC/UniRef-style catalogs are built). Looser thresholds than Tier 1 are
   conventional (~90-95% identity / ~80-90% coverage) since protein clustering usually
   tolerates synonymous/near-synonymous variation, but given the same-population context
   here, tighter thresholds closer to "same gene copy" rather than "same gene family" are
   worth considering too — an actual decision point, not obvious which is right without
   looking at real cluster size distributions at both settings.
3. Annotate representatives only with the full funcscan + MAP AMR/AMP/CAZyme stack.
4. Propagate back to every member ORF via cluster membership.

**Real engineering gap, doesn't exist yet**: the cluster-membership → (sample, gene_id)
→ annotation propagation itself. Needs to preserve enough structure that downstream
coverage/quantification can still work per original sample — the whole point of tracking
is *longitudinal* signal per participant, so the mapping table has to carry
`cluster_id ↔ {sample, gene_caller_id}`, not just collapse it away.

## Relationship between the two tiers

Independent, parallel dereplication passes over the same underlying (sample, contig, gene)
objects — not sequential. A contig's chosen Tier-1 representative doesn't determine which
of its genes get chosen as Tier-2 representatives (gene clustering pools *all* genes across
*all* samples, not just those on Tier-1's chosen contigs). Two separate cluster-membership
tables, each funneling its own tool set's results back to the same real per-sample objects.
This matches how real gene-catalog studies structure things — a gene catalog is normally
built independently of any contig-level analysis, not derived from it.

## Open questions worth resolving before building anything

1. **Threshold tuning** — no obvious universally-right identity/coverage cutoff for
   either tier without actually looking at real cluster-size distributions at a few
   candidate settings on real data first.
2. **Clustering compute cost itself** — potentially millions of contigs/genes across 296
   samples; MMseqs2 clustering at that scale is far cheaper than the annotation it replaces,
   but not free, and worth planning resources for deliberately rather than discovering the
   hard way (same OOM/timeout pattern this whole session already hit repeatedly for the
   annotation tools themselves).
3. **Retroactive vs. forward-only** — co-assembly work is done, single-sample MAP is
   already running. Does this apply only to *future* annotation work, or is it worth a
   retroactive redundancy analysis (how much of what was already run was actually
   duplicate sequence) even without redoing anything already complete? The number itself
   would be informative regardless of whether we act on it.
4. **Per-tool propagation audit** (Tier 1's open risk, generalized) — before trusting
   propagation for *any* tool's output, confirm it's a pure function of sequence, not
   something that could reasonably vary by sample/coverage context even for identical
   input sequence.
