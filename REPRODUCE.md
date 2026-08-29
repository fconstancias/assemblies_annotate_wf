# Reproducing this annotation workflow on a new cluster

Procedural setup/run guide — not a narrative. For the *why* behind any of this, see
`CLAUDE.md` (full decision history) and `coassembly_production/README.md` (production-scale
run specifically). This document assumes you already have, per sample/assembly: an anvi'o
**contigs-db with gene calls already done** (`anvi-gen-contigs-database` + prodigal/
pyrodigal-gv gene calling). That's the one hard prerequisite everything else builds on.

## 1. Software prerequisites

- **conda/mamba** (mamba strongly preferred — plain `conda create`/`conda search` can be
  10-100x slower on some HPC filesystems; genuinely worth installing if not already present).
- **Nextflow >= 25.10.4** (funcscan 4.0.0's floor). If your conda-packaged nextflow is older,
  `nextflow self-update` fixes it in place (jumps to whatever the latest release is).
- **Singularity or Apptainer**, with a writable cache directory somewhere with real space
  (funcscan+MAP+their DBs run into the tens of GB of pulled containers).
- **SLURM** (or adapt `nextflow.config`'s `executor`/`queue`/`clusterOptions` to your
  scheduler — everything else is scheduler-agnostic).
- A shell that actually initializes conda/modules for **non-interactive** sessions — many
  HPC `.bashrc` files gate conda-init behind an `if [[ $- == *i* ]]` interactive-only check
  that a script/CI/agent shell never satisfies. Source directly instead of relying on
  `.bashrc`: `source <conda_root>/etc/profile.d/conda.sh`.

## 2. Conda environments needed

| Env | Purpose | Build |
|---|---|---|
| `anvio-9` (or your anvi'o version) | contigs-db work, gene export, KEGG/KOfam metabolism | however you normally provision anvi'o — not covered here, assumed already present |
| `env_nf` | dedicated Nextflow env | `mamba create -n env_nf -c bioconda -c conda-forge nextflow`, then `nextflow self-update` if its floor is too old |
| `gff2gbk` | gff3+fna → real GenBank conversion (only needed if AMPCOMBI2's cross-sample clustering matters to you — see §6) | `mamba create -n gff2gbk -c conda-forge python=3.10 biopython` then `pip install --no-deps bcbio-gff six` — **do not** `conda install bcbio-gff` directly, its bioconda package pins an incompatible old biopython that lacks the `SimpleLocation` API its own code needs |

## 3. Pipeline checkouts + required patches

Both pipelines are used via **remote-form invocation** (`nextflow run ORG/repo -r TAG`) —
never a local clone path (switching between the two invalidates the entire resume cache).
Nextflow stages the real checkout at `~/.nextflow/assets/ORG/repo/` on first run; patch that
directly (it's a real git repo — `git diff`/`git stash`/`git apply --check` all work
normally).

**funcscan** (`nf-core/funcscan`, pin `4.0.0`):
1. `assets/schema_input.json` — drop `"protein": ["gbk"]` from `dependentRequired` (schema
   wrongly demands a gbk file even when gff is also supplied; the actual workflow code
   treats gbk as genuinely optional).
2. `modules/nf-core/ampcombi2/parsetables/main.nf` — `ampcombi`'s own validation rejects a
   fully-absent or bare `--gbk`; needs a placeholder `.gbff` generated *with Biopython
   itself* (a hand-written LOCUS line fails Biopython's own strict parser) when no real gbk
   is supplied, but passed straight through unmodified when a real gbk file *is* supplied
   (an earlier version of this patch always treated `--gbk`'s value as a directory to
   `mkdir`/`ls` into, which crashes on a real file input).
3. `bin/ampcombi_download.py` — DRAMP database download crashes on any row with an empty/NaN
   `Sequence` cell (`re.match()` on a float). Skip non-string rows.

Reference patches (regenerate against your own fresh checkout rather than assuming these
apply byte-for-byte — pipeline point-releases can drift):
`funcscan_gff_test/patches/{01_schema_drop_protein_requires_gbk,02_ampcombi2_parsetables_optional_gbk,03_ampcombi_download_dramp_nan_fix}.patch`

**MAP** (`EBI-Metagenomics/mobilome-annotation-pipeline`, pin `v5.0.0`):
- No code patches needed, but two upstream bugs need `errorStrategy = 'ignore'` overrides in
  your own `nextflow.config` (not the checkout itself) for `--download_dbs` to complete —
  see `map_run/nextflow.config` and `map_run/README.md` for both, plus the by-hand DB rescue
  each needs (`VFDB_setB_pro.dmnd` built manually; geNomad's tarball extracts to the wrong
  folder name than the pipeline expects).
- `params.run_native_signalp` must be explicitly set (even to `false`) in your config or the
  pipeline errors on an undefined param — this only applies if your checkout still carries
  any custom SignalP patches from prior work; a clean checkout may not need it.

## 4. Reference databases

- **anvi'o KEGG/KOfam** (for `anvi-run-kegg-kofams`/`anvi-estimate-metabolism`): one-time
  `anvi-setup-kegg-data --mode all --kegg-data-dir <path>` — ~10GB, no licensing gate, no
  login required. Point `--kegg-data-dir` at a real path yourself; if you pre-`mkdir` it,
  anvi'o refuses to touch a pre-existing directory as a safety check — let it create the
  directory itself.
- **funcscan's own DBs** (AMRFinderPlus, CARD/RGI, DeepARG, dbCAN, DRAMP): auto-downloaded
  on first run, no action needed beyond disk space.
- **MAP's own DBs** (`--download_dbs`): mostly auto-downloaded; two need the by-hand rescue
  in §3.
- **DRAM2** (only if pursuing it — optional, see `CLAUDE.md`'s DRAM2-vs-anvi'o-KEGG
  writeup): databases are Globus-only (UUID `97ed64b9-dea0-4eb5-a7e0-0b50ac94e889`,
  `https://app.globus.org/file-manager?origin_id=97ed64b9-dea0-4eb5-a7e0-0b50ac94e889`),
  requires a human to complete the browser-based login — not automatable. Licensed KEGG is
  explicitly excluded from that bundle (KOfam + the other free databases only).

## 5. Per-sample gene export (from an existing contigs-db)

```bash
conda activate anvio-9
DB=/path/to/sample.db          # already has gene calls
OUT=/path/to/export_dir

# GFF3 (default deflines are fine here — GFF3 mode rejects --defline-format anyway)
anvi-get-sequences-for-gene-calls -c "$DB" --export-gff3 -o "$OUT/sample.gff3"

# FAA — the ID-matching fix that everything downstream depends on
anvi-get-sequences-for-gene-calls -c "$DB" --get-aa-sequences \
    --defline-format "{contigs_db_project_name}___{gene_caller_id}" -o "$OUT/sample.faa"

# FNA (nucleotide gene sequences — only needed for DRAM2-style pre-called-gene input)
anvi-get-sequences-for-gene-calls -c "$DB" \
    --defline-format "{contigs_db_project_name}___{gene_caller_id}" -o "$OUT/sample.fna"
```

**Why the `--defline-format` matters**: `--export-gff3`'s `ID=` attribute and
`--get-aa-sequences`'s default defline use *different* ID schemes by default (bare
`{gene_caller_id}` vs `{contigs_db_project_name}___{gene_caller_id}`) — silently mismatched
IDs between your gff3 and faa otherwise, with no error, just wrong/missing joins downstream.
`--export-gff3` itself rejects `--defline-format` (incompatible), so only the FAA/FNA calls
need it; GFF3's `ID=` already uses the `{project}___{gene_id}` scheme by default.

Real per-file CDS-count = FAA-sequence-count is the check that actually catches a drift here
— don't just trust "the command didn't error."

### Optional: real GenBank (`.gbff`) for full AMPCOMBI2 harmonization

Only needed if AMPCOMBI2's cross-sample clustering step (`AMPCOMBI2_CLUSTER`, MMseqs2-based
AMP dedup) matters — per-sample AMP prediction works fine without it. Requires the `gff2gbk`
env from §2:

```python
# gff3_to_gbk.py — see assembly_annotation_wf/scripts/ for the working version
from BCBio import GFF
from Bio import SeqIO

base_dict = SeqIO.to_dict(SeqIO.parse(fasta_path, "fasta"))
records = []
for record in GFF.parse(gff3_path, base_dict=base_dict):
    record.annotations["molecule_type"] = "DNA"
    for feature in record.features:
        if feature.type == "CDS" and "ID" in feature.qualifiers:
            feature.qualifiers["locus_tag"] = feature.qualifiers["ID"]  # must match FAA header
    records.append(record)
SeqIO.write(records, output_path, "genbank")
```

Verify it round-trips through the *actual* `ampcombi` container's Biopython before trusting
it (`singularity exec <ampcombi_container> python3 -c "from Bio import SeqIO; ..."`) — a
Biopython version mismatch between whatever wrote the file and whatever reads it is exactly
the kind of thing that fails silently otherwise.

## 6. Samplesheets + launch commands

**funcscan** (`sample,fasta,protein,gff,gff_type[,gbk]` — `gbk` column optional, see §5):

```bash
nextflow run nf-core/funcscan -r 4.0.0 \
    --input samplesheet.csv --outdir results \
    --run_arg_screening --run_cazyme_screening \
    --run_amp_screening --amp_skip_amplify --amp_skip_macrel \
    -c nextflow.config -profile singularity -resume <run_name>
```
`gff_type` is a fixed enum (`NCBI_prok`/`prodigal`/`NCBI_euk`/`JGI`, no generic parser) —
use `prodigal` for anvi'o's minimal GFF3 export, verified working even with the bare `ID=`
attribute set (none of real Prodigal's `partial=`/`start_type=`/etc. needed). **No BGC
flag** — funcscan's BGC tools (antiSMASH/DeepBGC/GECCO) all require real GenBank content our
export doesn't produce; see `CLAUDE.md` for the full reasoning if BGC is ever reconsidered.

**MAP** (`sample,assembly,proteins_gff,proteins_faa,virify_gff,interproscan_tsv` — last two
optional/empty):

```bash
nextflow run EBI-Metagenomics/mobilome-annotation-pipeline -r v5.0.0 \
    --input samplesheet.csv --outdir results \
    -c nextflow.config -c my_paths.config \
    -profile singularity -resume <run_name>
```
Set `skip_sanntis`/`skip_gecco`/`skip_antismash` as real booleans in `nextflow.config`'s
`params {}` block, **not** as `--skip_x true/false` CLI flags — depending on your Nextflow
version, the CLI parser may reject these as `[string]` instead of `[boolean]` (a real
regression hit mid-session after a `nextflow self-update`). `--skip_sanntis false` (needs
real InterProScan as a prerequisite — that IPS run also gives real SignalP4.1-based export
prediction for virulence factors "for free", see `CLAUDE.md`), `--skip_gecco true
--skip_antismash true` (BGC is funcscan's job in this split — see above).

## 7. anvi'o-KEGG (metabolism estimation)

```bash
conda activate anvio-9
anvi-run-kegg-kofams -c sample.db --kegg-data-dir <kegg_dir> -T <threads>
anvi-estimate-metabolism -c sample.db --kegg-data-dir <kegg_dir> \
    --metagenome-mode -O sample_metabolism
```
`--metagenome-mode` is for an unbinned assembly (pools all gene annotations community-wide —
real signal, but completeness can be inflated since multiple populations contribute
enzymes). Once bins/MAGs exist as an anvi'o collection, switch to internal-genomes-file /
profile-db+collection mode for real per-MAG completeness — no re-export needed, same
contigs-db.

At real per-sample scale (hundreds of thousands of genes), run this via your scheduler
(`sbatch`, one task per sample/array index), not directly on a login/interactive node — see
`coassembly_production/anvio_kegg/run_kofams_array.sh` for a working SLURM array template.

## 8. Config gotchas to apply proactively, not rediscover

All confirmed the hard way this session — full detail in `CLAUDE.md`'s Nextflow operational
gotchas section. Apply these up front on a new cluster rather than waiting to hit them:

- `export TMPDIR=/tmp` (or any path Singularity actually auto-binds) before every
  `nextflow run` — an inherited `TMPDIR` pointing somewhere containers can't see breaks
  tools that create real temp files, silently and confusingly.
- Set a real `executor.queueSize` for your actual scale — don't copy a small-test value into
  a production run without checking real idle cluster capacity (`sinfo -p <queue> -o "%D
  %C"`). This is safe to change and `-resume` freely (executor-level, not part of the cache
  hash).
- Build in a retry policy from the start for production-scale runs — transient container/
  infra issues (loop-device exhaustion) and real OOMs (undersized defaults at real scale)
  both showed up only once actual production-scale data was used, never during small-subset
  validation. Some OOMs surface as a masked plain exit 1 (SIGKILL caught and re-raised by
  the tool's own wrapper), not the expected 137 — broaden `errorStrategy` accordingly if a
  run needs to be resilient, e.g.:
  ```groovy
  process {
      errorStrategy = { task.exitStatus in ([1] + (130..145) + [104, 255]) ? 'retry' : 'finish' }
      maxRetries    = 4
      withName: '.*' { memory = { [16.GB * task.attempt, 128.GB].min() } }  // fixed base only, never task.memory — see below
  }
  ```
- **Never reference `task.memory` inside a closure that's scaling `memory` itself** —
  Nextflow resolves `task.memory` by calling that same closure, so it recurses infinitely
  and crashes the run (`java.lang.StackOverflowError`). Always use a fixed literal base.
- Editing any `process{}`-level directive (as opposed to `executor{}`-level) invalidates the
  resume cache for *every* process, not just the one changed — expect a full recompute after
  such an edit, not a true resume.
