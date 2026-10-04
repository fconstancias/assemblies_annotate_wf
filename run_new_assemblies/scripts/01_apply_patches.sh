#!/usr/bin/env bash
# One-time (idempotent): pull funcscan 4.0.0 + MAP v5.0.0 and apply patches/ to the checkouts Nextflow uses
# (~/.nextflow/assets/...). Never use a separate clone: launching from a local path invalidates the resume cache.
# Prints "applied" / "already applied"; stops on a patch that applies neither way (checkout drifted).
set -euo pipefail
D=$(cd "$(dirname "$0")/.." && pwd)
source "$D/scripts/00_env.sh"
apply() {  # checkout dir, patch dir
  for p in "$2"/*.patch; do
    if git -C "$1" apply --reverse --check "$p" 2>/dev/null; then echo "  already applied  $(basename "$p")"
    elif git -C "$1" apply --check "$p" 2>/dev/null; then git -C "$1" apply "$p"; echo "  applied          $(basename "$p")"
    else echo "  FAILED           $(basename "$p") (checkout differs from the pinned version?)"; exit 1; fi
  done
}
FS=$HOME/.nextflow/assets/nf-core/funcscan
MAP=$HOME/.nextflow/assets/EBI-Metagenomics/mobilome-annotation-pipeline
[ -d "$FS" ]  || nextflow pull nf-core/funcscan -r 4.0.0
[ -d "$MAP" ] || nextflow pull EBI-Metagenomics/mobilome-annotation-pipeline -r v5.0.0
echo "funcscan ($FS, $(git -C "$FS" describe --tags 2>/dev/null))"; apply "$FS" "$D/patches/funcscan"
echo "MAP ($MAP, $(git -C "$MAP" describe --tags 2>/dev/null))";     apply "$MAP" "$D/patches/map"
