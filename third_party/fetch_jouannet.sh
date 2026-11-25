#!/usr/bin/env bash
# =============================================================================
# fetch_jouannet.sh — obtain the Jouannet benchmark repository
#
#   bash third_party/fetch_jouannet.sh                # clone beside this repo
#   bash third_party/fetch_jouannet.sh ~/benchmarks   # clone somewhere else
#   bash third_party/fetch_jouannet.sh --check <dir>  # verify an existing copy
#
# It is NOT vendored into SPACE-TCR: it carries no licence, it is 382 MB, and
# 292 MB of that is third-party databases with their own terms of use. See
# JOUANNET.md.
#
# Idempotent. Exits non-zero if the result is not usable, so it can gate a
# script. Verifies the eight paths the pipeline actually reads, so a partial
# or interrupted download cannot look like a good one.
# =============================================================================
set -uo pipefail

URL="https://github.com/i3-unit/TCR_Unsupervised_Benchmark"
NAME="TCR_Unsupervised_Benchmark"

# The paths sptce_pipeline.R actually reads.
REQUIRED=(
  "Data/Output_methods"
  "Functions/purity_functions.R"
  "Functions/similarity_functions.R"
  "Functions/subset_functions.R"
  "Analysis/benchmark_tools_pooled_database_CD8.Rmd"
)
# Needed only for some analyses; absence is a warning, not a failure.
OPTIONAL=(
  "Data/Database/database_pooled_human_2023_03_15.txt"   # cohort cross-check
  "Data/10X_Output_methods"                              # Benchmark A
  "Data/10X/Donor1"                                      # Benchmark A
)

say()  { printf '  %s\n' "$*"; }
fail() { printf '\n  FAILED: %s\n\n' "$*" >&2; exit 1; }

check_only=0
if [[ "${1:-}" == "--check" ]]; then check_only=1; shift; fi

# ---- where ------------------------------------------------------------------
here="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
parent="${1:-$(dirname "$here")}"

if [[ $check_only -eq 1 ]]; then
  # --check takes the repo directory itself, or its parent.
  if [[ -d "$parent/$NAME" ]]; then repo="$parent/$NAME"; else repo="$parent"; fi
else
  repo="$parent/$NAME"
fi

echo
echo "================================================================"
echo "  Jouannet benchmark repository"
echo "================================================================"
say "target: $repo"

# ---- fetch ------------------------------------------------------------------
if [[ $check_only -eq 0 ]]; then
  if [[ -d "$repo/.git" ]]; then
    say "already cloned; pulling"
    git -C "$repo" pull --ff-only || say "pull failed, keeping the existing copy"
  elif [[ -d "$repo" ]]; then
    say "directory exists and is not a git clone; leaving it alone"
  else
    command -v git >/dev/null || fail "git is not installed"
    mkdir -p "$parent" || fail "cannot write to $parent"
    say "cloning $URL  (382 MB, several minutes)"
    git clone --depth 1 "$URL" "$repo" || fail \
      "clone failed. Download the ZIP from $URL instead, unpack it, and rerun with --check <dir>."
  fi
fi

[[ -d "$repo" ]] || fail "$repo does not exist. Run without --check to clone it."

# ---- verify -----------------------------------------------------------------
echo
say "verifying the paths the pipeline reads:"
missing=0
for p in "${REQUIRED[@]}"; do
  if [[ -e "$repo/$p" ]]; then say "  [ok] $p"
  else                         say "  [--] $p   REQUIRED"; missing=$((missing + 1)); fi
done
for p in "${OPTIONAL[@]}"; do
  if [[ -e "$repo/$p" ]]; then say "  [ok] $p"
  else                         say "  [--] $p   optional"; fi
done

# Output_methods should hold the nine methods' cluster output.
if [[ -d "$repo/Data/Output_methods" ]]; then
  n=$(find "$repo/Data/Output_methods" -maxdepth 1 -type f | wc -l | tr -d ' ')
  say ""
  say "  Data/Output_methods: $n files"
  [[ "$n" -lt 9 ]] && say "  expected at least 9 (one per published method) — this copy looks partial"
fi

if [[ $missing -gt 0 ]]; then
  echo
  say "$missing required path(s) absent. This copy is not usable."
  fail "incomplete copy at $repo"
fi

sz=$(du -sh "$repo" 2>/dev/null | cut -f1)
echo
say "OK — complete copy, $sz"
echo
echo "  Now point the pipeline at it:"
echo
echo "    sptce_pipeline.R    CFG\$jouannet_repo <- \"$repo\""
echo
echo "  Cite Jouannet et al. (2025), and IEDB, VDJdb and McPAS-TCR, whose data"
echo "  the pooled database is built from."
echo "================================================================"
echo
