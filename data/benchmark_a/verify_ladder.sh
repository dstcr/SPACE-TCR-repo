#!/usr/bin/env bash
# =============================================================================
# verify_ladder.sh — check how many DISTINCT rungs a noise ladder actually has
#
#   bash verify_ladder.sh              # the tracked tables in this directory
#   bash verify_ladder.sh <dir>        # a fresh set of donor1_subset*_sptce.tsv
#
# Run it on any ladder output before plotting it. Two rungs built from the same
# receptors are one rung, however different their row order makes them look:
# identical file sizes are easy to miss and md5 of the unsorted file will not
# catch it.
#
# Exits 0 if every rung is distinct, 1 if any two coincide.
# =============================================================================
set -uo pipefail

dir="${1:-$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)}"
cd "$dir" || { echo "no such directory: $dir" >&2; exit 2; }

shopt -s nullglob
files=(donor1_subset*_sptce.tsv)
if [[ ${#files[@]} -le 1 ]]; then
  cat <<'EOF' >&2
This directory holds the ladder as ONE table plus a membership column, which is
the tracked layout: donor1_subset6_sptce.tsv is the last rung and
ladder_membership.tsv gives each receptor the first rung it appears in. There
are no separate per-rung files for this script to compare.

Materialise them first, then point this script at the result:

    Rscript rebuild_rungs.R build && bash verify_ladder.sh build

That round trip is itself the check worth running: rebuild_rungs.R asserts the
rung sizes, and this script then asserts the rungs are distinct.
EOF
  exit 2
fi

echo
echo "================================================================"
echo "  noise ladder: $dir"
echo "================================================================"
echo
printf '  %-28s %9s %9s  %s\n' file rows uniq-ids digest
echo "  ----------------------------------------------------------------"

declare -a names digests
for f in "${files[@]}"; do
  rows=$(($(wc -l < "$f") - 1))
  # Hash the receptor id set, order-independent. Sorting the whole table would
  # also catch a column reorder, but the id set is what defines a rung.
  dg=$(tail -n +2 "$f" | cut -f1 | sort -u | md5sum | cut -c1-16)
  uq=$(tail -n +2 "$f" | cut -f1 | sort -u | wc -l | tr -d ' ')
  printf '  %-28s %9s %9s  %s\n' "$f" "$rows" "$uq" "$dg"
  [[ "$rows" != "$uq" ]] && printf '  %-28s %s\n' "" "WARNING: $((rows - uq)) duplicate id(s)"
  names+=("$f"); digests+=("$dg")
done

# ---- collisions -------------------------------------------------------------
echo
collisions=0
n=${#names[@]}
for ((i = 0; i < n; i++)); do
  for ((j = i + 1; j < n; j++)); do
    if [[ "${digests[i]}" == "${digests[j]}" ]]; then
      echo "  IDENTICAL: ${names[i]}  ==  ${names[j]}"
      collisions=$((collisions + 1))
    fi
  done
done

# distinct rungs = distinct digests
distinct=$(printf '%s\n' "${digests[@]}" | sort -u | wc -l | tr -d ' ')

echo
echo "  $n table(s) on disk, $distinct distinct rung(s)"

if [[ $collisions -gt 0 ]]; then
  cat <<'EOF'

  Some rungs are built from the same receptors. Plotting them as separate
  noise levels shows a plateau the data does not demonstrate, and their
  spread is a reshuffling artefact rather than sampling variation.

  Plot the distinct rungs, or — only if the CLUSTERING was re-run per rung
  with a different seed — label them as replicates measuring algorithmic
  stability at one noise level. The subset tables alone cannot tell you
  which of those is true.

  See ../../BENCHMARK_A.md
EOF
  echo "================================================================"
  echo
  exit 1
fi

echo "  Every rung is distinct."
echo "================================================================"
echo
