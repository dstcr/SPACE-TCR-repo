#!/bin/bash
# =============================================================================
# prepare_figure_data.sh — shrink the figure inputs to something portable
#
#   bash prepare_figure_data.sh <sptce_output_dir> <labelled_input.tsv> [out_dir]
#
# e.g.
#   bash prepare_figure_data.sh ~/sptce_output data/jouannet_labelled_input.tsv
#
# The distance tables hold every pairwise RMSD in every CDR-length bin, which is
# quadratic in bin size. The pipeline reads only pairs at or below its loosest
# threshold (max(CFG$thresholds), 1.25 A), so dropping the rest changes no
# figure or number. Set RMSD_MAX to keep a higher ceiling. Structural
# near-neighbours are rare, so this typically shrinks the tables by one to
# three orders of magnitude.
# =============================================================================
set -euo pipefail

SRC="${1:?usage: prepare_figure_data.sh <sptce_output_dir> <labelled_input.tsv> [out_dir]}"
LAB="${2:?missing labelled input TSV}"
OUT="${3:-./figure_data}"
RMSD_MAX="${RMSD_MAX:-1.25}"

[[ -d "$SRC" ]] || { echo "not a directory: $SRC" >&2; exit 1; }
[[ -f "$LAB" ]] || { echo "not a file: $LAB" >&2; exit 1; }

DIST_SRC="$SRC/distance_tables"
[[ -d "$DIST_SRC" ]] || { echo "no distance_tables under $SRC" >&2; exit 1; }

mkdir -p "$OUT/sptce_out/distance_tables"
echo "  source:   $SRC"
echo "  ceiling:  rmsd <= $RMSD_MAX A"
echo "  output:   $OUT"
echo

before=$(du -sk "$DIST_SRC" | cut -f1)
n=0; kept_total=0; read_total=0

shopt -s nullglob
for f in "$DIST_SRC"/*.tsv.gz; do
  b="$(basename "$f")"
  # Keep the header, then only rows whose rmsd column is within the ceiling.
  # The rmsd column is located by name, so column order does not matter.
  out_rows=$(
    gzip -dc "$f" | awk -v max="$RMSD_MAX" -F'\t' '
      NR==1 { for (i=1;i<=NF;i++) if ($i=="rmsd") c=i; print; next }
      c && $c != "" && $c+0 <= max { print; k++ }
      END { print k+0 > "/dev/stderr" }
    ' 2>>"$OUT/.counts" | gzip -c > "$OUT/sptce_out/distance_tables/$b"; echo done
  )
  rows=$(gzip -dc "$f" | wc -l)
  kept=$(gzip -dc "$OUT/sptce_out/distance_tables/$b" | wc -l)
  read_total=$((read_total + rows)); kept_total=$((kept_total + kept))
  n=$((n+1))
  printf "  %-34s %9d -> %8d rows\n" "$b" "$((rows-1))" "$((kept-1))"
done
rm -f "$OUT/.counts"

[[ $n -gt 0 ]] || { echo "no .tsv.gz files found in $DIST_SRC" >&2; exit 1; }

# The dihedral table is per-residue, not pairwise, so it is already small.
if [[ -f "$SRC/dihedrals.tsv.gz" ]]; then
  cp "$SRC/dihedrals.tsv.gz" "$OUT/sptce_out/"
  echo "  copied dihedrals.tsv.gz"
else
  echo "  [warn] no dihedrals.tsv.gz — the gated modes will be skipped"
fi

cp "$LAB" "$OUT/labelled_input.tsv"
echo "  copied $(basename "$LAB") -> labelled_input.tsv"

after=$(du -sk "$OUT" | cut -f1)
echo
echo "================================================================"
printf "  distance tables: %d files, %d -> %d rows (%.1f%% kept)\n" \
       "$n" "$((read_total-n))" "$((kept_total-n))" \
       "$(awk -v a="$kept_total" -v b="$read_total" 'BEGIN{printf "%.1f", 100*a/b}')"
printf "  on disk:         %s -> %s\n" \
       "$(numfmt --to=iec $((before*1024)) 2>/dev/null || echo "${before}K")" \
       "$(numfmt --to=iec $((after*1024))  2>/dev/null || echo "${after}K")"
echo
echo "  Now make one archive and send that:"
echo "    tar czf figure_data.tar.gz -C \"$OUT\" ."
echo
echo "  Every figure in the paper can be rendered from it: pairs above"
echo "  $RMSD_MAX A are dropped at read time in every run anyway."
echo "================================================================"
