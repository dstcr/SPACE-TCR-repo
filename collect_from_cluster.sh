#!/bin/bash
# =============================================================================
# collect_from_cluster.sh — gather what the repository is missing
#
# Run this ON ARC (a login node is fine; it only reads and writes a tarball).
#
#   bash collect_from_cluster.sh                 # uses config.sh for paths
#   bash collect_from_cluster.sh /some/other/out # write the tarball elsewhere
#
# It collects:
#   1. the SPACE2TCR package      — the clustering library Stage C imports
#   2. conda env specs            — the four environments config.sh names
#   3. module versions            — HMMER, Anaconda, as loaded
#   4. an inventory               — what was found, what was not
#
# It collects NO data: no PDBs, no FASTAs, no distance tables, no input TSVs.
# Nothing here is large, and nothing here is patient data.
#
# INSTALL.md (Track B) explains what each item is used for.
# =============================================================================
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT_BASE="${1:-$HOME}"
STAMP="$(date +%Y%m%d_%H%M%S)"
STAGE="$(mktemp -d "${TMPDIR:-/tmp}/sptce_collect_XXXXXX")"
TARBALL="$OUT_BASE/sptce_missing_${STAMP}.tar.gz"
INV="$STAGE/INVENTORY.txt"

cleanup() { rm -rf "$STAGE"; }
trap cleanup EXIT

mkdir -p "$STAGE/env"
: > "$INV"
note() { echo "$*" | tee -a "$INV"; }
warn() { echo "  [MISSING] $*" | tee -a "$INV"; }

note "SPTCE collection — $(date)"
note "host: $(hostname)   user: ${USER:-unknown}"
note ""

# ---- config -----------------------------------------------------------------
# Honour SPTCE_CONFIG like the job scripts do, and fall back to the repo copy.
: "${SPTCE_CONFIG:=$HERE/config.sh}"
if [[ -f "$SPTCE_CONFIG" ]]; then
  # shellcheck disable=SC1090
  source "$SPTCE_CONFIG"
  note "config: $SPTCE_CONFIG"
  cp "$SPTCE_CONFIG" "$STAGE/config.collected.sh"
else
  note "config: NOT FOUND at $SPTCE_CONFIG — falling back to built-in defaults"
  DATA="${DATA:-$HOME}"
fi
DATA="${DATA:-$HOME}"
note ""

# ---- 1. SPACE2TCR -----------------------------------------------------------
note "== 1. SPACE2TCR package =="
# run_C_slurm.sh uses $DATA/TCR/SPACE2TCR; a few other layouts are tried in case
# the install moved.
SPACE_CANDIDATES=(
  "$DATA/TCR/SPACE2TCR"
  "$DATA/SPACE2TCR"
  "$DATA/TCR/SPACE2TCR/SPACE2TCR"
  "${SPTCE_ROOT:-$HERE}/SPACE2TCR"
)
FOUND_SPACE=""
for c in "${SPACE_CANDIDATES[@]}"; do
  if [[ -d "$c" ]]; then FOUND_SPACE="$c"; break; fi
done

if [[ -n "$FOUND_SPACE" ]]; then
  note "  found: $FOUND_SPACE"
  # Copy the python package only — exclude caches, VCS and any stray data.
  mkdir -p "$STAGE/SPACE2TCR_src"
  tar cf - -C "$(dirname "$FOUND_SPACE")" \
      --exclude='__pycache__' --exclude='*.pyc' --exclude='.git' \
      --exclude='*.pdb' --exclude='*.tsv.gz' --exclude='*.npy' \
      "$(basename "$FOUND_SPACE")" 2>/dev/null | tar xf - -C "$STAGE/SPACE2TCR_src"
  note "  files collected: $(find "$STAGE/SPACE2TCR_src" -type f | wc -l)"
  note "  python modules:"
  find "$STAGE/SPACE2TCR_src" -name '*.py' | sed "s|$STAGE/SPACE2TCR_src/|    |" | tee -a "$INV" >/dev/null
  find "$STAGE/SPACE2TCR_src" -name '*.py' | sed "s|$STAGE/SPACE2TCR_src/|    |"
  # Work out which directory is the importable package, so the extraction
  # instruction at the end names the right level. `from SPACE2TCR.util import`
  # means a directory literally called SPACE2TCR holding util.py must sit on the
  # path, and $DATA/TCR/SPACE2TCR may be either that package or the path root
  # containing it.
  PKG_DIR="$(find "$STAGE/SPACE2TCR_src" -type d -name SPACE2TCR \
               -exec test -f '{}/util.py' \; -print 2>/dev/null | head -1)"
  if [[ -n "$PKG_DIR" ]]; then
    MOVE_FROM="${PKG_DIR#$STAGE/}"
    note "  importable package: $MOVE_FROM  (this is the directory to place at the repo root)"
  else
    MOVE_FROM=""
    warn "no directory named SPACE2TCR containing util.py was found in the copy"
    note "    the import 'from SPACE2TCR.util import ...' will not resolve as collected;"
    note "    check the layout by hand before committing."
  fi
  # Record the interface the pipeline depends on, as it actually resolves here.
  if [[ -n "${PY_SPACE2:-}" && -x "${PY_SPACE2}" ]]; then
    note "  verifying the import against $PY_SPACE2 ..."
    PYTHONPATH="$FOUND_SPACE:${PYTHONPATH:-}" "$PY_SPACE2" - >"$STAGE/env/space2tcr_interface.txt" 2>&1 <<'PY'
import inspect
try:
    from SPACE2TCR import util
except Exception as e:
    print("IMPORT FAILED:", type(e).__name__, e); raise SystemExit(0)
names = ["rmsd", "reg_def", "reg_def_fw_all", "reg_def_CDR_all",
         "cluster_antibodies_by_CDR_length", "same_length_only"]
for n in names:
    ok = hasattr(util, n)
    print(f"{n:34s} {'OK' if ok else 'ABSENT'}")
    if ok:
        o = getattr(util, n)
        if callable(o):
            try: print(f"    signature: {n}{inspect.signature(o)}")
            except Exception: pass
rd = getattr(util, "reg_def", None)
if isinstance(rd, dict):
    print("\nreg_def keys:", sorted(rd.keys()))
print("\nmodule file:", getattr(util, "__file__", "?"))
PY
    sed 's/^/    /' "$STAGE/env/space2tcr_interface.txt" | tee -a "$INV"
  else
    warn "PY_SPACE2 not set or not executable — skipped the import check"
  fi
else
  warn "SPACE2TCR not found. Tried:"
  for c in "${SPACE_CANDIDATES[@]}"; do note "    $c"; done
  note "  Find it with:  find $DATA -maxdepth 4 -type d -name 'SPACE2TCR*'"
  note "  Stage C CANNOT run without this."
fi
note ""

# ---- 2. environments --------------------------------------------------------
note "== 2. conda environments =="
dump_env() {
  local label="$1" py="$2"
  if [[ -z "$py" ]]; then warn "$label: not set in config"; return; fi
  if [[ ! -x "$py" ]]; then warn "$label: $py not executable"; return; fi
  note "  $label -> $py"
  "$py" -V                              > "$STAGE/env/${label}.python-version.txt" 2>&1
  "$py" -m pip freeze                   > "$STAGE/env/${label}.pip-freeze.txt"     2>&1
  local prefix; prefix="$("$py" -c 'import sys; print(sys.prefix)' 2>/dev/null)"
  if [[ -n "$prefix" ]] && command -v conda >/dev/null 2>&1; then
    conda env export -p "$prefix" > "$STAGE/env/${label}.conda-env.yml" 2>/dev/null \
      || warn "$label: conda env export failed (not a conda prefix?)"
  fi
  note "    $(grep -c . "$STAGE/env/${label}.pip-freeze.txt" 2>/dev/null || echo 0) pip packages"
}
dump_env "space2tcr"   "${PY_SPACE2:-}"
dump_env "myenv2TCR"   "${PY_TCRBUILDER:-}"
dump_env "clustcr_env" "${PY_CLUSTCR:-}"

# The stitchr env is named rather than given as a python path, so resolve it.
if [[ -n "${STITCHR_ENV:-}" ]]; then
  note "  stitchr -> ${STITCHR_ENV}"
  if [[ -x "${STITCHR_ENV}/bin/python" ]]; then
    dump_env "stitchr" "${STITCHR_ENV}/bin/python"
  else
    note "    (not a prefix with bin/python; recording the name only)"
    echo "STITCHR_ENV=${STITCHR_ENV}" > "$STAGE/env/stitchr.name.txt"
  fi
  for cli in stitchr thimble; do
    if command -v "$cli" >/dev/null 2>&1; then
      { echo "$cli: $(command -v "$cli")"; "$cli" --version 2>&1 | head -2; } \
        >> "$STAGE/env/stitchr.cli.txt"
    fi
  done
else
  warn "STITCHR_ENV not set in config"
fi
note ""

# ---- 3. modules ------------------------------------------------------------
note "== 3. modules =="
{ module list; } >"$STAGE/env/modules_loaded.txt" 2>&1 || true
for m in HMMER Anaconda3; do
  { module avail "$m"; } >>"$STAGE/env/modules_avail.txt" 2>&1 || true
done
command -v hmmscan >/dev/null 2>&1 && hmmscan -h 2>&1 | head -2 > "$STAGE/env/hmmer_version.txt"
note "  recorded module list and HMMER version where available"
note ""

# ---- 4. inventory ----------------------------------------------------------
note "== 4. repository self-containment check =="
if [[ -n "$FOUND_SPACE" ]]; then
  note "  SPACE2TCR collected -> after extracting, Stage C should import cleanly"
else
  note "  SPACE2TCR NOT collected -> the repository remains incomplete"
fi

tar czf "$TARBALL" -C "$STAGE" . 2>/dev/null
echo
echo "================================================================"
echo "  Wrote: $TARBALL"
echo "  Size:  $(du -h "$TARBALL" | cut -f1)"
echo
if [[ -n "${MOVE_FROM:-}" ]]; then
  echo "  Copy it down, then from the repo root:"
  echo "    tar xzf $(basename "$TARBALL") -C /tmp/sptce_in"
  echo "    cp -r /tmp/sptce_in/$MOVE_FROM  ./SPACE2TCR"
  echo "    cp -r /tmp/sptce_in/env         ./env"
  echo "    PYTHONPATH=. python -c 'from SPACE2TCR.util import rmsd, reg_def; print(sorted(reg_def))'"
  echo "    git add SPACE2TCR env && git commit -m 'Vendor SPACE2TCR and pin environments'"
else
  echo "  SPACE2TCR was not collected in an importable layout — read INVENTORY.txt"
  echo "  before committing anything."
fi
echo
echo "  Read INVENTORY.txt inside the tarball first — it says what was"
echo "  found and what was not."
echo "================================================================"
