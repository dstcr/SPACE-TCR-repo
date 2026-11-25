#!/usr/bin/env bash
# =============================================================================
# SPACE-TCR — one-command install
#
#   ./install.sh                 reproduce the paper (R only, ~5 min)
#   ./install.sh --own-data      the above, plus structural clustering
#   ./install.sh --everything    the above, plus the ML distillation
#   ./install.sh --check         verify an existing install, change nothing
#
# This wraps the installers under install/, which remain usable on their own.
# Everything here is idempotent: re-running is safe. Any failure exits non-zero
# with the command that failed, so this can gate a CI job.
#
# NOT installed by this script, because each needs a conda environment and a
# licence acceptance that cannot be automated:
#   * TCRBuilder2  (structure prediction, Stage B)  -> install/TCRBUILDER2.md
#   * stitchr      (sequence reconstruction, Stage A) -> INSTALL.md
# Neither is needed to reproduce the paper: the predicted geometry is committed.
# =============================================================================
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE"

MODE=paper
CHECK_ONLY=0
for arg in "$@"; do
  case "$arg" in
    --own-data)   MODE=own ;;
    --everything) MODE=all ;;
    --check)      CHECK_ONLY=1 ;;
    -h|--help)    sed -n '2,18p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
    *) echo "unknown option: $arg (try --help)" >&2; exit 2 ;;
  esac
done

bold() { printf '\n\033[1m%s\033[0m\n' "$*"; }
ok()   { printf '  \033[32m[ok]\033[0m   %s\n' "$*"; }
die()  { printf '  \033[31m[fail]\033[0m %s\n' "$*" >&2; exit 1; }

need() {
  command -v "$1" >/dev/null 2>&1 || die "$1 not found. $2"
}

bold "SPACE-TCR install  (mode: $MODE)"

# ---- prerequisites ---------------------------------------------------------
need Rscript "Install R >= 4.1 — see INSTALL.md"
R_VER=$(Rscript -e 'cat(paste(R.version$major, R.version$minor, sep="."))')
ok "R $R_VER"

if [ "$MODE" != paper ]; then
  need python3 "Install Python 3.10-3.13 — see INSTALL.md"
  PY_VER=$(python3 -c 'import sys; print("%d.%d"%sys.version_info[:2])')
  ok "Python $PY_VER"
fi

# ---- verify-only path ------------------------------------------------------
if [ "$CHECK_ONLY" = 1 ]; then
  bold "Verifying (nothing will be installed)"
  Rscript check_install.R || die "check_install.R failed"
  if [ "$MODE" != paper ]; then
    bash check_install.sh || die "check_install.sh failed"
  fi
  bold "All checks passed."
  exit 0
fi

# ---- 1. R ------------------------------------------------------------------
bold "1/3  R packages"
if [ "$MODE" = paper ]; then
  Rscript install/r_packages.R || die "Rscript install/r_packages.R"
else
  Rscript install/r_packages.R --figures || die "Rscript install/r_packages.R --figures"
fi
ok "R packages installed"

# ---- 2. structural clustering ---------------------------------------------
if [ "$MODE" != paper ]; then
  bold "2/3  SPACE2-TCR (structural clustering)"
  bash install/space2tcr.sh || die "bash install/space2tcr.sh"
  ok "SPACE2-TCR installed and RMSD invariance verified"
else
  bold "2/3  SPACE2-TCR — skipped"
  echo "     Not needed to reproduce the paper; the geometry is committed."
  echo "     Run ./install.sh --own-data to add it."
fi

# ---- 3. ML distillation ----------------------------------------------------
if [ "$MODE" = all ]; then
  bold "3/3  Distillation model"
  bash install/distill.sh --venv .venv-ml || die "bash install/distill.sh"
  ok "distillation environment installed"
else
  bold "3/3  Distillation model — skipped"
  echo "     Run ./install.sh --everything to add it."
fi

# ---- verify ----------------------------------------------------------------
bold "Verifying"
Rscript check_install.R || die "check_install.R failed after install"
if [ "$MODE" != paper ]; then
  bash check_install.sh || die "check_install.sh failed after install"
fi

bold "Done. Next:"
cat <<'EOF'
  Rscript validate_outputs.R        reproduce the paper's headline numbers
                                    (under a minute, 36 assertions, exits non-zero on
                                    any mismatch)

  Then see README.md for where to go next.
EOF
