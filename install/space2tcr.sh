#!/bin/bash
# =============================================================================
# install/space2tcr.sh — install and verify SPACE2TCR (Stages C and D)
#
#   bash install/space2tcr.sh                 # install into the current python
#   bash install/space2tcr.sh --venv .venv    # make and use a venv
#   bash install/space2tcr.sh --python /path/to/python
#   bash install/space2tcr.sh --check         # verify only, install nothing
#
# SPACE2TCR itself is vendored in this repository at SPACE2TCR/ — there is
# nothing to download and nothing to pip install. This script installs its
# five runtime dependencies and then proves the install works by computing real
# RMSDs, not merely importing the module.
#
# Idempotent: safe to re-run.
# =============================================================================
set -uo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO"

PY="${PY_SPACE2:-}"
VENV=""
CHECK_ONLY=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --venv)   VENV="$2"; shift 2 ;;
    --python) PY="$2"; shift 2 ;;
    --check)  CHECK_ONLY=1; shift ;;
    -h|--help) sed -n '2,18p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0 ;;
    *) echo "unknown argument: $1" >&2; exit 2 ;;
  esac
done

say()  { printf '  %s\n' "$*"; }
ok()   { printf '  [ok]   %s\n' "$*"; }
bad()  { printf '  [FAIL] %s\n' "$*"; FAILED=1; }
FAILED=0

echo
echo "================================================================"
echo "  SPACE2TCR install"
echo "================================================================"

# ---- 1. the vendored package -------------------------------------------------
if [[ -f SPACE2TCR/util.py ]]; then
  ok "vendored package found at $REPO/SPACE2TCR"
else
  bad "SPACE2TCR/util.py is missing — is this the repository root?"
  exit 1
fi

# ---- 2. python ---------------------------------------------------------------
if [[ -n "$VENV" ]]; then
  if [[ ! -d "$VENV" ]]; then
    say "creating virtualenv at $VENV"
    "${PY:-python3}" -m venv "$VENV" || { bad "could not create venv"; exit 1; }
  fi
  PY="$VENV/bin/python"
fi
PY="${PY:-$(command -v python3)}"
[[ -x "$PY" ]] || { bad "no usable python (tried '$PY')"; exit 1; }
ok "python: $("$PY" -V 2>&1)  ($PY)"

PYV=$("$PY" -c 'import sys; print("%d.%d" % sys.version_info[:2])')
case "$PYV" in
  3.9|3.10|3.11|3.12|3.13) ok "python $PYV is supported" ;;
  *) say "[warn] python $PYV is untested here; 3.10-3.13 are known good" ;;
esac

# ---- 3. dependencies ---------------------------------------------------------
if (( ! CHECK_ONLY )); then
  say "installing dependencies from install/requirements-space2tcr.txt ..."
  PIP_ARGS=(--quiet --disable-pip-version-check -r install/requirements-space2tcr.txt)
  # A system python on Debian/Ubuntu refuses to install without this flag; a
  # venv neither needs nor accepts it.
  if [[ -z "$VENV" ]] && "$PY" -c 'import sysconfig,sys; sys.exit(0 if sysconfig.get_config_var("EXT_SUFFIX") and __import__("os").path.exists("/usr/lib/python3/dist-packages") else 1)' 2>/dev/null; then
    "$PY" -m pip install "${PIP_ARGS[@]}" --break-system-packages 2>/dev/null \
      || "$PY" -m pip install "${PIP_ARGS[@]}"
  else
    "$PY" -m pip install "${PIP_ARGS[@]}"
  fi
  [[ $? -eq 0 ]] && ok "dependencies installed" || bad "pip install failed"
fi

MISSING=$("$PY" - <<'EOF'
import importlib.util as u
need = ["numpy", "pandas", "scipy", "numba", "joblib", "sklearn"]
print(" ".join(m for m in need if u.find_spec(m) is None))
EOF
)
if [[ -z "$MISSING" ]]; then
  ok "all 6 runtime modules importable"
  "$PY" - <<'EOF'
import importlib.metadata as md
for p in ("numpy","pandas","scipy","numba","joblib","scikit-learn"):
    try: print(f"         {p:14s} {md.version(p)}")
    except Exception: pass
EOF
else
  bad "missing modules: $MISSING"
fi

# ---- 4. the import the pipeline actually makes --------------------------------
if PYTHONPATH="$REPO" "$PY" -c 'from SPACE2TCR.util import (rmsd, reg_def, reg_def_fw_all, reg_def_CDR_all, cluster_antibodies_by_CDR_length, same_length_only)' 2>/dev/null; then
  ok "all six names import from SPACE2TCR.util"
else
  bad "the SPACE2TCR import failed — Stage C will not run"
  PYTHONPATH="$REPO" "$PY" -c 'import SPACE2TCR.util' 2>&1 | tail -3 | sed 's/^/         /'
fi

# ---- 5. prove it computes, not just imports ----------------------------------
if [[ -f tests/verify_space2tcr.py ]]; then
  say "running tests/verify_space2tcr.py (real RMSDs) ..."
  if OUT=$(PYTHONPATH="$REPO" "$PY" tests/verify_space2tcr.py "$REPO" 2>&1); then
    ok "verification passed: $(grep -c '  PASS  ' <<<"$OUT") checks"
    grep -E 'rigid-body|against itself|perturbing CDR3B' <<<"$OUT" | sed 's/^/       /'
  else
    bad "verification FAILED"
    grep -E 'FAIL' <<<"$OUT" | head -5 | sed 's/^/         /'
  fi
else
  say "[skip] tests/verify_space2tcr.py not present"
fi

# ---- 6. the numba cache note -------------------------------------------------
cat <<'EOF'

  Note on numba: SPACE2TCR compiles two hot functions with
  @nb.njit(cache=True), so SPACE2TCR/__pycache__ holds an AOT cache tied to the
  numba and python versions. It is gitignored and regenerated on first run. If
  you change either version and see a cache error, delete it:
      rm -rf SPACE2TCR/__pycache__
EOF

echo
echo "================================================================"
if (( FAILED )); then
  echo "  Install INCOMPLETE — see the failures above."
  echo "  Dependency list: install/requirements-space2tcr.txt"
else
  echo "  SPACE2TCR is installed and verified."
  echo "  Stage C/D need PYTHONPATH to include the repo root; sbatch/run_C_slurm.sh"
  echo "  sets that for you and preflights the import."
fi
echo "================================================================"
echo
exit $FAILED
