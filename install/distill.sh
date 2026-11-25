#!/bin/bash
# =============================================================================
# install/distill.sh — install and verify the distillation model (ml/)
#
#   bash install/distill.sh                 # install into the current python
#   bash install/distill.sh --venv .venv-ml # make and use a venv (recommended)
#   bash install/distill.sh --python /path/to/python
#   bash install/distill.sh --check         # verify only, install nothing
#   bash install/distill.sh --no-selftest   # skip the training self-test
#
# Installs torch, numpy, pandas and rapidfuzz, then proves the install by
# training the model on synthetic data and checking it beats the Levenshtein
# baseline. CPU is enough — the self-test takes about a minute.
#
# A venv is worth it here: torch is large and pins numpy, which can disturb
# other things sharing a system python.
#
# Idempotent: safe to re-run.
# =============================================================================
set -uo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO"

PY=""
VENV=""
CHECK_ONLY=0
SELFTEST=1
while [[ $# -gt 0 ]]; do
  case "$1" in
    --venv)        VENV="$2"; shift 2 ;;
    --python)      PY="$2"; shift 2 ;;
    --check)       CHECK_ONLY=1; shift ;;
    --no-selftest) SELFTEST=0; shift ;;
    -h|--help) sed -n '2,20p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0 ;;
    *) echo "unknown argument: $1" >&2; exit 2 ;;
  esac
done

ok()  { printf '  [ok]   %s\n' "$*"; }
bad() { printf '  [FAIL] %s\n' "$*"; FAILED=1; }
say() { printf '  %s\n' "$*"; }
FAILED=0

echo
echo "================================================================"
echo "  Distillation model install (ml/)"
echo "================================================================"

[[ -f ml/distill.py ]] || { bad "ml/distill.py missing — is this the repo root?"; exit 1; }
ok "ml/ found at $REPO/ml"

# ---- python ------------------------------------------------------------------
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

# ---- dependencies ------------------------------------------------------------
if (( ! CHECK_ONLY )); then
  say "installing from install/requirements-distill.txt (torch is ~250 MB) ..."
  ARGS=(--quiet --disable-pip-version-check -r install/requirements-distill.txt)
  if [[ -n "$VENV" ]]; then
    "$PY" -m pip install "${ARGS[@]}"
  else
    "$PY" -m pip install "${ARGS[@]}" --break-system-packages 2>/dev/null \
      || "$PY" -m pip install "${ARGS[@]}"
  fi
  [[ $? -eq 0 ]] && ok "dependencies installed" || bad "pip install failed"
fi

MISSING=$("$PY" - <<'EOF'
import importlib.util as u
need = ["torch", "numpy", "pandas", "rapidfuzz"]
print(" ".join(m for m in need if u.find_spec(m) is None))
EOF
)
if [[ -z "$MISSING" ]]; then
  ok "all 4 modules importable"
  "$PY" - <<'EOF'
import torch, importlib.metadata as md
for p in ("torch", "numpy", "pandas", "rapidfuzz"):
    try: print(f"         {p:12s} {md.version(p)}")
    except Exception: pass
print(f"         CUDA available: {torch.cuda.is_available()} "
      f"(CPU is sufficient at this scale)")
print(f"         threads: {torch.get_num_threads()}")
EOF
else
  bad "missing modules: $MISSING"
  say "torch wheels: https://pytorch.org/get-started/locally/ if the default fails"
fi

# ---- ml/metrics.py must agree with sptce_metrics.R ------------------------
# A learned model scored by a different equation cannot go in the same table as
# SPACE-TCR, so this is not optional.
if [[ -z "$MISSING" && -f ml/tests/test_metrics_match_r.py ]]; then
  if command -v Rscript >/dev/null 2>&1 && \
     Rscript -e 'quit(status = !requireNamespace("jsonlite", quietly = TRUE))' 2>/dev/null; then
    say "checking ml/metrics.py against sptce_metrics.R ..."
    if OUT=$("$PY" ml/tests/test_metrics_match_r.py 2>&1); then
      ok "agrees with sptce_metrics.R: $(grep -oE '[0-9]+ agreed, [0-9]+ disagreed' <<<"$OUT" | tail -1)"
    else
      bad "ml/metrics.py disagrees with sptce_metrics.R — results would not be comparable to Table 1/2"
      grep FAIL <<<"$OUT" | head -5 | sed 's/^/         /'
    fi
  else
    say "[skip] the R comparison needs R with jsonlite (Rscript install/r_packages.R --figures)"
  fi
fi

# ---- self-test: does it actually learn? --------------------------------------
if (( SELFTEST )) && [[ -z "$MISSING" ]]; then
  say "running the synthetic self-test (trains a model, ~1 min) ..."
  TMP=$(mktemp -d)
  if OUT=$("$PY" ml/distill.py --synthetic --epochs 25 --quiet --out "$TMP" 2>&1); then
    M=$(grep -A2 'AUROC' <<<"$OUT" | grep 'distilled model' | head -1 | awk '{print $3}')
    L=$(grep 'Levenshtein CDR3b' <<<"$OUT" | head -1 | awk '{print $3}')
    if [[ -n "$M" && -n "$L" ]] && "$PY" -c "import sys; sys.exit(0 if float('$M') > float('$L') else 1)"; then
      ok "self-test: model AUROC $M beats Levenshtein $L"
    else
      bad "self-test: model AUROC $M did NOT beat Levenshtein $L"
      say "       on the synthetic fixture the signal is in CDR3-alpha and V/J,"
      say "       where edit distance is blind by construction — if the model"
      say "       loses here, the install or the build is wrong, not the data."
    fi
  else
    bad "self-test failed to run"
    tail -5 <<<"$OUT" | sed 's/^/         /'
  fi
  rm -rf "$TMP"
fi

echo
echo "================================================================"
if (( FAILED )); then
  echo "  Install INCOMPLETE — see the failures above."
else
  echo "  Distillation model ready."
  echo
  echo "  Run on your own data (needs only the small files):"
  echo "    $PY ml/distill.py \\"
  echo "        --assignments cluster_assignments.tsv \\"
  echo "        --labels jouannet_tcrbuilder_input_exact_labelled.tsv \\"
  echo "        --out ml_out"
  echo
  echo "  Read ml/README.md first — especially what counts as a real result."
fi
echo "================================================================"
echo
exit $FAILED
