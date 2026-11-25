#!/bin/bash
# =============================================================================
# check_install.sh — verify a SPTCE installation
#
#   bash check_install.sh            # check whatever is installed
#   bash check_install.sh --analysis # R analysis layer only (laptop)
#   bash check_install.sh --pipeline # cluster pipeline only (HPC)
#
# Exits 0 if everything checked passed, 1 otherwise. Needs no data: every test
# builds its own fixtures. See INSTALL.md.
# =============================================================================
set -uo pipefail
cd "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

MODE="${1:---all}"
PASS=0; FAIL=0; SKIP=0
ok()   { echo "  [ok]   $*"; PASS=$((PASS+1)); }
bad()  { echo "  [FAIL] $*"; FAIL=$((FAIL+1)); }
skip() { echo "  [skip] $*"; SKIP=$((SKIP+1)); }
hdr()  { echo; echo "=== $* ==="; }

want_analysis=0; want_pipeline=0
R_OK=0
case "$MODE" in
  --analysis) want_analysis=1 ;;
  --pipeline) want_pipeline=1 ;;
  *)          want_analysis=1; want_pipeline=1 ;;
esac

# -----------------------------------------------------------------------------
hdr "0. repository layout"
for f in sptce_metrics.R sptce_pipeline.R config.sh \
         SPACE2TCR/util.py scripts/make_bins.py sbatch/run_C_slurm.sh \
         tests/test_metrics.R tests/test_e2e.R tests/verify_space2tcr.py; do
  [[ -e "$f" ]] && ok "$f" || bad "$f is missing"
done
if [[ -e SPACE2TCR/_upstream/LICENSE ]]; then
  ok "SPACE2TCR/_upstream/LICENSE present"
else
  skip "SPACE2TCR/_upstream/LICENSE absent — required before publishing (see _upstream/NOTICE.md)"
fi

# -----------------------------------------------------------------------------
if (( want_analysis )); then
hdr "1. R toolchain"
if command -v Rscript >/dev/null 2>&1; then
  ok "Rscript: $(Rscript -e 'cat(R.version.string)' 2>/dev/null)"
  MISSING_R="$(Rscript -e 'p<-c("readr","dplyr","tidyr","purrr","stringr","tibble","stringdist");m<-p[!sapply(p,requireNamespace,quietly=TRUE)];cat(paste(m,collapse=" "))' 2>/dev/null)"
  if [[ -z "$MISSING_R" ]]; then
    R_OK=1
    ok "all 7 required R packages present"
  else
    bad "missing R packages: $MISSING_R"
    echo "         install.packages(c($(echo "$MISSING_R" | sed 's/\([^ ]*\)/"\1"/g; s/ /,/g')))"
  fi
  # Optional, only needed for the figure code in sptce_pipeline.R.
  OPT_R="$(Rscript -e 'p<-c("ggplot2","patchwork","scales","igraph","ggraph","tidygraph","cowplot","ggseqlogo","ggforce","ggrepel","reshape2","forcats");m<-p[!sapply(p,requireNamespace,quietly=TRUE)];cat(paste(m,collapse=" "))' 2>/dev/null)"
  [[ -z "$OPT_R" ]] && ok "optional figure packages present (figures will render)" \
                    || skip "figure packages absent: $OPT_R  (validation and the tests do not need them)"
else
  bad "Rscript not found — install R (see INSTALL.md)"
fi

hdr "2. R metrics and pipeline"
if [[ "${R_OK:-0}" == "1" ]]; then
  if Rscript -e 'invisible(parse("sptce_metrics.R"))' >/dev/null 2>&1; then
    ok "sptce_metrics.R parses"
  else bad "sptce_metrics.R does not parse"; fi
  for f in sptce_pipeline.R; do
    [[ -f "$f" ]] || continue
    Rscript -e "invisible(parse('$f'))" >/dev/null 2>&1 && ok "$f parses" || bad "$f does not parse"
  done
  # The pipeline must be complete: the stage functions it needs, plus the
  # dihedral gate setting, which must actually be read.
  if out="$(Rscript -e '
      src <- readLines("sptce_pipeline.R", warn=FALSE)
      code <- src[!grepl("^\\s*#", src)]
      d <- gsub("`","", trimws(sub("\\s*(<-|=)\\s*(function|\\\\\\().*","",
           grep("^\\s*`?\\.?[A-Za-z._][A-Za-z0-9._]*`?\\s*(<-|=)\\s*(function|\\\\\\()", code, value=TRUE))))
      need <- c("run_benchmark_a","jouannet_10x_ladder",".size_col",".sptce_tag",
                ".focal_universe",".assign_overlay_factor","run_stage5",
                "run_stage5d_cohort_accounting","main")
      miss <- setdiff(need, d)
      if (length(miss)) { cat("MISSING:", paste(miss, collapse=" "), "\n"); quit(status=1) }
      n <- length(grep("CFG\\$focal_dihedral", src, fixed=FALSE))
      if (n < 2) { cat("focal_dihedral read only", n, "time(s) - the knob is dead\n"); quit(status=1) }
      cat(length(d), "definitions; focal_dihedral read", n, "times\n")' 2>&1)"; then
    ok "sptce_pipeline.R is complete: $out"
  else
    bad "sptce_pipeline.R incomplete: $out"
  fi

  echo "  running tests/test_metrics.R ..."
  if out="$(Rscript tests/test_metrics.R 2>&1)"; then
    ok "metric unit tests: $(echo "$out" | grep -oE '[0-9]+ passed, [0-9]+ failed' | tail -1)"
  else
    bad "metric unit tests failed"; echo "$out" | grep -E 'FAIL|Error' | head -10 | sed 's/^/         /'
  fi

  echo "  running tests/test_pipeline.R ..."
  if out="$(Rscript tests/test_pipeline.R 2>&1)"; then
    ok "pipeline checks: $(echo "$out" | grep -oE '[0-9]+ passed, [0-9]+ failed' | tail -1)"
  else
    bad "pipeline checks failed"; echo "$out" | grep -E 'FAIL' | head -10 | sed 's/^/         /'
  fi

  echo "  running tests/test_e2e.R ..."
  if out="$(Rscript tests/test_e2e.R 2>&1)"; then
    ok "end-to-end tests: $(echo "$out" | grep -oE 'e2e: [0-9]+ failures' | tail -1)"
  else
    bad "end-to-end tests failed"; echo "$out" | grep -E 'FAIL|Error' | head -10 | sed 's/^/         /'
  fi
else
  skip "R tests (toolchain incomplete)"
fi
fi

# -----------------------------------------------------------------------------
if (( want_pipeline )); then
hdr "3. Python toolchain and SPACE2TCR"
PY="${PY_SPACE2:-}"
if [[ -z "$PY" && -f config.sh ]]; then
  # shellcheck disable=SC1091
  PY="$(bash -c 'source ./config.sh >/dev/null 2>&1; echo "${PY_SPACE2:-}"' 2>/dev/null)"
fi
[[ -x "${PY:-}" ]] || PY="$(command -v python3 || true)"

if [[ -n "$PY" ]]; then
  ok "python: $("$PY" -V 2>&1)  ($PY)"
  MISS_PY="$("$PY" - <<'EOF' 2>/dev/null
import importlib
need = ["numpy","pandas","scipy","numba","joblib","sklearn"]
print(" ".join(m for m in need if not importlib.util.find_spec(m)))
EOF
)"
  if [[ -z "$MISS_PY" ]]; then ok "SPACE2TCR dependencies present (numpy pandas scipy numba joblib sklearn)"
  else bad "missing python packages: $MISS_PY"; echo "         pip install $MISS_PY"; fi

  if PYTHONPATH=. "$PY" -c 'import SPACE2TCR.util' >/dev/null 2>&1; then
    ok "SPACE2TCR imports (PYTHONPATH=.)"
    echo "  running tests/verify_space2tcr.py ..."
    if out="$(PYTHONPATH=. "$PY" tests/verify_space2tcr.py . 2>&1)"; then
      ok "SPACE2TCR verification: $(echo "$out" | grep -c '  PASS  ') checks passed"
    else
      bad "SPACE2TCR verification failed"; echo "$out" | grep -E 'FAIL' | head -10 | sed 's/^/         /'
    fi
  else
    bad "cannot import SPACE2TCR — Stage C will not run"
  fi

  if [[ -f ml/metrics.py ]]; then
    if "$PY" -c 'import torch, rapidfuzz' >/dev/null 2>&1; then
      echo "  running ml/tests/test_metrics_match_r.py ..."
      if out="$("$PY" ml/tests/test_metrics_match_r.py 2>&1)"; then
        ok "ml/metrics.py agrees with sptce_metrics.R: $(echo "$out" | grep -oE '[0-9]+ agreed, [0-9]+ disagreed' | tail -1)"
      else
        bad "ml/metrics.py disagrees with sptce_metrics.R"; echo "$out" | grep FAIL | head -5 | sed 's/^/         /'
      fi
    else
      skip "ml metrics check (needs torch and rapidfuzz: pip install torch rapidfuzz)"
    fi
  fi

  for f in scripts/*.py; do
    "$PY" -m py_compile "$f" 2>/dev/null || bad "does not compile: $f"
  done
  ok "all $(ls scripts/*.py | wc -l) pipeline scripts compile"
else
  bad "no python3 found"
fi

hdr "4. shell scripts and config"
for f in sbatch/*.sh collect_from_cluster.sh check_install.sh; do
  bash -n "$f" 2>/dev/null || bad "syntax error: $f"
done
ok "all shell scripts parse"
if [[ -f config.sh ]]; then
  MISSING_VARS=""
  for v in INPUT_CSV WORK_DIR SPECIES CUTOFF SPTCE_ROOT PY_TCRBUILDER PY_SPACE2 \
           STITCHR_ENV DATA PY_CLUSTCR; do
    grep -qE "^[[:space:]]*${v}=" config.sh || MISSING_VARS="$MISSING_VARS $v"
  done
  [[ -z "$MISSING_VARS" ]] && ok "config.sh defines all 10 core variables" \
                           || bad "config.sh is missing:$MISSING_VARS"
else
  bad "config.sh not found"
fi

hdr "5. cluster-only prerequisites (informational)"
command -v sbatch  >/dev/null 2>&1 && ok "sbatch present"  || skip "sbatch absent — not on a Slurm cluster"
command -v hmmscan >/dev/null 2>&1 && ok "hmmscan present" || skip "hmmscan absent — Stage B needs HMMER"
"${PY}" -c 'import ImmuneBuilder' >/dev/null 2>&1 && ok "ImmuneBuilder present" \
  || skip "ImmuneBuilder absent — Stage B only, lives in the myenv2TCR env"
"${PY}" -c 'import clustcr' >/dev/null 2>&1 && ok "clustcr present" \
  || skip "clustcr absent — Stage E only, lives in the clustcr_env env"
fi

# -----------------------------------------------------------------------------
echo
echo "================================================================"
printf "  %d passed, %d failed, %d skipped\n" "$PASS" "$FAIL" "$SKIP"
if (( FAIL == 0 )); then
  echo "  Install looks good."
  (( SKIP > 0 )) && echo "  Skips are optional components — see INSTALL.md if you need them."
else
  echo "  Fix the failures above. INSTALL.md has the install steps."
fi
echo "================================================================"
(( FAIL == 0 ))
