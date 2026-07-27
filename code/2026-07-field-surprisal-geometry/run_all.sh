#!/usr/bin/env bash
# Run the entire verification suite in dependency order and report pass counts.
# Every harness is fail-first: it exits non-zero on the first failed check.
# Ordering note: t1_engine.py writes P_coeffs.pkl, which t1_windowproof.py reads
# for its final (462-coefficient dictionary) cross-check. Run engine first.
set -u
cd "$(dirname "$0")/harnesses"
pass=0; fail=0; total_checks=0

run () {  # run <file> <timeout-seconds>
  local f="$1" t="${2:-600}"
  printf '%-34s ' "$f"
  if out=$(timeout "$t" python3 "$f" 2>&1); then
    n=$(printf '%s' "$out" | grep -oE 'ALL [0-9]+ CHECKS PASSED|[0-9]+/[0-9]+ passed' | tail -1)
    echo "OK   $n"; pass=$((pass+1))
  else
    code=$?
    if [ "$code" -eq 124 ]; then echo "TIMEOUT (${t}s) — see note in MANIFEST.md"; else echo "FAIL (exit $code)"; fi
    fail=$((fail+1))
  fi
}

echo "== foundational =="
run field_surprisal_v2.py
run suspension_theorem.py
run field_surprisal_tier2.py
run field_surprisal_classification.py
echo "== core classification (k=2) =="
run t1_core.py
run t1_reduction.py
run t1_engine.py            # writes P_coeffs.pkl
run t1_branches.py
run t1_windowproof.py       # reads P_coeffs.pkl -> full 44/44
run t3c_partC_exact.py
echo "== thermodynamics / suspension =="
run t2_temperature.py
run t3_suspension.py
echo "== higher-k, census, selection =="
run t4_kwindows.py 300      # lane-1 checks fast; block-3 symbolic census is slow (see note)
run t4b_census_fast.py 900  # practical census lane
run t5_catalog_census.py 900
run t6_selection.py
echo "== last two fronts =="
run t7_knecessity.py
run t8_compositum.py
run t9_landscape.py 900
run t10_coupled.py

echo
echo "harnesses OK: $pass   not-clean: $fail"
echo "(full-suite expected total: 422 checks across 20 harnesses, all exit 0;"
echo " t4_kwindows block-3 census may TIMEOUT in a small sandbox — it is the same"
echo " result the fast lane t4b_census_fast proves, with sampled two-lane agreement.)"
