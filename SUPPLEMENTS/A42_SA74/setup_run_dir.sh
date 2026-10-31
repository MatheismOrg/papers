#!/bin/bash
# Make a folder where the SA74 scripts find the module names they import, then run them from there.
set -e
S="$(cd "$(dirname "$0")/.." && pwd)"   # the SUPPLEMENTS folder
mkdir -p run && cd run
cp "$S/A42_Rule_Search.py" a42_rule_search.py
cp "$S/A42_SA74/A42_SA74_Classify.py" a42_classify.py
cp "$S/A42_SA74/A42_SA74_Extensions.py" a42x.py
cp "$S"/A42_SA74/*.py "$S"/A42_SA74/*.json "$S"/A42_SA74/*.jsonl "$S"/A42_SA74/*.c . 2>/dev/null || true
cp "$S/A38_Symbolic_Benchmark_Spec.txt" "$S/A42_Search_Results.json" "$S/A42_All_Rules.csv" .
echo "Ready in $(pwd). Example: A42_RUNS=\$(pwd) python3 A42_SA74_Numbers.py"
