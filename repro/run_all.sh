#!/bin/sh
# End-to-end robustness pipeline:
#   1. build the boundary scan harness
#   2. run the full boundary scan (1592 cases, one process per case)
#   3. compare against the reference implementation's recorded scores
#
# Artifacts written:
#   artifacts/boundary_scan.csv        per-case scores and sanity flags
#   artifacts/reference_diff.csv       candidate vs reference deltas
#   artifacts/robustness_summary.json  pass/fail summary
set -e
cd "$(dirname "$0")/.."

CXX="${CXX:-g++}"
"$CXX" -std=c++11 -O2 -I. repro/boundary_scan.cc butteraugli/butteraugli.cc \
    -o repro/boundary_scan

python3 repro/run_scan.py --binary repro/boundary_scan \
    --out artifacts/boundary_scan.csv --summary /tmp/scan_summary.json
python3 repro/compare_reference.py
