#!/usr/bin/env python3
"""Boundary scan driver for the butteraugli robustness harness.

Runs repro/boundary_scan once per case (subprocess isolation: a crash on
one case is recorded as CRASH and does not abort the scan) and writes
artifacts/boundary_scan.csv plus artifacts/robustness_summary.json.

Usage:
  python3 repro/run_scan.py [--binary repro/boundary_scan] \
      [--out artifacts/boundary_scan.csv] \
      [--summary artifacts/robustness_summary.json]
"""

import argparse
import csv
import json
import math
import subprocess
import sys

SIZES = [1, 2, 3, 5, 7, 8, 9, 15, 16, 17, 31, 32, 33, 64]

# (case name, pattern_a, pattern_b)
PATTERN_CASES = [
    ("black_black", 0, 0),
    ("white_white", 1, 1),
    ("black_white", 0, 1),
    ("white_black", 1, 0),
    ("black_near_white", 0, 5),
    ("checker_black", 2, 0),
    ("hot_pixel_black", 3, 0),
    ("ramp_black", 4, 0),
]

ULTRAWIDE_SIZES = [
    (1, 4096), (4096, 1), (2, 2048), (2048, 2), (1, 65536),
    (65536, 1), (3, 4096), (4096, 3), (16, 4096), (4096, 16),
]

SPECIAL_CASES = [
    ("fltmax_black", 16, 16, 6, 0),
    ("posinf_black", 16, 16, 7, 0),
    ("nan_black", 16, 16, 8, 0),
]


def run_case(binary, width, height, pattern_a, pattern_b):
    """Returns (status, fields) where status is 'ok', 'crash' or 'error'."""
    cmd = [binary,
           "--width=%d" % width, "--height=%d" % height,
           "--pattern_a=%d" % pattern_a, "--pattern_b=%d" % pattern_b]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        return ("crash", None)
    fields = proc.stdout.strip().split(",")
    if len(fields) != 7:
        return ("error", None)
    return ("ok", fields)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", default="repro/boundary_scan")
    parser.add_argument("--out", default="artifacts/boundary_scan.csv")
    parser.add_argument("--summary", default="artifacts/robustness_summary.json")
    args = parser.parse_args()

    cases = []
    for width in SIZES:
        for height in SIZES:
            for name, pattern_a, pattern_b in PATTERN_CASES:
                cases.append((name, width, height, pattern_a, pattern_b))
    for width, height in ULTRAWIDE_SIZES:
        cases.append(("ultrawide_black_white", width, height, 0, 1))
        cases.append(("ultrawide_black_black", width, height, 0, 0))
    cases.extend(SPECIAL_CASES)
    cases.append(("zero_size", 0, 0, 0, 0))

    header = ["case", "width", "height", "pattern_a", "pattern_b", "ok",
              "score", "diffmap_width", "diffmap_height", "finite",
              "nonnegative", "max_diffmap"]
    rows = []
    crashes = []
    for name, width, height, pattern_a, pattern_b in cases:
        status, fields = run_case(args.binary, width, height, pattern_a,
                                  pattern_b)
        if status != "ok":
            crashes.append((name, width, height, pattern_a, pattern_b))
            rows.append([name, width, height, pattern_a, pattern_b,
                         "CRASH", "", "", "", "", "", ""])
            continue
        rows.append([name, width, height, pattern_a, pattern_b] + fields)

    with open(args.out, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerows(rows)

    # Summarize.
    valid_rows = [r for r in rows if r[5] == "1"]
    finite_rows = [r for r in valid_rows if r[9] == "1"]
    invalid_scores = [r for r in valid_rows
                      if r[9] != "1" or r[10] != "1"]
    identical_nonzero = [r for r in valid_rows
                         if r[3] == r[4] and float(r[6]) != 0.0]
    scores = [float(r[6]) for r in valid_rows]
    ultrawide = [r for r in valid_rows if r[0].startswith("ultrawide")]
    ultrawide_nonfinite = [r for r in ultrawide if r[9] != "1"]
    zero = rows[-1]
    summary = {
        "scan": {
            "rows": len(rows),
            "valid_finite_rows": len(finite_rows),
            "crashes": len(crashes),
            "crash_cases": [list(c) for c in crashes],
            "invalid_scores": len(invalid_scores),
            "identical_images_nonzero_score": len(identical_nonzero),
            "zero_size_case": dict(zip(header, zero)),
            "score_min": min(scores) if scores else None,
            "score_max": max(scores) if scores else None,
            "ultrawide_cases": len(ultrawide),
            "ultrawide_nonfinite": len(ultrawide_nonfinite),
        },
    }
    with open(args.summary, "w") as f:
        json.dump(summary, f, indent=2)
        f.write("\n")

    print("rows=%d crashes=%d invalid=%d identical_nonzero=%d" % (
        len(rows), len(crashes), len(invalid_scores),
        len(identical_nonzero)))
    if crashes:
        print("crash cases (first 10):")
        for c in crashes[:10]:
            print("  ", c)
    return 1 if crashes or invalid_scores else 0


if __name__ == "__main__":
    sys.exit(main())
