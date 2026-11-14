#!/usr/bin/env python3
"""Compares a boundary scan against the reference implementation's scores.

Reads artifacts/boundary_scan.csv (candidate) and
repro/reference_scores.csv (recorded output of the reference
implementation), then writes:
  artifacts/reference_diff.csv       per-case candidate/reference deltas
  artifacts/robustness_summary.json  scan + reference summary

The pass criterion is that every case differs from the reference by less
than one butteraugli quantization level (kButteraugliQuantLow = 0.26).
"""

import csv
import json
import sys

QUANTIZATION_LEVEL = 0.26  # kButteraugliQuantLow

SCAN_PATH = "artifacts/boundary_scan.csv"
REFERENCE_PATH = "repro/reference_scores.csv"
DIFF_PATH = "artifacts/reference_diff.csv"
SUMMARY_PATH = "artifacts/robustness_summary.json"


def load_rows(path):
    with open(path) as f:
        reader = csv.reader(f)
        header = next(reader)
        return header, [row for row in reader]


def main():
    header, scan_rows = load_rows(SCAN_PATH)
    _, ref_rows = load_rows(REFERENCE_PATH)
    reference = {tuple(r[:5]): r[5] for r in ref_rows}
    candidate = {tuple(r[:5]): r for r in scan_rows}

    diff_header = ["case", "width", "height", "pattern_a", "pattern_b",
                   "candidate_score", "reference_score", "score_abs_delta",
                   "max_diffmap_abs_delta"]
    diff_rows = []
    score_abs_delta_max = 0.0
    diffmap_peak_abs_delta_max = 0.0
    missing = []
    for key in sorted(reference):
        row = candidate.get(key)
        if row is None or row[5] != "1":
            missing.append(key)
            continue
        candidate_score = float(row[6])
        candidate_peak = float(row[11])
        reference_score = float(reference[key])
        # The reference score is the peak of its diffmap, so the diffmap
        # peak delta compares candidate peak against the reference score.
        score_delta = abs(candidate_score - reference_score)
        peak_delta = abs(candidate_peak - reference_score)
        score_abs_delta_max = max(score_abs_delta_max, score_delta)
        diffmap_peak_abs_delta_max = max(diffmap_peak_abs_delta_max,
                                         peak_delta)
        diff_rows.append(list(key) + [
            "%.9g" % candidate_score,
            reference[key],
            repr(score_delta),
            repr(peak_delta),
        ])

    with open(DIFF_PATH, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(diff_header)
        writer.writerows(diff_rows)

    valid_rows = [r for r in scan_rows if r[5] == "1"]
    finite_rows = [r for r in valid_rows if r[9] == "1"]
    invalid_scores = [r for r in valid_rows if r[9] != "1" or r[10] != "1"]
    identical_nonzero = [r for r in valid_rows
                         if r[3] == r[4] and float(r[6]) != 0.0]
    scores = [float(r[6]) for r in valid_rows]
    ultrawide = [r for r in valid_rows if r[0].startswith("ultrawide")]
    ultrawide_nonfinite = [r for r in ultrawide if r[9] != "1"]
    zero_size = next(r for r in scan_rows if r[0] == "zero_size")

    passed = (not missing and not invalid_scores and
              not identical_nonzero and not ultrawide_nonfinite and
              score_abs_delta_max < QUANTIZATION_LEVEL and
              diffmap_peak_abs_delta_max < QUANTIZATION_LEVEL)
    summary = {
        "scan": {
            "rows": len(scan_rows),
            "valid_finite_rows": len(finite_rows),
            "invalid_scores": len(invalid_scores),
            "identical_images_nonzero_score": len(identical_nonzero),
            "zero_size_case": dict(zip(header, zero_size)),
            "score_min": min(scores) if scores else None,
            "score_max": max(scores) if scores else None,
            "ultrawide_cases": len(ultrawide),
            "ultrawide_nonfinite": len(ultrawide_nonfinite),
        },
        "reference": {
            "cases": len(diff_rows),
            "missing_cases": [list(k) for k in missing],
            "score_abs_delta_max": score_abs_delta_max,
            "diffmap_peak_abs_delta_max": diffmap_peak_abs_delta_max,
            "quantization_level": QUANTIZATION_LEVEL,
            "pass": passed,
        },
    }
    with open(SUMMARY_PATH, "w") as f:
        json.dump(summary, f, indent=2)
        f.write("\n")

    print("reference cases=%d missing=%d score_delta_max=%.3e "
          "peak_delta_max=%.3e pass=%s" % (
              len(diff_rows), len(missing), score_abs_delta_max,
              diffmap_peak_abs_delta_max, passed))
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())
