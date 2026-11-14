# Differential regression harness

This harness compares the working tree against a reference revision in the
same process. The reference source is extracted from git and its namespace is
renamed at build time; both implementations receive byte-for-byte identical
linear RGB planes.

Run the default campaign (5000 groups, five distortion severities):

```sh
tests/run_differential_regression.sh
```

Useful overrides:

```sh
REFERENCE_REV=<git-rev> tests/run_differential_regression.sh
tests/run_differential_regression.sh --groups 5000 --seed 12345
tests/run_differential_regression.sh --min-spearman 0.95 --rank-tolerance 0.02
```

The harness checks:

- Spearman rank correlation against the reference implementation.
- Non-decreasing Butteraugli scores as the distortion blend increases.
- Bitwise-identical candidate scores on two evaluations of the same inputs.
- Candidate crashes on degenerate and extreme image dimensions, isolated in
  child processes.
- Per-sample rank outliers (`|rank_candidate - rank_reference|` greater than
  `--rank-tolerance` times the cohort size).

Every monotonicity or rank-outlier sample is minimized by greedy crop/severity
reduction. A predicate is accepted only when the failure remains reproducible
on the reduced input. Reports are written to
`artifacts/differential-regression/`:

- `summary.json`, `summary.md`
- `spearman_by_level.csv`, `spearman_by_family.csv`
- `counterexamples.csv`, `counterexamples.md`
- `extreme_sizes.csv`
- cropped PPM pairs for the first `--max-artifacts` counterexamples
