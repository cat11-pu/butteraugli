#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD_DIR="${BUILD_DIR:-${ROOT}/build/differential-regression}"
REFERENCE_REV="${REFERENCE_REV:-HEAD}"
CXX="${CXX:-c++}"

OUT_DIR="artifacts/differential-regression"
prev=""
for arg in "$@"; do
  if [[ "${prev}" == "--out" ]]; then
    OUT_DIR="${arg}"
  fi
  prev="${arg}"
done
if [[ "${OUT_DIR}" != /* ]]; then
  OUT_DIR="${ROOT}/${OUT_DIR}"
fi
mkdir -p "${OUT_DIR}"

REFERENCE_SHA="$(git -C "${ROOT}" rev-parse "${REFERENCE_REV}")"
CANDIDATE_SHA="$(
  cd "${ROOT}"
  sha256sum butteraugli/butteraugli.cc butteraugli/butteraugli.h |
    sha256sum | cut -d' ' -f1
)"

rm -rf "${BUILD_DIR}"
mkdir -p "${BUILD_DIR}/reference"
git -C "${ROOT}" archive "${REFERENCE_REV}" \
  butteraugli/butteraugli.cc butteraugli/butteraugli.h |
  tar -x -C "${BUILD_DIR}/reference"

# The reference is linked into the same test process. Rename its public
# namespace and include guard so both implementations can coexist.
sed -i \
  -e 's/BUTTERAUGLI_BUTTERAUGLI_H_/BUTTERAUGLI_REF_BUTTERAUGLI_H_/g' \
  -e 's/namespace butteraugli/namespace butteraugli_ref/g' \
  -e 's/butteraugli::/butteraugli_ref::/g' \
  "${BUILD_DIR}/reference/butteraugli/butteraugli.h" \
  "${BUILD_DIR}/reference/butteraugli/butteraugli.cc"

COMMON_FLAGS=(-std=c++11 -O2 -g -Wall -Wextra -Wno-sign-compare
              -Wno-ignored-qualifiers -Wno-unused-parameter
              -Wno-unused-function -fno-omit-frame-pointer)

echo "Building reference ${REFERENCE_REV} (${REFERENCE_SHA:0:12}) and candidate..."
"${CXX}" "${COMMON_FLAGS[@]}" -I"${BUILD_DIR}/reference" \
  -c "${BUILD_DIR}/reference/butteraugli/butteraugli.cc" \
  -o "${BUILD_DIR}/reference_butteraugli.o"
"${CXX}" "${COMMON_FLAGS[@]}" -I"${ROOT}" \
  -c "${ROOT}/butteraugli/butteraugli.cc" \
  -o "${BUILD_DIR}/candidate_butteraugli.o"
"${CXX}" "${COMMON_FLAGS[@]}" -I"${ROOT}" -I"${ROOT}/tests" \
  -c "${ROOT}/tests/candidate_scorer.cc" \
  -o "${BUILD_DIR}/candidate_scorer.o"
"${CXX}" "${COMMON_FLAGS[@]}" -I"${BUILD_DIR}/reference" \
  -I"${ROOT}/tests" \
  -c "${ROOT}/tests/reference_scorer.cc" \
  -o "${BUILD_DIR}/reference_scorer.o"
"${CXX}" "${COMMON_FLAGS[@]}" -I"${ROOT}/tests" \
  -c "${ROOT}/tests/differential_regression.cc" \
  -o "${BUILD_DIR}/differential_regression.o"
"${CXX}" "${COMMON_FLAGS[@]}" \
  "${BUILD_DIR}/differential_regression.o" \
  "${BUILD_DIR}/candidate_scorer.o" \
  "${BUILD_DIR}/candidate_butteraugli.o" \
  "${BUILD_DIR}/reference_scorer.o" \
  "${BUILD_DIR}/reference_butteraugli.o" \
  -o "${BUILD_DIR}/differential_regression"

echo "Reports: ${OUT_DIR}"
set +e
"${BUILD_DIR}/differential_regression" \
  --out "${OUT_DIR}" \
  --reference-label "git ${REFERENCE_REV}" \
  --candidate-label "working tree" \
  --reference-sha "${REFERENCE_SHA}" \
  --candidate-sha "${CANDIDATE_SHA}" \
  "$@" 2>&1 | tee "${OUT_DIR}/run.log"
STATUS=${PIPESTATUS[0]}
set -e
exit "${STATUS}"
