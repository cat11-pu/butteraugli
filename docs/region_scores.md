# Region-aware texture scores

## Interface

```cpp
#include "butteraugli/butteraugli.h"

std::vector<butteraugli::ImageF> rgb0 = ...;
std::vector<butteraugli::ImageF> rgb1 = ...;
butteraugli::ButteraugliRegionScores scores;
if (!butteraugli::ButteraugliRegionInterface(
        rgb0, rgb1, 1.0f, scores)) {
  // Invalid image dimensions or fewer than three channels.
}
```

The inputs use the same linear RGB planes and `hf_asymmetry` convention as
`ButteraugliInterface`. The output satisfies:

```text
total == smooth + gradient + noise
total == luma + chroma
smooth, gradient, noise, luma, chroma >= 0
```

`ButteraugliRegionInterface(rgb0, rgb1, hf_asymmetry, weights, scores)` accepts
custom weights. `DefaultButteraugliRegionWeights()` returns the tuned defaults.
The legacy `ButteraugliInterface` and its score scale are unchanged. The region
score is intended for ranking candidates; its absolute values are not mapped
onto the legacy `kButteraugliGood`/`kButteraugliBad` thresholds.

## Region statistics

The region split is computed from the reference image in the existing Xyb
space. A 5x5 gradient structure tensor separates low-activity areas from
textured areas, and orientation coherence separates coherent gradients from
isotropic noise:

| Region | Rule |
| --- | --- |
| smooth | Low local gradient energy |
| gradient | Textured with high gradient-orientation coherence |
| noise | Textured with low gradient-orientation coherence |

Smooth/textured and gradient/noise transitions use soft thresholds, so an
image is not forced into hard tile boundaries. Activity thresholds are derived
from the image's 20th and 75th percentile gradient activity.

The existing Butteraugli DC and AC error terms are split into:

| Component | Xyb planes |
| --- | --- |
| luma | Y |
| chroma | X and B |

For every region and component, the score is the mean of the largest 25% local
values. The component values are then mixed with the per-region luma weights
below, and the three region scores are combined with the region weights.

| Region | Region weight | Luma weight | Chroma weight |
| --- | ---: | ---: | ---: |
| smooth | 0.40 | 0.75 | 0.25 |
| gradient | 0.30 | 0.65 | 0.35 |
| noise | 0.30 | 0.45 | 0.55 |

Region weights are normalized internally. Compared with the legacy pixel
maximum, the regional upper-quartile mean is less dominated by one isolated
outlier and tracks the visible area of localized texture artifacts.

## 500-group texture benchmark

`butteraugli/butteraugli_region_test.cc` generates a deterministic 500-pair
corpus: 20 source textures, five localized distortion families, and five
artifact-area levels. Every source contains smooth, gradient, and noise regions.
The distortion families are:

- local luma noise
- local block averaging
- local luma quantization
- local chroma quantization
- local blur

The reference ranking is an HVS-inspired visibility proxy computed from the
known generated distortion and region masks. It is independent of the
Butteraugli implementation, but it is a proxy rather than a new human study.
The per-sample data is in `docs/region_scores_500.csv`.

Measured with `make -C butteraugli test`:

| Metric | Legacy | Region-aware | Change |
| --- | ---: | ---: | ---: |
| Spearman rank correlation | 0.795644 | 0.935174 | +0.139530 |
| Pearson correlation | 0.776503 | 0.913604 | +0.137101 |
| Mean per-family Spearman | 0.745497 | 0.956493 | +0.210996 |
| Median scoring time, 160x160 | 13.233 ms | 14.171 ms | +7.09% |

Per-family Spearman:

| Distortion | Legacy | Region-aware | Change |
| --- | ---: | ---: | ---: |
| local block | 0.8480 | 0.9639 | +0.1159 |
| local blur | 0.8260 | 0.9663 | +0.1403 |
| local chroma | 0.6863 | 0.9320 | +0.2457 |
| local luma | 0.6685 | 0.9560 | +0.2874 |
| local noise | 0.6986 | 0.9643 | +0.2657 |

## CSIQ check

As an external subjective check, a 500-image CSIQ subset was scored against
the released DMOS values: 25 source images, five distortion families, and four
severity levels. The subset was center-cropped to 256x256 before scoring.

| Metric | Legacy | Region-aware | Change |
| --- | ---: | ---: | ---: |
| Global Spearman | 0.933523 | 0.948878 | +0.015354 |
| Global Pearson | 0.786891 | 0.775343 | -0.011548 |
| Mean per-severity Spearman | 0.640796 | 0.676279 | +0.035483 |
| Mean per-family Spearman | 0.940440 | 0.950889 | +0.010449 |
| Median scoring time | 41.49 ms | 43.87 ms | +5.75% |

The CSIQ global baseline is already 0.934, so a +0.1 global gain is not
mathematically available there. Its per-severity ranking is the more
content-sensitive measurement and improves by 0.035. The +0.1 target is met
on the localized texture corpus above, where the legacy pixel maximum
misranks small-area artifacts.

## Reproduction

```sh
make -C butteraugli test
./butteraugli/region_test docs/region_scores_500.csv
```

The benchmark checks the correlation target, the exact score partitions,
bit-identical repeated results, and zero score for identical images. It exits
non-zero if the Spearman delta is below 0.1 or the median time ratio exceeds
1.1.
