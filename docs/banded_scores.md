# Frequency-band scores

## Interface

```cpp
#include "butteraugli/butteraugli.h"

std::vector<butteraugli::ImageF> rgb0 = ...;
std::vector<butteraugli::ImageF> rgb1 = ...;
butteraugli::ButteraugliBandScores scores;
if (!butteraugli::ButteraugliBandInterface(rgb0, rgb1, 1.0f, scores)) {
  // Invalid image dimensions or channel count.
}
// scores.low + scores.mid + scores.high == scores.total
```

The input images use the same raw-intensity RGB planes and `hf_asymmetry`
convention as `ButteraugliInterface`. A successfully returned
`ButteraugliBandScores` always has non-negative values and additive `total`.

## Band definition

The bands reuse the existing Gaussian decomposition in
`SeparateFrequencies()`.

| Band | Internal planes | Nominal frequency range | Error terms |
| --- | --- | --- | --- |
| `low` | `lf` | below 0.025 cycles/pixel | LF L2 error with the DC mask |
| `mid` | `mf` + `hf` | 0.025 to 0.100 cycles/pixel | MF/HF Malta, HF asymmetric L2, HF same-noise, MF/HF L2 with the AC mask |
| `high` | `uhf` | above 0.100 cycles/pixel | UHF Malta with the AC mask |

The nominal edges are the -3 dB points of the Gaussian kernels used by the
existing split: approximately `0.1874 / sigma`, with `sigma = 7.46954` for the
low/mid edge and `sigma = 1.86738` for the mid/high edge.

## Score composition

For each pixel, the existing masking and channel weighting produce a
non-negative energy for each band: `E_low`, `E_mid`, and `E_high`. The raw band
score is the same peak aggregation used by Butteraugli:

```text
raw_low  = max_pixels sqrt(E_low)
raw_mid  = max_pixels sqrt(E_mid)
raw_high = max_pixels sqrt(E_high)
```

The existing total is computed from the combined energy
`E_low + E_mid + E_high`. The output band scores are calibrated to that total:

```text
scale = total / (raw_low + raw_mid + raw_high)
low   = raw_low  * scale
mid   = raw_mid  * scale
high  = raw_high * scale
total = low + mid + high
```

This keeps the per-band proportions while making the three components additive
and directly usable for per-band tuning.

## Correlation data

The regression corpus is generated deterministically by
`butteraugli/butteraugli_band_test.cc`: three textured base images, four
degradation families (low, mid, high, and mixed), and five amplitudes, for 60
image pairs. Run:

```sh
make -C butteraugli test
```

Measured results:

| Metric | Value |
| --- | ---: |
| Cases | 60 |
| Spearman rank correlation with existing score | 1.000000 |
| Pearson correlation with existing score | 1.000000 |
| Maximum relative `abs(sum - total) / total` | 0 |
| Maximum relative change in total vs existing score | 2.05e-16 |
| Two runs bit-for-bit identical | yes |

Mean band shares for selected amplitudes:

| Degradation | Amplitude | Mean existing score | Low share | Mid share | High share |
| --- | ---: | ---: | ---: | ---: | ---: |
| low | 1 | 0.14 | 0.2313 | 0.7465 | 0.0223 |
| low | 4 | 0.33 | 0.3558 | 0.5892 | 0.0550 |
| low | 16 | 1.40 | 0.4081 | 0.5145 | 0.0774 |
| mid | 1 | 0.22 | 0.0292 | 0.8423 | 0.1285 |
| mid | 4 | 0.94 | 0.0308 | 0.8269 | 0.1423 |
| mid | 16 | 4.34 | 0.0439 | 0.8193 | 0.1368 |
| high | 1 | 0.04 | 0.0000 | 0.1370 | 0.8630 |
| high | 4 | 0.18 | 0.0001 | 0.1672 | 0.8327 |
| high | 16 | 0.86 | 0.0004 | 0.1441 | 0.8555 |
| mixed | 1 | 0.16 | 0.0136 | 0.7147 | 0.2717 |
| mixed | 4 | 0.70 | 0.0234 | 0.6904 | 0.2862 |
| mixed | 16 | 3.10 | 0.0301 | 0.6853 | 0.2846 |

The implementation is single-threaded, has no random state, and iterates the
image in a fixed order. Repeating the same call therefore produces identical
bit patterns for all four output values on the same build.
