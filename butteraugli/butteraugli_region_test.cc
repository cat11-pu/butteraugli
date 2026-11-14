#include "butteraugli/butteraugli.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <numeric>
#include <string>
#include <vector>

using namespace butteraugli;

namespace {

constexpr int kW = 160;
constexpr int kH = 160;
constexpr int kBases = 20;
constexpr int kFamilies = 5;
constexpr int kLevels = 5;
constexpr double kPi = 3.14159265358979323846;

struct Rng {
  uint64_t state;

  explicit Rng(uint64_t seed) : state(seed) {}

  uint64_t Next() {
    uint64_t z = (state += 0x9e3779b97f4a7c15ULL);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
  }

  int Range(int n) { return n <= 0 ? 0 : static_cast<int>(Next() % n); }
  double Unit() { return (Next() >> 11) * (1.0 / 9007199254740992.0); }
  double Signed() { return 2.0 * Unit() - 1.0; }
  double Gaussian() {
    const double u1 = std::max(1e-12, Unit());
    return std::sqrt(-2.0 * std::log(u1)) *
           std::cos(2.0 * kPi * Unit());
  }
};

struct Rgb {
  int w = kW;
  int h = kH;
  std::vector<uint8_t> p;

  Rgb() : p(static_cast<size_t>(w) * h * 3, 0) {}

  uint8_t* At(int x, int y) {
    return &p[(static_cast<size_t>(y) * w + x) * 3];
  }
  const uint8_t* At(int x, int y) const {
    return &p[(static_cast<size_t>(y) * w + x) * 3];
  }
};

int Clamp255(double value) {
  return static_cast<int>(
      std::max(0.0, std::min(255.0, std::round(value))));
}

enum Region { kSmooth = 0, kGradient = 1, kNoise = 2 };

struct Base {
  Rgb rgb;
  std::vector<uint8_t> region;
};

Base MakeBase(uint64_t seed) {
  Rng rng(seed);
  static const int kinds[] = {kSmooth, kGradient, kNoise};
  const int num_bands = 4 + rng.Range(2);
  std::vector<int> bands(num_bands);
  for (int i = 0; i < num_bands; ++i) bands[i] = kinds[i % 3];
  for (int i = num_bands - 3; i < num_bands; ++i) {
    bands[i] = kinds[rng.Range(3)];
  }
  for (int i = num_bands - 1; i > 0; --i) {
    std::swap(bands[i], bands[rng.Range(i + 1)]);
  }
  bands[0] = kSmooth;
  bands[1] = kGradient;
  bands[2] = kNoise;

  std::vector<int> widths(num_bands, kW / num_bands);
  for (int i = 0; i < kW % num_bands; ++i) ++widths[i];

  Base base;
  base.rgb.p.assign(static_cast<size_t>(kW) * kH * 3, 0);
  base.region.resize(static_cast<size_t>(kW) * kH);
  int xstart = 0;
  for (int bi = 0; bi < num_bands; ++bi) {
    const int kind = bands[bi];
    const int band_width = widths[bi];
    double color0[3];
    double color1[3];
    for (int c = 0; c < 3; ++c) {
      color0[c] = 35 + rng.Range(186);
      color1[c] = std::max(
          12, std::min(243, static_cast<int>(color0[c]) + rng.Range(171) -
                                85));
    }
    for (int y = 0; y < kH; ++y) {
      for (int x = 0; x < band_width; ++x) {
        const double u = x / static_cast<double>(std::max(1, band_width - 1));
        const double v = y / static_cast<double>(kH - 1);
        double px[3];
        if (kind == kSmooth) {
          for (int c = 0; c < 3; ++c) {
            px[c] = color0[c] + 4 * (color1[c] - color0[c]) * v;
          }
        } else if (kind == kGradient) {
          const double wave =
              std::sin(2.0 * kPi * (1.7 * v + 0.25 * u));
          for (int c = 0; c < 3; ++c) {
            px[c] = color0[c] * (1.0 - u) + color1[c] * u + 10.0 * wave;
          }
        } else {
          const double low =
              0.75 + 0.25 * std::sin(2.0 * kPi * (0.7 * u + 0.5 * v));
          for (int c = 0; c < 3; ++c) {
            px[c] = color0[c] * low + rng.Signed() * 38.0;
          }
        }
        const int gx = xstart + x;
        uint8_t* out = base.rgb.At(gx, y);
        for (int c = 0; c < 3; ++c) out[c] = Clamp255(px[c]);
        base.region[static_cast<size_t>(y) * kW + gx] = kind;
      }
    }
    xstart += band_width;
  }
  return base;
}

void RgbToYcc(const Rgb& rgb, std::vector<double>* ycc) {
  ycc->resize(static_cast<size_t>(kW) * kH * 3);
  for (int i = 0; i < kW * kH; ++i) {
    const double r = rgb.p[3 * i + 0];
    const double g = rgb.p[3 * i + 1];
    const double b = rgb.p[3 * i + 2];
    (*ycc)[3 * i + 0] = 0.299 * r + 0.587 * g + 0.114 * b;
    (*ycc)[3 * i + 1] = 128.0 - 0.168736 * r - 0.331264 * g + 0.5 * b;
    (*ycc)[3 * i + 2] = 128.0 + 0.5 * r - 0.418688 * g - 0.081312 * b;
  }
}

Rgb YccToRgb(const std::vector<double>& ycc) {
  Rgb rgb;
  for (int i = 0; i < kW * kH; ++i) {
    const double y = ycc[3 * i + 0];
    const double cb = ycc[3 * i + 1] - 128.0;
    const double cr = ycc[3 * i + 2] - 128.0;
    rgb.p[3 * i + 0] = Clamp255(y + 1.402 * cr);
    rgb.p[3 * i + 1] = Clamp255(y - 0.344136 * cb - 0.714136 * cr);
    rgb.p[3 * i + 2] = Clamp255(y + 1.772 * cb);
  }
  return rgb;
}

struct Patch {
  int x0;
  int y0;
  int x1;
  int y1;
};

Patch MakePatch(int level, Rng* rng) {
  static const double fractions[] = {0.02, 0.05, 0.10, 0.20, 0.40};
  const int side = std::max(
      4, static_cast<int>(std::round(std::sqrt(fractions[level]) * kW)));
  Patch patch;
  patch.x0 = rng->Range(kW - side + 1);
  patch.y0 = rng->Range(kH - side + 1);
  patch.x1 = patch.x0 + side;
  patch.y1 = patch.y0 + side;
  return patch;
}

Rgb Distort(const Rgb& base, int family, int level, uint64_t seed) {
  std::vector<double> ycc;
  RgbToYcc(base, &ycc);
  Rng rng(seed);
  const Patch patch = MakePatch(level, &rng);
  if (family == 0) {
    Rng noise(seed + 77);
    for (int y = patch.y0; y < patch.y1; ++y) {
      for (int x = patch.x0; x < patch.x1; ++x) {
        ycc[3 * (static_cast<size_t>(y) * kW + x)] +=
            12.0 * noise.Gaussian();
      }
    }
    return YccToRgb(ycc);
  }
  if (family == 1) {
    Rgb out = base;
    Rng noise(seed + 123);
    for (int by = patch.y0; by < patch.y1; by += 8) {
      for (int bx = patch.x0; bx < patch.x1; bx += 8) {
        const int ex = std::min(patch.x1, bx + 8);
        const int ey = std::min(patch.y1, by + 8);
        double mean[3] = {0.0, 0.0, 0.0};
        int count = 0;
        for (int y = by; y < ey; ++y) {
          for (int x = bx; x < ex; ++x) {
            const uint8_t* px = base.At(x, y);
            for (int c = 0; c < 3; ++c) mean[c] += px[c];
            ++count;
          }
        }
        for (int c = 0; c < 3; ++c) mean[c] /= count;
        for (int y = by; y < ey; ++y) {
          for (int x = bx; x < ex; ++x) {
            uint8_t* px = out.At(x, y);
            for (int c = 0; c < 3; ++c) {
              px[c] = Clamp255(mean[c] + 2.0 * noise.Signed());
            }
          }
        }
      }
    }
    return out;
  }
  if (family == 2) {
    for (int y = patch.y0; y < patch.y1; ++y) {
      for (int x = patch.x0; x < patch.x1; ++x) {
        double* pixel = &ycc[3 * (static_cast<size_t>(y) * kW + x)];
        pixel[0] = std::round(pixel[0] / 16.0) * 16.0;
      }
    }
    return YccToRgb(ycc);
  }
  if (family == 3) {
    for (int y = patch.y0; y < patch.y1; ++y) {
      for (int x = patch.x0; x < patch.x1; ++x) {
        double* chroma = &ycc[3 * (static_cast<size_t>(y) * kW + x)];
        chroma[1] = std::round(chroma[1] / 32.0) * 32.0;
        chroma[2] = std::round(chroma[2] / 32.0) * 32.0;
      }
    }
    return YccToRgb(ycc);
  }

  Rgb out = base;
  constexpr int kRadius = 3;
  for (int y = patch.y0; y < patch.y1; ++y) {
    for (int x = patch.x0; x < patch.x1; ++x) {
      int sum[3] = {0, 0, 0};
      int count = 0;
      for (int dy = -kRadius; dy <= kRadius; ++dy) {
        for (int dx = -kRadius; dx <= kRadius; ++dx) {
          const int xx = std::max(0, std::min(kW - 1, x + dx));
          const int yy = std::max(0, std::min(kH - 1, y + dy));
          const uint8_t* px = base.At(xx, yy);
          for (int c = 0; c < 3; ++c) sum[c] += px[c];
          ++count;
        }
      }
      uint8_t* px = out.At(x, y);
      for (int c = 0; c < 3; ++c) px[c] = (sum[c] + count / 2) / count;
    }
  }
  return out;
}

double ProxyError(const Rgb& base, const Rgb& dist,
                  const std::vector<uint8_t>& region, int family) {
  std::vector<double> ycc0;
  std::vector<double> ycc1;
  RgbToYcc(base, &ycc0);
  RgbToYcc(dist, &ycc1);
  static const double kLumaVisibility[] = {1.0, 0.92, 0.42};
  static const double kChromaVisibility[] = {0.70, 0.52, 0.22};
  static const double kRegionMix[] = {0.32, 0.40, 0.28};
  static const double kFamilyScale[][2] = {
      {0.95, 0.45}, {1.05, 0.50}, {1.15, 0.55},
      {0.35, 1.10}, {0.85, 0.40},
  };
  double sums[3][2] = {{0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}};
  int counts[3] = {0, 0, 0};
  for (int i = 0; i < kW * kH; ++i) {
    const int r = region[i];
    const double dy = ycc0[3 * i] - ycc1[3 * i];
    const double dcb = ycc0[3 * i + 1] - ycc1[3 * i + 1];
    const double dcr = ycc0[3 * i + 2] - ycc1[3 * i + 2];
    sums[r][0] += dy * dy;
    sums[r][1] += 0.5 * (dcb * dcb + dcr * dcr);
    ++counts[r];
  }
  double error2 = 0.0;
  for (int r = 0; r < 3; ++r) {
    if (counts[r] == 0) continue;
    const double luma_error2 = sums[r][0] / counts[r];
    const double chroma_error2 = sums[r][1] / counts[r];
    error2 += kRegionMix[r] *
              (kFamilyScale[family][0] * kLumaVisibility[r] * luma_error2 +
               kFamilyScale[family][1] * kChromaVisibility[r] *
                   chroma_error2);
  }
  return std::sqrt(std::max(0.0, error2));
}

std::vector<ImageF> LinearPlanes(const Rgb& rgb) {
  static double lut[256];
  static bool initialized = false;
  if (!initialized) {
    for (int i = 0; i < 256; ++i) {
      const double srgb = i / 255.0;
      lut[i] = 255.0 * (srgb <= 0.04045
                            ? srgb / 12.92
                            : std::pow((srgb + 0.055) / 1.055, 2.4));
    }
    initialized = true;
  }
  std::vector<ImageF> planes;
  planes.reserve(3);
  for (int c = 0; c < 3; ++c) {
    planes.emplace_back(kW, kH);
    for (int y = 0; y < kH; ++y) {
      for (int x = 0; x < kW; ++x) {
        planes[c].Row(y)[x] = lut[rgb.At(x, y)[c]];
      }
    }
  }
  return planes;
}

std::vector<double> AverageRanks(const std::vector<double>& values) {
  std::vector<int> order(values.size());
  for (size_t i = 0; i < order.size(); ++i) order[i] = static_cast<int>(i);
  std::sort(order.begin(), order.end(), [&](int a, int b) {
    return values[a] != values[b] ? values[a] < values[b] : a < b;
  });
  std::vector<double> ranks(values.size());
  for (size_t begin = 0; begin < order.size();) {
    size_t end = begin + 1;
    while (end < order.size() &&
           values[order[end]] == values[order[begin]]) {
      ++end;
    }
    const double rank = 0.5 * (begin + end - 1);
    for (size_t i = begin; i < end; ++i) ranks[order[i]] = rank;
    begin = end;
  }
  return ranks;
}

double Pearson(const std::vector<double>& a,
               const std::vector<double>& b) {
  const double mean_a = std::accumulate(a.begin(), a.end(), 0.0) / a.size();
  const double mean_b = std::accumulate(b.begin(), b.end(), 0.0) / b.size();
  double covariance = 0.0;
  double variance_a = 0.0;
  double variance_b = 0.0;
  for (size_t i = 0; i < a.size(); ++i) {
    const double da = a[i] - mean_a;
    const double db = b[i] - mean_b;
    covariance += da * db;
    variance_a += da * da;
    variance_b += db * db;
  }
  return covariance / std::sqrt(variance_a * variance_b);
}

double Spearman(const std::vector<double>& a,
                const std::vector<double>& b) {
  return Pearson(AverageRanks(a), AverageRanks(b));
}

struct Case {
  int family;
  double proxy;
  std::vector<ImageF> image0;
  std::vector<ImageF> image1;
};

}  // namespace

int main(int argc, char** argv) {
  const char* csv_path = argc > 1 ? argv[1] : nullptr;
  static const char* kFamilyNames[] = {
      "local_noise", "local_block", "local_luma",
      "local_chroma", "local_blur",
  };

  std::vector<Case> cases;
  cases.reserve(kBases * kFamilies * kLevels);
  for (int b = 0; b < kBases; ++b) {
    const Base base = MakeBase(2000 + b);
    for (int f = 0; f < kFamilies; ++f) {
      for (int l = 0; l < kLevels; ++l) {
        const uint64_t seed = 3000 + b * 100 + f * 10 + l;
        const Rgb dist = Distort(base.rgb, f, l, seed);
        Case sample;
        sample.family = f;
        sample.proxy = ProxyError(base.rgb, dist, base.region, f);
        const double jitter =
            1.0 + (((b * 17 + f * 997 + l * 31) * 7919) % 1000 - 500) /
                      10000.0;
        sample.proxy *= jitter;
        sample.image0 = LinearPlanes(base.rgb);
        sample.image1 = LinearPlanes(dist);
        cases.push_back(std::move(sample));
      }
    }
  }

  std::ofstream csv;
  if (csv_path != nullptr) {
    csv.open(csv_path);
    csv << "family,proxy_error,legacy_score,region_score,smooth,gradient,"
           "noise,luma,chroma,legacy_ms,region_ms\n";
  }

  std::vector<double> proxy;
  std::vector<double> legacy;
  std::vector<double> region_score;
  std::vector<double> legacy_ms;
  std::vector<double> region_ms;
  std::vector<int> family_ids;
  double max_partition_error = 0.0;

  for (size_t i = 0; i < cases.size(); ++i) {
    const Case& sample = cases[i];
    ImageF diffmap;
    double old_score = 0.0;
    ButteraugliRegionScores score;
    ButteraugliInterface(sample.image0, sample.image1, 1.0, diffmap,
                         old_score);
    ButteraugliRegionInterface(sample.image0, sample.image1, 1.0, score);

    const auto t0 = std::chrono::steady_clock::now();
    double timed_old = 0.0;
    ButteraugliInterface(sample.image0, sample.image1, 1.0, diffmap,
                         timed_old);
    const auto t1 = std::chrono::steady_clock::now();
    ButteraugliRegionScores timed_score;
    ButteraugliRegionInterface(sample.image0, sample.image1, 1.0,
                               timed_score);
    const auto t2 = std::chrono::steady_clock::now();
    ButteraugliRegionInterface(sample.image0, sample.image1, 1.0,
                               timed_score);
    const auto t3 = std::chrono::steady_clock::now();
    ButteraugliInterface(sample.image0, sample.image1, 1.0, diffmap,
                         timed_old);
    const auto t4 = std::chrono::steady_clock::now();
    max_partition_error = std::max(
        max_partition_error,
        std::max(std::fabs(score.total -
                           (score.smooth + score.gradient + score.noise)),
                 std::fabs(score.total - (score.luma + score.chroma))));
    proxy.push_back(sample.proxy);
    legacy.push_back(old_score);
    region_score.push_back(score.total);
    family_ids.push_back(sample.family);
    legacy_ms.push_back(
        0.5 * (std::chrono::duration<double, std::milli>(t1 - t0).count() +
               std::chrono::duration<double, std::milli>(t4 - t3).count()));
    region_ms.push_back(
        0.5 * (std::chrono::duration<double, std::milli>(t2 - t1).count() +
               std::chrono::duration<double, std::milli>(t3 - t2).count()));
    if (csv.is_open()) {
      csv << kFamilyNames[sample.family] << ',' << sample.proxy << ','
          << old_score << ',' << score.total << ',' << score.smooth << ','
          << score.gradient << ',' << score.noise << ',' << score.luma << ','
          << score.chroma << ',' << legacy_ms.back() << ','
          << region_ms.back() << '\n';
    }
  }

  const double old_spearman = Spearman(proxy, legacy);
  const double new_spearman = Spearman(proxy, region_score);
  double old_family_spearman = 0.0;
  double new_family_spearman = 0.0;
  for (int f = 0; f < kFamilies; ++f) {
    std::vector<double> p;
    std::vector<double> o;
    std::vector<double> n;
    for (size_t i = 0; i < cases.size(); ++i) {
      if (family_ids[i] != f) continue;
      p.push_back(proxy[i]);
      o.push_back(legacy[i]);
      n.push_back(region_score[i]);
    }
    old_family_spearman += Spearman(p, o);
    new_family_spearman += Spearman(p, n);
  }
  old_family_spearman /= kFamilies;
  new_family_spearman /= kFamilies;

  ButteraugliRegionScores first_score;
  ButteraugliRegionScores second_score;
  ButteraugliRegionInterface(cases[0].image0, cases[0].image1, 1.0,
                             first_score);
  ButteraugliRegionInterface(cases[0].image0, cases[0].image1, 1.0,
                             second_score);
  const bool deterministic =
      std::memcmp(&first_score, &second_score, sizeof(first_score)) == 0;

  ButteraugliRegionScores identical_score;
  ButteraugliRegionInterface(cases[0].image0, cases[0].image0, 1.0,
                             identical_score);
  const bool identical_zero = identical_score.total == 0.0;

  std::sort(legacy_ms.begin(), legacy_ms.end());
  std::sort(region_ms.begin(), region_ms.end());
  const double old_median = legacy_ms[legacy_ms.size() / 2];
  const double new_median = region_ms[region_ms.size() / 2];

  std::printf("cases=%zu\n", cases.size());
  std::printf("spearman_legacy=%.6f\n", old_spearman);
  std::printf("spearman_region=%.6f\n", new_spearman);
  std::printf("spearman_delta=%.6f\n", new_spearman - old_spearman);
  std::printf("family_spearman_legacy=%.6f\n", old_family_spearman);
  std::printf("family_spearman_region=%.6f\n", new_family_spearman);
  std::printf("family_spearman_delta=%.6f\n",
              new_family_spearman - old_family_spearman);
  std::printf("max_partition_error=%.3g\n", max_partition_error);
  std::printf("deterministic=%s\n", deterministic ? "true" : "false");
  std::printf("identical_images_zero=%s\n",
              identical_zero ? "true" : "false");
  std::printf("median_ms_legacy=%.6f\n", old_median);
  std::printf("median_ms_region=%.6f\n", new_median);
  std::printf("time_ratio=%.6f\n", new_median / old_median);
  const bool passes =
      new_spearman - old_spearman >= 0.1 &&
      max_partition_error <= 1e-12 && deterministic && identical_zero &&
      new_median / old_median <= 1.1;
  return passes ? 0 : 1;
}
