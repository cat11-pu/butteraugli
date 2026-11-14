#include "scorer_api.h"

#include <algorithm>
#include <array>
#include <cerrno>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#include <utility>
#include <vector>

namespace {

constexpr int kDefaultGroups = 5000;
constexpr int kLevels = 5;
constexpr double kAlphas[kLevels] = {0.0, 0.125, 0.25, 0.5, 0.75};

const char *const kBaseKindNames[] = {
    "smooth", "gradient", "checker", "waves", "noise", "shapes"};
const char *const kFamilyNames[] = {
    "noise", "blur", "quantize", "colorcast", "dropout", "checker"};

struct Config {
  int groups = kDefaultGroups;
  uint64_t seed = 0x3692d1ffULL;
  double min_spearman = 0.95;
  double rank_tolerance = 0.02;
  bool shrink = true;
  bool thorough_shrink = true;
  int max_artifacts = 200;
  std::string out_dir = "artifacts/differential-regression";
  std::string reference_label = "git HEAD";
  std::string candidate_label = "working tree";
  std::string reference_sha;
  std::string candidate_sha;
};

struct Rng {
  uint64_t state;

  explicit Rng(uint64_t seed) : state(seed) {}

  uint64_t Next() {
    uint64_t z = (state += 0x9e3779b97f4a7c15ULL);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
  }

  int Range(int n) {
    return n <= 0 ? 0 : static_cast<int>(Next() % static_cast<uint64_t>(n));
  }

  double Unit() {
    return static_cast<double>(Next() >> 11) *
           (1.0 / 9007199254740992.0);
  }

  double Signed() { return 2.0 * Unit() - 1.0; }
};

int Clamp255(int value) {
  return value < 0 ? 0 : (value > 255 ? 255 : value);
}

struct RgbImage {
  int w = 0;
  int h = 0;
  std::vector<uint8_t> data;

  RgbImage() = default;
  RgbImage(int width, int height)
      : w(width), h(height), data(static_cast<size_t>(width) * height * 3) {}

  uint8_t *At(int x, int y) {
    return data.data() + (static_cast<size_t>(y) * w + x) * 3;
  }

  const uint8_t *At(int x, int y) const {
    return data.data() + (static_cast<size_t>(y) * w + x) * 3;
  }
};

struct LinearImage {
  int w = 0;
  int h = 0;
  std::array<std::vector<float>, 3> plane;

  LinearImage() = default;
  LinearImage(int width, int height) : w(width), h(height) {
    const size_t n = static_cast<size_t>(width) * height;
    for (auto &p : plane) p.resize(n);
  }
};

struct Crop {
  int x = 0;
  int y = 0;
  int w = 0;
  int h = 0;

  Crop() = default;
  Crop(int x_value, int y_value, int w_value, int h_value)
      : x(x_value), y(y_value), w(w_value), h(h_value) {}

  int64_t Area() const { return static_cast<int64_t>(w) * h; }
};

struct GroupSpec {
  uint64_t seed = 0;
  uint64_t target_seed = 0;
  int w = 0;
  int h = 0;
  int base_kind = 0;
  int family = 0;
  int amplitude = 0;
};

struct ScoreResult {
  bool ok = false;
  double value = 0.0;
  uint64_t bits = 0;
};

struct Record {
  int group = 0;
  int level = 0;
  int family = 0;
  ScoreResult reference;
  ScoreResult candidate;
  ScoreResult candidate_repeat;
  bool deterministic = false;
};

struct SpearmanResult {
  size_t n = 0;
  double rho = 0.0;
  bool defined = false;
};

struct Counterexample {
  std::string kind;
  int group = 0;
  int level = 0;
  std::string family;
  uint64_t group_seed = 0;
  uint64_t target_seed = 0;
  std::string base_kind;
  int amplitude = 0;
  double alpha_lo = 0.0;
  double alpha_hi = 0.0;
  Crop crop;
  int original_w = 0;
  int original_h = 0;
  double reference_lo = 0.0;
  double reference_hi = 0.0;
  double candidate_lo = 0.0;
  double candidate_hi = 0.0;
  double metric = 0.0;
  double threshold = 0.0;
  std::string note;
};

struct ExtremeResult {
  int w = 0;
  int h = 0;
  bool candidate_crash = false;
  bool reference_crash = false;
  int candidate_signal = 0;
  int reference_signal = 0;
};

const std::array<double, 256> &SrgbToLinearTable() {
  static const std::array<double, 256> table = [] {
    std::array<double, 256> values{};
    for (int i = 0; i < 256; ++i) {
      const double srgb = i / 255.0;
      values[i] =
          255.0 * (srgb <= 0.04045
                       ? srgb / 12.92
                       : std::pow((srgb + 0.055) / 1.055, 2.4));
    }
    return values;
  }();
  return table;
}

LinearImage ToLinear(const RgbImage &image) {
  LinearImage linear(image.w, image.h);
  const auto &table = SrgbToLinearTable();
  for (int y = 0; y < image.h; ++y) {
    for (int x = 0; x < image.w; ++x) {
      const uint8_t *pixel = image.At(x, y);
      const size_t index = static_cast<size_t>(y) * image.w + x;
      for (int c = 0; c < 3; ++c) {
        linear.plane[c][index] = static_cast<float>(table[pixel[c]]);
      }
    }
  }
  return linear;
}

LinearImage CropLinear(const LinearImage &image, const Crop &crop) {
  LinearImage out(crop.w, crop.h);
  for (int c = 0; c < 3; ++c) {
    for (int y = 0; y < crop.h; ++y) {
      const float *src =
          image.plane[c].data() + static_cast<size_t>(crop.y + y) * image.w +
          crop.x;
      float *dst = out.plane[c].data() + static_cast<size_t>(y) * crop.w;
      std::copy(src, src + crop.w, dst);
    }
  }
  return out;
}

RgbImage CropRgb(const RgbImage &image, const Crop &crop) {
  RgbImage out(crop.w, crop.h);
  for (int y = 0; y < crop.h; ++y) {
    for (int x = 0; x < crop.w; ++x) {
      const uint8_t *src = image.At(crop.x + x, crop.y + y);
      uint8_t *dst = out.At(x, y);
      dst[0] = src[0];
      dst[1] = src[1];
      dst[2] = src[2];
    }
  }
  return out;
}

int ScoreImages(const LinearImage &a, const LinearImage &b, bool reference,
                ScoreResult *result) {
  double value = 0.0;
  const int ok =
      reference
          ? ButteraugliReferenceScore(
                a.plane[0].data(), a.plane[1].data(), a.plane[2].data(),
                b.plane[0].data(), b.plane[1].data(), b.plane[2].data(), a.w,
                a.h, &value)
          : ButteraugliCandidateScore(
                a.plane[0].data(), a.plane[1].data(), a.plane[2].data(),
                b.plane[0].data(), b.plane[1].data(), b.plane[2].data(), a.w,
                a.h, &value);
  result->ok = ok != 0;
  result->value = ok != 0 ? value : 0.0;
  result->bits = 0;
  if (ok != 0) std::memcpy(&result->bits, &value, sizeof(value));
  return ok;
}

void FillGrid(Rng *rng, int grid_w, int grid_h,
              std::vector<uint8_t> *grid) {
  grid->resize(static_cast<size_t>(grid_w) * grid_h * 3);
  for (auto &value : *grid) value = static_cast<uint8_t>(rng->Range(256));
}

double Bilinear(const std::vector<uint8_t> &grid, int grid_w, int grid_h,
                int channel, double gx, double gy) {
  const int x0 =
      std::max(0, std::min(grid_w - 1, static_cast<int>(std::floor(gx))));
  const int y0 =
      std::max(0, std::min(grid_h - 1, static_cast<int>(std::floor(gy))));
  const int x1 = std::min(grid_w - 1, x0 + 1);
  const int y1 = std::min(grid_h - 1, y0 + 1);
  const double fx = gx - x0;
  const double fy = gy - y0;
  const auto at = [&](int x, int y) {
    return static_cast<double>(
        grid[(static_cast<size_t>(y) * grid_w + x) * 3 + channel]);
  };
  const double top = at(x0, y0) * (1.0 - fx) + at(x1, y0) * fx;
  const double bottom = at(x0, y1) * (1.0 - fx) + at(x1, y1) * fx;
  return top * (1.0 - fy) + bottom * fy;
}

RgbImage MakeBase(const GroupSpec &group) {
  RgbImage image(group.w, group.h);
  if (group.w <= 0 || group.h <= 0) return image;
  Rng rng(group.seed ^ 0x51ed270d2f83a1c9ULL);

  switch (group.base_kind) {
    case 0: {  // Smooth, natural-looking low-frequency field.
      const int gw = 4 + rng.Range(3);
      const int gh = 4 + rng.Range(3);
      std::vector<uint8_t> grid;
      FillGrid(&rng, gw, gh, &grid);
      const double sx = group.w > 1 ? static_cast<double>(group.w - 1) /
                                          std::max(1, gw - 1)
                                    : 0.0;
      const double sy = group.h > 1 ? static_cast<double>(group.h - 1) /
                                          std::max(1, gh - 1)
                                    : 0.0;
      for (int y = 0; y < group.h; ++y) {
        for (int x = 0; x < group.w; ++x) {
          uint8_t *pixel = image.At(x, y);
          for (int c = 0; c < 3; ++c) {
            const double gx = group.w > 1 ? x / sx : 0.0;
            const double gy = group.h > 1 ? y / sy : 0.0;
            const double v =
                Bilinear(grid, gw, gh, c, gx, gy) + rng.Signed() * 2.0;
            pixel[c] = static_cast<uint8_t>(Clamp255(
                static_cast<int>(std::lround(v))));
          }
        }
      }
      break;
    }
    case 1: {  // Two-corner color gradient with a mild diagonal modulation.
      std::array<int, 3> c0;
      std::array<int, 3> c1;
      for (int c = 0; c < 3; ++c) {
        c0[c] = rng.Range(256);
        c1[c] = rng.Range(256);
      }
      for (int y = 0; y < group.h; ++y) {
        for (int x = 0; x < group.w; ++x) {
          const double tx = group.w > 1 ? static_cast<double>(x) /
                                              (group.w - 1)
                                        : 0.0;
          const double ty = group.h > 1 ? static_cast<double>(y) /
                                              (group.h - 1)
                                        : 0.0;
          const double t = std::min(1.0, 0.7 * tx + 0.3 * ty);
          uint8_t *pixel = image.At(x, y);
          for (int c = 0; c < 3; ++c) {
            const double wave =
                10.0 * std::sin(6.283185307179586 * (tx + 0.7 * ty + c * 0.2));
            pixel[c] =
                static_cast<uint8_t>(Clamp255(static_cast<int>(std::lround(
                    c0[c] * (1.0 - t) + c1[c] * t + wave))));
          }
        }
      }
      break;
    }
    case 2: {  // Checkerboard, with independent channel colors.
      const int tile = 1 + rng.Range(12);
      std::array<int, 3> c0;
      std::array<int, 3> c1;
      for (int c = 0; c < 3; ++c) {
        c0[c] = 24 + rng.Range(208);
        c1[c] = 24 + rng.Range(208);
      }
      for (int y = 0; y < group.h; ++y) {
        for (int x = 0; x < group.w; ++x) {
          const int parity = ((x / tile) + (y / tile)) & 1;
          uint8_t *pixel = image.At(x, y);
          for (int c = 0; c < 3; ++c) {
            pixel[c] = static_cast<uint8_t>(
                Clamp255((parity ? c1[c] : c0[c]) + rng.Range(7) - 3));
          }
        }
      }
      break;
    }
    case 3: {  // Superposition of waves.
      std::array<double, 3> fx;
      std::array<double, 3> fy;
      std::array<double, 3> phase;
      for (int c = 0; c < 3; ++c) {
        fx[c] = 1.0 + rng.Range(5);
        fy[c] = 1.0 + rng.Range(5);
        phase[c] = rng.Unit() * 6.283185307179586;
      }
      for (int y = 0; y < group.h; ++y) {
        for (int x = 0; x < group.w; ++x) {
          uint8_t *pixel = image.At(x, y);
          for (int c = 0; c < 3; ++c) {
            const double tx = group.w > 1 ? static_cast<double>(x) /
                                                (group.w - 1)
                                          : 0.0;
            const double ty = group.h > 1 ? static_cast<double>(y) /
                                                (group.h - 1)
                                          : 0.0;
            const double v =
                128.0 + 55.0 * std::sin(6.283185307179586 *
                                        (fx[c] * tx + fy[c] * ty) +
                                        phase[c]);
            pixel[c] =
                static_cast<uint8_t>(Clamp255(static_cast<int>(std::lround(v))));
          }
        }
      }
      break;
    }
    case 4: {  // Low-frequency background plus per-pixel texture.
      const int gw = 3;
      const int gh = 3;
      std::vector<uint8_t> grid;
      FillGrid(&rng, gw, gh, &grid);
      const double sx =
          group.w > 1 ? static_cast<double>(group.w - 1) / (gw - 1) : 0.0;
      const double sy =
          group.h > 1 ? static_cast<double>(group.h - 1) / (gh - 1) : 0.0;
      for (int y = 0; y < group.h; ++y) {
        for (int x = 0; x < group.w; ++x) {
          uint8_t *pixel = image.At(x, y);
          for (int c = 0; c < 3; ++c) {
            const double gx = group.w > 1 ? x / sx : 0.0;
            const double gy = group.h > 1 ? y / sy : 0.0;
            const double smooth = Bilinear(grid, gw, gh, c, gx, gy);
            const double v = smooth + rng.Signed() * 45.0;
            pixel[c] = static_cast<uint8_t>(
                Clamp255(static_cast<int>(std::lround(v))));
          }
        }
      }
      break;
    }
    default: {  // Random rectangles.
      std::array<int, 3> background;
      for (int c = 0; c < 3; ++c) background[c] = rng.Range(256);
      for (int y = 0; y < group.h; ++y) {
        for (int x = 0; x < group.w; ++x) {
          uint8_t *pixel = image.At(x, y);
          for (int c = 0; c < 3; ++c) pixel[c] = background[c];
        }
      }
      const int rectangles = 2 + rng.Range(5);
      for (int r = 0; r < rectangles; ++r) {
        const int x0 = rng.Range(group.w);
        const int y0 = rng.Range(group.h);
        const int rw = 1 + rng.Range(std::max(1, group.w - x0));
        const int rh = 1 + rng.Range(std::max(1, group.h - y0));
        std::array<int, 3> color;
        for (int c = 0; c < 3; ++c) color[c] = rng.Range(256);
        for (int y = y0; y < std::min(group.h, y0 + rh); ++y) {
          for (int x = x0; x < std::min(group.w, x0 + rw); ++x) {
            uint8_t *pixel = image.At(x, y);
            for (int c = 0; c < 3; ++c) pixel[c] = color[c];
          }
        }
      }
      break;
    }
  }
  return image;
}

RgbImage BoxBlur(const RgbImage &input, int radius) {
  RgbImage output(input.w, input.h);
  if (input.w <= 0 || input.h <= 0) return output;
  radius = std::max(1, radius);
  for (int y = 0; y < input.h; ++y) {
    for (int x = 0; x < input.w; ++x) {
      int sum[3] = {0, 0, 0};
      int count = 0;
      for (int dy = -radius; dy <= radius; ++dy) {
        const int yy = std::min(input.h - 1, std::max(0, y + dy));
        for (int dx = -radius; dx <= radius; ++dx) {
          const int xx = std::min(input.w - 1, std::max(0, x + dx));
          const uint8_t *pixel = input.At(xx, yy);
          for (int c = 0; c < 3; ++c) sum[c] += pixel[c];
          ++count;
        }
      }
      uint8_t *pixel = output.At(x, y);
      for (int c = 0; c < 3; ++c) {
        pixel[c] =
            static_cast<uint8_t>((sum[c] + count / 2) / std::max(1, count));
      }
    }
  }
  return output;
}

RgbImage MakeTarget(const RgbImage &base, const GroupSpec &group) {
  RgbImage target(base.w, base.h);
  Rng rng(group.target_seed ^ 0xa0761d6478bd642fULL);
  switch (group.family) {
    case 0: {  // Additive Gaussian noise.
      for (int y = 0; y < base.h; ++y) {
        for (int x = 0; x < base.w; ++x) {
          uint8_t *dst = target.At(x, y);
          const uint8_t *src = base.At(x, y);
          for (int c = 0; c < 3; ++c) {
            const double u1 = std::max(1e-12, rng.Unit());
            const double u2 = rng.Unit();
            const double gaussian =
                std::sqrt(-2.0 * std::log(u1)) *
                std::cos(6.283185307179586 * u2);
            dst[c] = static_cast<uint8_t>(Clamp255(static_cast<int>(
                std::lround(src[c] + gaussian * group.amplitude))));
          }
        }
      }
      break;
    }
    case 1:
      return BoxBlur(base, group.amplitude);
    case 2: {  // Posterization.
      const int step = std::max(2, group.amplitude);
      for (int y = 0; y < base.h; ++y) {
        for (int x = 0; x < base.w; ++x) {
          const uint8_t *src = base.At(x, y);
          uint8_t *dst = target.At(x, y);
          for (int c = 0; c < 3; ++c) {
            dst[c] =
                static_cast<uint8_t>(((src[c] + step / 2) / step) * step);
          }
        }
      }
      break;
    }
    case 3: {  // A three-channel color cast.
      const int delta[3] = {group.amplitude, -group.amplitude / 2,
                            group.amplitude / 3};
      for (int y = 0; y < base.h; ++y) {
        for (int x = 0; x < base.w; ++x) {
          const uint8_t *src = base.At(x, y);
          uint8_t *dst = target.At(x, y);
          for (int c = 0; c < 3; ++c) {
            dst[c] = static_cast<uint8_t>(Clamp255(src[c] + delta[c]));
          }
        }
      }
      break;
    }
    case 4: {  // Random pixel dropout.
      target = base;
      const int percent = std::min(90, std::max(1, group.amplitude));
      for (int y = 0; y < base.h; ++y) {
        for (int x = 0; x < base.w; ++x) {
          if (rng.Range(100) >= percent) continue;
          uint8_t *dst = target.At(x, y);
          const int value = rng.Range(2) ? 0 : 255;
          dst[0] = dst[1] = dst[2] = static_cast<uint8_t>(value);
        }
      }
      break;
    }
    default: {  // Independent high-contrast checker pattern.
      const int tile = std::max(1, group.amplitude);
      std::array<int, 3> c0;
      std::array<int, 3> c1;
      for (int c = 0; c < 3; ++c) {
        c0[c] = rng.Range(256);
        c1[c] = rng.Range(256);
      }
      for (int y = 0; y < base.h; ++y) {
        for (int x = 0; x < base.w; ++x) {
          const int parity = ((x / tile) + (y / tile)) & 1;
          uint8_t *dst = target.At(x, y);
          for (int c = 0; c < 3; ++c) dst[c] = parity ? c1[c] : c0[c];
        }
      }
      break;
    }
  }
  return target;
}

RgbImage Mix(const RgbImage &base, const RgbImage &target, double alpha) {
  if (alpha <= 0.0) return base;
  if (alpha >= 1.0) return target;
  RgbImage mixed(base.w, base.h);
  for (size_t i = 0; i < mixed.data.size(); ++i) {
    const double value =
        base.data[i] * (1.0 - alpha) + target.data[i] * alpha;
    mixed.data[i] = static_cast<uint8_t>(
        Clamp255(static_cast<int>(std::lround(value))));
  }
  return mixed;
}

GroupSpec MakeGroup(int index, Rng *rng) {
  GroupSpec group;
  group.seed = rng->Next();
  group.target_seed = rng->Next();
  const bool narrow = rng->Range(16) == 0;
  group.w = narrow ? 1 + rng->Range(4) : 4 + rng->Range(61);
  group.h = narrow ? 1 + rng->Range(4) : 4 + rng->Range(61);
  group.base_kind = rng->Range(6);
  group.family = rng->Range(6);
  switch (group.family) {
    case 0:
      group.amplitude = 24 + rng->Range(73);
      break;
    case 1:
      group.amplitude = 2 + rng->Range(7);
      break;
    case 2:
      group.amplitude = 8 << rng->Range(3);
      break;
    case 3:
      group.amplitude = 32 + rng->Range(97);
      break;
    case 4:
      group.amplitude = 5 + rng->Range(26);
      break;
    default:
      group.amplitude = 2 + rng->Range(11);
      break;
  }
  (void)index;
  return group;
}

struct PairContext {
  RgbImage base;
  RgbImage distorted;
  LinearImage base_linear;
  LinearImage distorted_linear;
};

PairContext MakePairContext(const GroupSpec &group, int level) {
  PairContext context;
  context.base = MakeBase(group);
  const RgbImage target = MakeTarget(context.base, group);
  context.distorted = Mix(context.base, target, kAlphas[level]);
  context.base_linear = ToLinear(context.base);
  context.distorted_linear = ToLinear(context.distorted);
  return context;
}

double MonotonicTolerance(double previous) {
  return 1e-7 + 1e-5 * std::max(1.0, std::fabs(previous));
}

std::vector<double> AverageRanks(const std::vector<double> &values,
                                 const std::vector<size_t> &indices) {
  std::vector<size_t> order = indices;
  std::sort(order.begin(), order.end(), [&](size_t a, size_t b) {
    if (values[a] != values[b]) return values[a] < values[b];
    return a < b;
  });
  std::vector<double> ranks(values.size(),
                            std::numeric_limits<double>::quiet_NaN());
  size_t i = 0;
  while (i < order.size()) {
    size_t j = i + 1;
    while (j < order.size() &&
           values[order[j]] == values[order[i]]) {
      ++j;
    }
    const double rank = (static_cast<double>(i) +
                         static_cast<double>(j - 1)) *
                        0.5;
    for (size_t k = i; k < j; ++k) ranks[order[k]] = rank;
    i = j;
  }
  return ranks;
}

double Pearson(const std::vector<double> &a, const std::vector<double> &b) {
  if (a.empty() || a.size() != b.size()) return 0.0;
  double mean_a = 0.0;
  double mean_b = 0.0;
  for (size_t i = 0; i < a.size(); ++i) {
    mean_a += a[i];
    mean_b += b[i];
  }
  mean_a /= a.size();
  mean_b /= b.size();
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
  if (variance_a == 0.0 || variance_b == 0.0) {
    return variance_a == 0.0 && variance_b == 0.0 ? 1.0 : 0.0;
  }
  return covariance / std::sqrt(variance_a * variance_b);
}

SpearmanResult SpearmanFor(const std::vector<Record> &records,
                           const std::vector<size_t> &indices,
                           bool reference_first) {
  SpearmanResult result;
  std::vector<double> a;
  std::vector<double> b;
  a.reserve(indices.size());
  b.reserve(indices.size());
  for (size_t index : indices) {
    const Record &record = records[index];
    const double x = reference_first ? record.reference.value
                                     : record.candidate.value;
    const double y = reference_first ? record.candidate.value
                                     : record.reference.value;
    a.push_back(x);
    b.push_back(y);
  }
  std::vector<size_t> all(a.size());
  for (size_t i = 0; i < all.size(); ++i) all[i] = i;
  const std::vector<double> ranks_a = AverageRanks(a, all);
  const std::vector<double> ranks_b = AverageRanks(b, all);
  result.n = all.size();
  result.defined = all.size() >= 2;
  result.rho = result.defined ? Pearson(ranks_a, ranks_b) : 1.0;
  return result;
}

double InsertionAverageRank(const std::vector<double> &sorted_values,
                            double value) {
  const auto lower =
      std::lower_bound(sorted_values.begin(), sorted_values.end(), value);
  const auto upper =
      std::upper_bound(sorted_values.begin(), sorted_values.end(), value);
  const double first = static_cast<double>(lower - sorted_values.begin());
  const double equal_count = static_cast<double>(upper - lower) + 1.0;
  return first + (equal_count - 1.0) * 0.5;
}

uint64_t CropKey(const Crop &crop) {
  return (static_cast<uint64_t>(crop.x) << 48) |
         (static_cast<uint64_t>(crop.y) << 32) |
         (static_cast<uint64_t>(crop.w) << 16) |
         static_cast<uint64_t>(crop.h);
}

std::vector<Crop> CropCandidates(const Crop &current) {
  std::vector<Crop> candidates;
  const int w = current.w;
  const int h = current.h;
  if (w <= 0 || h <= 0 || w * h <= 1) return candidates;
  const int hw = std::max(1, w / 2);
  const int hh = std::max(1, h / 2);
  const int x = current.x;
  const int y = current.y;

  const Crop quadrants[] = {
      {x, y, hw, hh},
      {x + hw, y, w - hw, hh},
      {x, y + hh, hw, h - hh},
      {x + hw, y + hh, w - hw, h - hh},
      {x, y, hw, h},
      {x + hw, y, w - hw, h},
      {x, y, w, hh},
      {x, y + hh, w, h - hh},
      {x + std::max(0, (w - hw) / 2), y + std::max(0, (h - hh) / 2), hw, hh},
  };
  for (const Crop &candidate : quadrants) {
    if (candidate.w > 0 && candidate.h > 0 &&
        candidate.Area() < current.Area()) {
      candidates.push_back(candidate);
    }
  }

  const int sample_w = std::min(w, 16);
  const int sample_h = std::min(h, 16);
  for (int sy = 0; sy < sample_h; ++sy) {
    for (int sx = 0; sx < sample_w; ++sx) {
      const int px =
          x + static_cast<int>((static_cast<int64_t>(sx) * (w - 1)) /
                               std::max(1, sample_w - 1));
      const int py =
          y + static_cast<int>((static_cast<int64_t>(sy) * (h - 1)) /
                               std::max(1, sample_h - 1));
      candidates.push_back({px, py, 1, 1});
    }
  }
  candidates.push_back({x + w / 2, y + h / 2, 1, 1});

  if (w > 1) {
    candidates.push_back({x + 1, y, w - 1, h});
    candidates.push_back({x, y, w - 1, h});
  }
  if (h > 1) {
    candidates.push_back({x, y + 1, w, h - 1});
    candidates.push_back({x, y, w, h - 1});
  }
  return candidates;
}

template <typename Predicate>
Crop MinimizeCrop(int full_w, int full_h, bool thorough,
                  Predicate predicate) {
  Crop best{0, 0, full_w, full_h};
  if (full_w <= 0 || full_h <= 0 || !predicate(best)) return best;
  std::set<uint64_t> tested;
  tested.insert(CropKey(best));

  for (int iteration = 0; iteration < 32; ++iteration) {
    const std::vector<Crop> candidates = CropCandidates(best);
    Crop next = best;
    for (const Crop &candidate : candidates) {
      if (candidate.Area() >= next.Area()) continue;
      const uint64_t key = CropKey(candidate);
      if (!tested.insert(key).second) continue;
      if (predicate(candidate)) next = candidate;
    }
    if (next.Area() >= best.Area()) break;
    best = next;
    if (best.Area() == 1) break;
  }

  if (thorough && full_w > 0 && full_h > 0) {
    for (int y = 0; y < full_h && best.Area() > 1; ++y) {
      for (int x = 0; x < full_w; ++x) {
        const Crop pixel{x, y, 1, 1};
        const uint64_t key = CropKey(pixel);
        if (!tested.insert(key).second) continue;
        if (predicate(pixel)) {
          best = pixel;
          break;
        }
      }
    }
  }
  return best;
}

bool WritePpm(const std::string &path, const RgbImage &image) {
  std::ofstream out(path, std::ios::binary);
  if (!out) return false;
  out << "P6\n" << image.w << " " << image.h << "\n255\n";
  out.write(reinterpret_cast<const char *>(image.data.data()),
            static_cast<std::streamsize>(image.data.size()));
  return static_cast<bool>(out);
}

std::string JoinPath(const std::string &a, const std::string &b) {
  if (a.empty()) return b;
  if (a.back() == '/') return a + b;
  return a + "/" + b;
}

std::string Timestamp() {
  const auto now = std::chrono::system_clock::now();
  const std::time_t value = std::chrono::system_clock::to_time_t(now);
  std::tm tm_value{};
#ifdef _WIN32
  localtime_s(&tm_value, &value);
#else
  localtime_r(&value, &tm_value);
#endif
  char buffer[64];
  std::strftime(buffer, sizeof(buffer), "%Y-%m-%dT%H:%M:%S%z", &tm_value);
  return buffer;
}

void PrintUsage(const char *program) {
  std::cout
      << "Usage: " << program << " [options]\n"
      << "  --groups N             number of random groups (default 5000)\n"
      << "  --seed N               deterministic master seed\n"
      << "  --min-spearman X       minimum rank correlation (default 0.95)\n"
      << "  --rank-tolerance X     per-sample rank outlier fraction (0.02)\n"
      << "  --max-artifacts N      write at most N cropped PPM pairs (200)\n"
      << "  --out DIR              report directory\n"
      << "  --no-shrink            disable counterexample minimization\n"
      << "  --reference-label STR  human readable reference description\n"
      << "  --candidate-label STR  human readable candidate description\n"
      << "  --reference-sha STR    reference revision\n"
      << "  --candidate-sha STR    candidate source hash\n";
}

bool ParseConfig(int argc, char **argv, Config *config) {
  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    const auto next = [&](const char *name) -> const char * {
      if (i + 1 >= argc) {
        std::cerr << "Missing value for " << name << "\n";
        return nullptr;
      }
      return argv[++i];
    };
    if (arg == "--help" || arg == "-h") {
      PrintUsage(argv[0]);
      std::exit(0);
    } else if (arg == "--groups") {
      const char *value = next("--groups");
      if (value == nullptr) return false;
      config->groups = std::max(1, std::atoi(value));
    } else if (arg == "--seed") {
      const char *value = next("--seed");
      if (value == nullptr) return false;
      config->seed = std::strtoull(value, nullptr, 0);
    } else if (arg == "--min-spearman") {
      const char *value = next("--min-spearman");
      if (value == nullptr) return false;
      config->min_spearman = std::atof(value);
    } else if (arg == "--rank-tolerance") {
      const char *value = next("--rank-tolerance");
      if (value == nullptr) return false;
      config->rank_tolerance = std::atof(value);
    } else if (arg == "--out") {
      const char *value = next("--out");
      if (value == nullptr) return false;
      config->out_dir = value;
    } else if (arg == "--no-shrink") {
      config->shrink = false;
    } else if (arg == "--fast-shrink") {
      config->thorough_shrink = false;
    } else if (arg == "--max-artifacts") {
      const char *value = next("--max-artifacts");
      if (value == nullptr) return false;
      config->max_artifacts = std::max(0, std::atoi(value));
    } else if (arg == "--reference-label") {
      const char *value = next("--reference-label");
      if (value == nullptr) return false;
      config->reference_label = value;
    } else if (arg == "--candidate-label") {
      const char *value = next("--candidate-label");
      if (value == nullptr) return false;
      config->candidate_label = value;
    } else if (arg == "--reference-sha") {
      const char *value = next("--reference-sha");
      if (value == nullptr) return false;
      config->reference_sha = value;
    } else if (arg == "--candidate-sha") {
      const char *value = next("--candidate-sha");
      if (value == nullptr) return false;
      config->candidate_sha = value;
    } else {
      std::cerr << "Unknown option: " << arg << "\n";
      PrintUsage(argv[0]);
      return false;
    }
  }
  return true;
}

int RunExtremeChild(const char *implementation, int width, int height) {
  const size_t n = static_cast<size_t>(std::max(0, width)) *
                   static_cast<size_t>(std::max(0, height));
  std::vector<float> image0(n * 3);
  std::vector<float> image1(n * 3);
  for (size_t i = 0; i < n; ++i) {
    for (int c = 0; c < 3; ++c) {
      image0[i * 3 + c] = 0.2f + 0.1f * static_cast<float>(c);
      image1[i * 3 + c] =
          image0[i * 3 + c] + 0.01f * static_cast<float>((i + c) % 7);
    }
  }
  double score = 0.0;
  if (std::strcmp(implementation, "reference") == 0) {
    ButteraugliReferenceScore(
        image0.data(), image0.data() + n, image0.data() + 2 * n,
        image1.data(), image1.data() + n, image1.data() + 2 * n, width, height,
        &score);
  } else {
    ButteraugliCandidateScore(
        image0.data(), image0.data() + n, image0.data() + 2 * n,
        image1.data(), image1.data() + n, image1.data() + 2 * n, width, height,
        &score);
  }
  return 0;
}

ExtremeResult RunExtremeCase(const char *program, int width, int height,
                             bool reference) {
  ExtremeResult result;
  result.w = width;
  result.h = height;
  std::cout.flush();
  std::cerr.flush();
  const pid_t pid = fork();
  if (pid == 0) {
    alarm(20);
    char executable[4096];
    const ssize_t length = readlink("/proc/self/exe", executable,
                                    sizeof(executable) - 1);
    const char *path = program;
    if (length > 0) {
      executable[length] = '\0';
      path = executable;
    }
    execl(path, path, "--extreme-child", reference ? "reference"
                                                   : "candidate",
          std::to_string(width).c_str(), std::to_string(height).c_str(),
          static_cast<char *>(nullptr));
    _exit(127);
  }
  if (pid < 0) {
    std::cerr << "fork failed for extreme case " << width << "x" << height
              << ": " << std::strerror(errno) << "\n";
    if (reference) {
      result.reference_crash = true;
      result.reference_signal = -1;
    } else {
      result.candidate_crash = true;
      result.candidate_signal = -1;
    }
    return result;
  }
  int status = 0;
  bool wait_failed = false;
  while (waitpid(pid, &status, 0) < 0) {
    if (errno != EINTR) {
      wait_failed = true;
      break;
    }
  }
  const bool crashed = wait_failed || !WIFEXITED(status) ||
                       WEXITSTATUS(status) != 0 || WIFSIGNALED(status);
  const int signal =
      WIFSIGNALED(status) ? WTERMSIG(status) : (crashed ? -1 : 0);
  if (reference) {
    result.reference_crash = crashed;
    result.reference_signal = signal;
  } else {
    result.candidate_crash = crashed;
    result.candidate_signal = signal;
  }
  return result;
}

bool ScoreCrop(const PairContext &context, const Crop &crop, bool reference,
               double *score) {
  const LinearImage a = CropLinear(context.base_linear, crop);
  const LinearImage b = CropLinear(context.distorted_linear, crop);
  ScoreResult result;
  if (!ScoreImages(a, b, reference, &result)) return false;
  *score = result.value;
  return true;
}

std::vector<Counterexample> ShrinkMonotonicFailures(
    const Config &config, const std::vector<GroupSpec> &groups,
    const std::vector<Record> &records) {
  std::vector<Counterexample> counterexamples;
  for (size_t index = 0; index < records.size(); ++index) {
    const Record &current = records[index];
    if (current.level == 0 || !current.candidate.ok) continue;
    const Record &previous = records[index - 1];
    if (!previous.candidate.ok ||
        previous.candidate.value - current.candidate.value <=
            MonotonicTolerance(previous.candidate.value)) {
      continue;
    }
    const GroupSpec &group = groups[current.group];
    PairContext lo = MakePairContext(group, current.level - 1);
    PairContext hi = MakePairContext(group, current.level);
    const Crop full{0, 0, group.w, group.h};
    auto predicate = [&](const Crop &crop) {
      double a = 0.0;
      double b = 0.0;
      if (!ScoreCrop(lo, crop, false, &a) ||
          !ScoreCrop(hi, crop, false, &b)) {
        return false;
      }
      return a - b > MonotonicTolerance(a);
    };
    Crop minimal = full;
    if (config.shrink) {
      minimal = MinimizeCrop(group.w, group.h, config.thorough_shrink,
                             predicate);
    }
    double ref_lo = 0.0;
    double ref_hi = 0.0;
    double cand_lo = 0.0;
    double cand_hi = 0.0;
    if (!ScoreCrop(lo, minimal, false, &cand_lo) ||
        !ScoreCrop(hi, minimal, false, &cand_hi)) {
      cand_lo = previous.candidate.value;
      cand_hi = current.candidate.value;
      minimal = full;
    }
    ScoreCrop(lo, minimal, true, &ref_lo);
    ScoreCrop(hi, minimal, true, &ref_hi);
    Counterexample counterexample;
    counterexample.kind = "monotonicity";
    counterexample.group = current.group;
    counterexample.level = current.level;
    counterexample.family = kFamilyNames[group.family];
    counterexample.group_seed = group.seed;
    counterexample.target_seed = group.target_seed;
    counterexample.base_kind = kBaseKindNames[group.base_kind];
    counterexample.amplitude = group.amplitude;
    counterexample.alpha_lo = kAlphas[current.level - 1];
    counterexample.alpha_hi = kAlphas[current.level];
    counterexample.crop = minimal;
    counterexample.original_w = group.w;
    counterexample.original_h = group.h;
    counterexample.reference_lo = ref_lo;
    counterexample.reference_hi = ref_hi;
    counterexample.candidate_lo = cand_lo;
    counterexample.candidate_hi = cand_hi;
    counterexample.metric = cand_lo - cand_hi;
    counterexample.threshold = MonotonicTolerance(cand_lo);
    counterexample.note =
        "candidate score decreased while distortion increased";
    if (ref_lo - ref_hi > MonotonicTolerance(ref_lo)) {
      counterexample.note += "; reference violates as well";
    }
    counterexamples.push_back(counterexample);
  }
  return counterexamples;
}

std::vector<Counterexample> ShrinkRankFailures(
    const Config &config, const std::vector<GroupSpec> &groups,
    const std::vector<Record> &records, const std::vector<size_t> &valid,
    double rank_threshold,
    const std::vector<size_t> &failure_indices) {
  std::vector<Counterexample> counterexamples;
  for (size_t record_index : failure_indices) {
    const Record &record = records[record_index];
    const GroupSpec &group = groups[record.group];
    PairContext context = MakePairContext(group, record.level);
    std::vector<double> background_reference;
    std::vector<double> background_candidate;
    background_reference.reserve(valid.size() - 1);
    background_candidate.reserve(valid.size() - 1);
    for (size_t candidate_index : valid) {
      if (candidate_index == record_index) continue;
      background_reference.push_back(records[candidate_index].reference.value);
      background_candidate.push_back(records[candidate_index].candidate.value);
    }
    std::sort(background_reference.begin(), background_reference.end());
    std::sort(background_candidate.begin(), background_candidate.end());
    const double original_rank_reference =
        InsertionAverageRank(background_reference, record.reference.value);
    const double original_rank_candidate =
        InsertionAverageRank(background_candidate, record.candidate.value);
    auto predicate = [&](const Crop &crop) {
      double ref = 0.0;
      double cand = 0.0;
      if (!ScoreCrop(context, crop, true, &ref) ||
          !ScoreCrop(context, crop, false, &cand)) {
        return false;
      }
      const double rank_reference =
          InsertionAverageRank(background_reference, ref);
      const double rank_candidate =
          InsertionAverageRank(background_candidate, cand);
      return std::fabs(rank_candidate - rank_reference) > rank_threshold;
    };
    const Crop full{0, 0, group.w, group.h};
    Crop minimal = full;
    if (config.shrink) {
      minimal = MinimizeCrop(group.w, group.h, config.thorough_shrink,
                             predicate);
    }
    double ref = record.reference.value;
    double cand = record.candidate.value;
    if (!ScoreCrop(context, minimal, true, &ref) ||
        !ScoreCrop(context, minimal, false, &cand)) {
      minimal = full;
      ref = record.reference.value;
      cand = record.candidate.value;
    }
    const double rank_reference =
        InsertionAverageRank(background_reference, ref);
    const double rank_candidate =
        InsertionAverageRank(background_candidate, cand);
    Counterexample counterexample;
    counterexample.kind = "rank";
    counterexample.group = record.group;
    counterexample.level = record.level;
    counterexample.family = kFamilyNames[group.family];
    counterexample.group_seed = group.seed;
    counterexample.target_seed = group.target_seed;
    counterexample.base_kind = kBaseKindNames[group.base_kind];
    counterexample.amplitude = group.amplitude;
    counterexample.alpha_lo = kAlphas[record.level];
    counterexample.alpha_hi = kAlphas[record.level];
    counterexample.crop = minimal;
    counterexample.original_w = group.w;
    counterexample.original_h = group.h;
    counterexample.reference_lo = ref;
    counterexample.reference_hi = ref;
    counterexample.candidate_lo = cand;
    counterexample.candidate_hi = cand;
    counterexample.metric =
        std::fabs(rank_candidate - rank_reference);
    counterexample.threshold = rank_threshold;
    std::ostringstream note;
    note << "reference rank " << std::fixed << std::setprecision(1)
         << rank_reference << ", candidate rank " << rank_candidate
         << " (original " << original_rank_reference << "/"
         << original_rank_candidate << ")";
    counterexample.note = note.str();
    counterexamples.push_back(counterexample);
  }
  return counterexamples;
}

void WriteCounterexampleCsv(const std::string &path,
                            const std::vector<Counterexample> &items) {
  std::ofstream out(path);
  out << "kind,group,level,family,original_w,original_h,crop_x,crop_y,"
         "minimal_w,minimal_h,reference_lo,reference_hi,candidate_lo,"
         "candidate_hi,metric,threshold,group_seed,target_seed,base_kind,"
         "amplitude,alpha_lo,alpha_hi,note\n";
  out << std::setprecision(17);
  for (const Counterexample &item : items) {
    out << item.kind << ',' << item.group << ',' << item.level << ','
        << item.family << ',' << item.original_w << ',' << item.original_h
        << ',' << item.crop.x << ',' << item.crop.y << ',' << item.crop.w << ','
        << item.crop.h << ',' << item.reference_lo << ',' << item.reference_hi
        << ',' << item.candidate_lo << ',' << item.candidate_hi << ','
        << item.metric << ',' << item.threshold << ',' << item.group_seed
        << ',' << item.target_seed << ',' << item.base_kind << ','
        << item.amplitude << ',' << item.alpha_lo << ',' << item.alpha_hi
        << ",\"" << item.note << "\"\n";
  }
}

void WriteCounterexampleMarkdown(const std::string &path,
                                 const std::vector<Counterexample> &items) {
  std::ofstream out(path);
  out << "# Counterexamples\n\n";
  if (items.empty()) {
    out << "No monotonicity or rank-outlier counterexamples were found.\n";
    return;
  }
  out << "All counterexamples are in `counterexamples.csv`. The table below "
         "shows the minimized cases.\n\n";
  out << "| # | kind | group | level | family | original | minimal crop | "
         "metric | threshold |\n";
  out << "|---:|---|---:|---:|---|---|---|---:|---:|\n";
  out << std::setprecision(8);
  for (size_t i = 0; i < items.size(); ++i) {
    const Counterexample &item = items[i];
    out << '|' << i + 1 << '|' << item.kind << '|' << item.group << '|'
        << item.level << '|' << item.family << '|' << item.original_w << 'x'
        << item.original_h << '|' << item.crop.w << 'x' << item.crop.h
        << " @ " << item.crop.x << ',' << item.crop.y << '|' << item.metric
        << '|' << item.threshold << "|\n";
  }
}

void WriteExtremeCsv(const std::string &path,
                     const std::vector<ExtremeResult> &results) {
  std::ofstream out(path);
  out << "width,height,candidate_crash,candidate_signal,reference_crash,"
         "reference_signal\n";
  for (const ExtremeResult &result : results) {
    out << result.w << ',' << result.h << ',' << result.candidate_crash << ','
        << result.candidate_signal << ',' << result.reference_crash << ','
        << result.reference_signal << '\n';
  }
}

void WriteSpearmanCsv(const std::string &path,
                      const std::vector<std::pair<std::string, SpearmanResult>>
                          &rows,
                      const char *metric) {
  std::ofstream out(path);
  out << metric << ",samples,spearman,defined\n";
  out << std::setprecision(17);
  for (const auto &row : rows) {
    out << row.first << ',' << row.second.n << ',' << row.second.rho << ','
        << row.second.defined << '\n';
  }
}

std::string JsonEscape(const std::string &value) {
  std::ostringstream out;
  for (char c : value) {
    switch (c) {
      case '\\':
        out << "\\\\";
        break;
      case '"':
        out << "\\\"";
        break;
      case '\n':
        out << "\\n";
        break;
      case '\r':
        out << "\\r";
        break;
      case '\t':
        out << "\\t";
        break;
      default:
        out << c;
        break;
    }
  }
  return out.str();
}

int Run(int argc, char **argv) {
  if (argc >= 2 && std::strcmp(argv[1], "--extreme-child") == 0) {
    if (argc != 5) return 2;
    return RunExtremeChild(argv[2], std::atoi(argv[3]), std::atoi(argv[4]));
  }

  Config config;
  if (!ParseConfig(argc, argv, &config)) return 2;
  if (config.rank_tolerance < 0.0 || config.rank_tolerance > 1.0) {
    std::cerr << "--rank-tolerance must be in [0, 1]\n";
    return 2;
  }
  const std::string mkdir_command = "mkdir -p '" + config.out_dir + "'";
  if (std::system(mkdir_command.c_str()) != 0) {
    std::cerr << "Cannot create output directory " << config.out_dir << "\n";
    return 2;
  }

  const auto start = std::chrono::steady_clock::now();
  Rng rng(config.seed);
  std::vector<GroupSpec> groups;
  groups.reserve(config.groups);
  for (int i = 0; i < config.groups; ++i) groups.push_back(MakeGroup(i, &rng));

  std::vector<Record> records;
  records.reserve(static_cast<size_t>(config.groups) * kLevels);
  size_t invalid_scores = 0;
  size_t deterministic_mismatches = 0;
  std::cout << "Running " << config.groups << " groups x " << kLevels
            << " levels against the reference implementation...\n";
  for (int group_index = 0; group_index < config.groups; ++group_index) {
    const GroupSpec &group = groups[group_index];
    const RgbImage base = MakeBase(group);
    const LinearImage base_linear = ToLinear(base);
    const RgbImage target = MakeTarget(base, group);
    for (int level = 0; level < kLevels; ++level) {
      const RgbImage distorted = Mix(base, target, kAlphas[level]);
      const LinearImage distorted_linear = ToLinear(distorted);
      Record record;
      record.group = group_index;
      record.level = level;
      record.family = group.family;
      ScoreImages(base_linear, distorted_linear, true, &record.reference);
      ScoreImages(base_linear, distorted_linear, false, &record.candidate);
      ScoreImages(base_linear, distorted_linear, false,
                  &record.candidate_repeat);
      record.deterministic =
          record.candidate.ok && record.candidate_repeat.ok &&
          record.candidate.bits == record.candidate_repeat.bits;
      if (!record.reference.ok || !record.candidate.ok ||
          !record.candidate_repeat.ok) {
        ++invalid_scores;
      }
      if (!record.deterministic) ++deterministic_mismatches;
      records.push_back(record);
    }
    if ((group_index + 1) % 250 == 0 || group_index + 1 == config.groups) {
      std::cout << "  scored " << group_index + 1 << '/' << config.groups
                << " groups\r" << std::flush;
    }
  }
  std::cout << "\n";

  std::vector<size_t> valid;
  valid.reserve(records.size());
  for (size_t i = 0; i < records.size(); ++i) {
    if (records[i].reference.ok && records[i].candidate.ok) valid.push_back(i);
  }
  const SpearmanResult overall = SpearmanFor(records, valid, true);
  std::vector<size_t> distorted;
  distorted.reserve(valid.size());
  for (size_t index : valid) {
    if (records[index].level > 0) distorted.push_back(index);
  }
  const SpearmanResult overall_distorted =
      SpearmanFor(records, distorted, true);
  std::vector<std::pair<std::string, SpearmanResult>> by_level;
  for (int level = 0; level < kLevels; ++level) {
    std::vector<size_t> subset;
    for (size_t index : valid) {
      if (records[index].level == level) subset.push_back(index);
    }
    by_level.emplace_back("level_" + std::to_string(level),
                          SpearmanFor(records, subset, true));
  }
  std::vector<std::pair<std::string, SpearmanResult>> by_family;
  for (int family = 0; family < 6; ++family) {
    std::vector<size_t> subset;
    for (size_t index : valid) {
      if (records[index].family == family) subset.push_back(index);
    }
    by_family.emplace_back(kFamilyNames[family],
                           SpearmanFor(records, subset, true));
  }

  std::vector<double> values_reference;
  std::vector<double> values_candidate;
  values_reference.reserve(valid.size());
  values_candidate.reserve(valid.size());
  for (size_t index : valid) {
    values_reference.push_back(records[index].reference.value);
    values_candidate.push_back(records[index].candidate.value);
  }
  std::vector<size_t> all(valid.size());
  for (size_t i = 0; i < all.size(); ++i) all[i] = i;
  const std::vector<double> ranks_reference =
      AverageRanks(values_reference, all);
  const std::vector<double> ranks_candidate =
      AverageRanks(values_candidate, all);
  const double rank_threshold =
      config.rank_tolerance * std::max(0.0, static_cast<double>(valid.size() - 1));

  std::vector<size_t> rank_failures;
  std::vector<std::pair<double, size_t>> rank_deltas;
  for (size_t i = 0; i < valid.size(); ++i) {
    const double delta = std::fabs(ranks_candidate[i] - ranks_reference[i]);
    rank_deltas.emplace_back(delta, valid[i]);
    if (delta > rank_threshold) rank_failures.push_back(valid[i]);
  }
  if (!overall.defined || overall.rho < config.min_spearman) {
    std::sort(rank_deltas.begin(), rank_deltas.end(),
              [](const std::pair<double, size_t> &a,
                 const std::pair<double, size_t> &b) {
                return a.first > b.first;
              });
    if (rank_failures.empty() && !rank_deltas.empty()) {
      const size_t count =
          std::max<size_t>(1, rank_deltas.size() / 1000);
      for (size_t i = 0; i < count; ++i) {
        rank_failures.push_back(rank_deltas[i].second);
      }
    }
  }
  std::sort(rank_failures.begin(), rank_failures.end());
  rank_failures.erase(std::unique(rank_failures.begin(), rank_failures.end()),
                      rank_failures.end());

  size_t monotonic_violations = 0;
  size_t reference_monotonic_violations = 0;
  size_t candidate_only_monotonic_violations = 0;
  for (size_t index = 1; index < records.size(); ++index) {
    const Record &previous = records[index - 1];
    const Record &current = records[index];
    if (previous.group != current.group || current.level == 0) continue;
    if (!previous.candidate.ok || !current.candidate.ok) continue;
    if (previous.candidate.value - current.candidate.value >
        MonotonicTolerance(previous.candidate.value)) {
      ++monotonic_violations;
      if (previous.reference.value - current.reference.value <=
          MonotonicTolerance(previous.reference.value)) {
        ++candidate_only_monotonic_violations;
      }
    }
    if (previous.reference.value - current.reference.value >
        MonotonicTolerance(previous.reference.value)) {
      ++reference_monotonic_violations;
    }
  }

  std::vector<Counterexample> monotonic_counterexamples;
  std::vector<Counterexample> rank_counterexamples;
  if (config.shrink && monotonic_violations > 0) {
    std::cout << "Shrinking " << monotonic_violations
              << " monotonicity failures...\n";
    monotonic_counterexamples =
        ShrinkMonotonicFailures(config, groups, records);
  }
  if (config.shrink && !rank_failures.empty()) {
    std::cout << "Shrinking " << rank_failures.size()
              << " rank-outlier samples...\n";
    rank_counterexamples =
        ShrinkRankFailures(config, groups, records, valid, rank_threshold,
                           rank_failures);
  }

  std::vector<Counterexample> counterexamples = monotonic_counterexamples;
  counterexamples.insert(counterexamples.end(), rank_counterexamples.begin(),
                         rank_counterexamples.end());
  WriteCounterexampleCsv(JoinPath(config.out_dir, "counterexamples.csv"),
                         counterexamples);
  WriteCounterexampleMarkdown(JoinPath(config.out_dir, "counterexamples.md"),
                              counterexamples);
  for (size_t i = 0;
       i < counterexamples.size() &&
       static_cast<int>(i) < config.max_artifacts;
       ++i) {
    const Counterexample &item = counterexamples[i];
    PairContext context = MakePairContext(groups[item.group], item.level);
    const RgbImage base = CropRgb(context.base, item.crop);
    const RgbImage distorted = CropRgb(context.distorted, item.crop);
    std::ostringstream prefix;
    prefix << "counterexample_" << std::setw(4) << std::setfill('0')
           << i + 1;
    WritePpm(JoinPath(config.out_dir, prefix.str() + "_orig.ppm"), base);
    WritePpm(JoinPath(config.out_dir, prefix.str() + "_dist.ppm"), distorted);
  }
  WriteSpearmanCsv(JoinPath(config.out_dir, "spearman_by_level.csv"),
                   by_level, "level");
  WriteSpearmanCsv(JoinPath(config.out_dir, "spearman_by_family.csv"),
                   by_family, "family");

  const std::array<std::pair<int, int>, 22> extreme_sizes = {{
      {0, 0},    {0, 1},    {1, 0},    {1, 1},    {1, 2},
      {2, 1},    {1, 7},    {7, 1},    {2, 2},    {3, 3},
      {5, 7},    {7, 5},    {1, 4096}, {4096, 1}, {2, 2048},
      {2048, 2}, {1, 8192}, {8192, 1}, {3, 1024}, {1024, 3},
      {16, 16},  {64, 64},
  }};
  std::vector<ExtremeResult> extreme_results;
  for (const auto &size : extreme_sizes) {
    const ExtremeResult candidate =
        RunExtremeCase(argv[0], size.first, size.second, false);
    const ExtremeResult reference =
        RunExtremeCase(argv[0], size.first, size.second, true);
    ExtremeResult result = candidate;
    result.reference_crash = reference.reference_crash;
    result.reference_signal = reference.reference_signal;
    extreme_results.push_back(result);
  }
  size_t candidate_extreme_crashes = 0;
  size_t reference_extreme_crashes = 0;
  for (const ExtremeResult &result : extreme_results) {
    candidate_extreme_crashes += result.candidate_crash ? 1 : 0;
    reference_extreme_crashes += result.reference_crash ? 1 : 0;
  }
  WriteExtremeCsv(JoinPath(config.out_dir, "extreme_sizes.csv"),
                  extreme_results);

  const double elapsed =
      std::chrono::duration<double>(std::chrono::steady_clock::now() - start)
          .count();
  const bool correlation_pass =
      overall.defined && overall.rho >= config.min_spearman &&
      overall_distorted.defined &&
      overall_distorted.rho >= config.min_spearman;
  const bool monotonic_pass = monotonic_violations == 0;
  const bool determinism_pass = deterministic_mismatches == 0;
  const bool extremes_pass = candidate_extreme_crashes == 0;
  const bool validity_pass = invalid_scores == 0;
  const bool counterexamples_pass = counterexamples.empty();
  const bool pass = correlation_pass && monotonic_pass && determinism_pass &&
                    extremes_pass && validity_pass && counterexamples_pass;
  const char *status = pass ? "pass" : "fail";

  std::ofstream json(JoinPath(config.out_dir, "summary.json"));
  json << std::setprecision(17);
  json << "{\n"
       << "  \"status\": \"" << status << "\",\n"
       << "  \"timestamp\": \"" << Timestamp() << "\",\n"
       << "  \"seed\": " << config.seed << ",\n"
       << "  \"groups\": " << config.groups << ",\n"
       << "  \"levels_per_group\": " << kLevels << ",\n"
       << "  \"samples\": " << records.size() << ",\n"
       << "  \"reference\": {\"label\": \""
       << JsonEscape(config.reference_label) << "\", \"sha\": \""
       << JsonEscape(config.reference_sha) << "\"},\n"
       << "  \"candidate\": {\"label\": \""
       << JsonEscape(config.candidate_label) << "\", \"sha\": \""
       << JsonEscape(config.candidate_sha) << "\"},\n"
       << "  \"spearman\": {\"overall\": " << overall.rho
       << ", \"distorted_only\": " << overall_distorted.rho
       << ", \"minimum\": " << config.min_spearman
       << ", \"defined\": " << overall.defined
       << ", \"pass\": " << correlation_pass << "},\n"
       << "  \"monotonicity\": {\"checks\": "
       << static_cast<size_t>(config.groups) * (kLevels - 1)
       << ", \"violations\": " << monotonic_violations
       << ", \"reference_violations\": " << reference_monotonic_violations
       << ", \"candidate_only_violations\": "
       << candidate_only_monotonic_violations
       << ", \"pass\": " << monotonic_pass << "},\n"
       << "  \"determinism\": {\"checks\": " << records.size()
       << ", \"bitwise_mismatches\": " << deterministic_mismatches
       << ", \"pass\": " << determinism_pass << "},\n"
       << "  \"extreme_sizes\": {\"cases\": " << extreme_results.size()
       << ", \"candidate_crashes\": " << candidate_extreme_crashes
       << ", \"reference_crashes\": " << reference_extreme_crashes
       << ", \"pass\": " << extremes_pass << "},\n"
       << "  \"invalid_scores\": " << invalid_scores << ",\n"
       << "  \"rank_outliers\": {\"fraction\": " << config.rank_tolerance
       << ", \"threshold_ranks\": " << rank_threshold
       << ", \"samples\": " << rank_failures.size() << "},\n"
       << "  \"counterexamples\": {\"monotonicity\": "
       << monotonic_counterexamples.size()
       << ", \"rank\": " << rank_counterexamples.size()
       << ", \"total\": " << counterexamples.size() << "},\n"
       << "  \"elapsed_seconds\": " << elapsed << "\n"
       << "}\n";

  std::ofstream markdown(JoinPath(config.out_dir, "summary.md"));
  markdown << std::setprecision(8);
  markdown << "# Differential regression summary\n\n"
           << "- Status: **" << status << "**\n"
           << "- Reference: " << config.reference_label;
  if (!config.reference_sha.empty()) markdown << " (`" << config.reference_sha << "`)";
  markdown << "\n- Candidate: " << config.candidate_label;
  if (!config.candidate_sha.empty()) markdown << " (`" << config.candidate_sha << "`)";
  markdown << "\n- Groups: " << config.groups << " x " << kLevels
           << " levels = " << records.size() << " samples\n"
           << "- Spearman overall: " << overall.rho << " (minimum "
           << config.min_spearman << ")\n"
           << "- Spearman on distorted samples only: "
           << overall_distorted.rho << "\n"
           << "- Monotonicity: " << monotonic_violations << " violations in "
           << static_cast<size_t>(config.groups) * (kLevels - 1)
           << " checks (reference: " << reference_monotonic_violations
           << ", candidate-only: " << candidate_only_monotonic_violations
           << ")\n"
           << "- Determinism: " << deterministic_mismatches
           << " bitwise mismatches in " << records.size() << " checks\n"
           << "- Extreme sizes: " << candidate_extreme_crashes
           << " candidate crashes in " << extreme_results.size() << " cases\n"
           << "- Invalid scores: " << invalid_scores << "\n"
           << "- Rank outliers: " << rank_failures.size() << " (|rank delta| > "
           << rank_threshold << ")\n"
           << "- Shrunk counterexamples: " << counterexamples.size() << "\n"
           << "- Elapsed: " << elapsed << " s\n\n"
           << "## Spearman by severity\n\n"
           << "| severity | samples | rho |\n"
           << "|---:|---:|---:|\n";
  for (const auto &row : by_level) {
    markdown << '|' << row.first << '|' << row.second.n << '|'
             << row.second.rho << "|\n";
  }
  markdown << "\n## Spearman by distortion family\n\n"
           << "| family | samples | rho |\n"
           << "|---|---:|---:|\n";
  for (const auto &row : by_family) {
    markdown << '|' << row.first << '|' << row.second.n << '|'
             << row.second.rho << "|\n";
  }

  std::cout << std::fixed << std::setprecision(6)
            << "Spearman overall: " << overall.rho << " (min "
            << config.min_spearman << ")\n"
            << "Spearman distorted-only: " << overall_distorted.rho << "\n"
            << "Monotonicity violations: " << monotonic_violations << "\n"
            << "Reference monotonicity violations: "
            << reference_monotonic_violations << "\n"
            << "Candidate-only monotonicity violations: "
            << candidate_only_monotonic_violations << "\n"
            << "Determinism mismatches: " << deterministic_mismatches << "\n"
            << "Extreme candidate crashes: " << candidate_extreme_crashes
            << " / " << extreme_results.size() << "\n"
            << "Counterexamples: " << counterexamples.size() << "\n"
            << "Status: " << status << "\n"
            << "Reports: " << config.out_dir << "\n";
  return pass ? 0 : 1;
}

}  // namespace

int main(int argc, char **argv) {
  return Run(argc, argv);
}
