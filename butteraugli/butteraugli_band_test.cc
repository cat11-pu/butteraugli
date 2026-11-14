#include "butteraugli/butteraugli.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <utility>
#include <vector>

namespace {

constexpr int kWidth = 72;
constexpr int kHeight = 72;
constexpr double kPi = 3.14159265358979323846;

std::vector<butteraugli::ImageF> MakeBase(int variant) {
  std::vector<butteraugli::ImageF> image =
      butteraugli::CreatePlanes<float>(kWidth, kHeight, 3);
  for (int y = 0; y < kHeight; ++y) {
    for (int x = 0; x < kWidth; ++x) {
      const double u = static_cast<double>(x) / kWidth;
      const double v = static_cast<double>(y) / kHeight;
      const double diagonal =
          std::sin(2.0 * kPi * (7.0 * u + 5.0 * v + variant * 0.11));
      const double mid = std::cos(2.0 * kPi * (23.0 * u - 19.0 * v));
      const double fine =
          std::sin(2.0 * kPi * (41.0 * u + 37.0 * v + variant * 0.07));
      image[0].Row(y)[x] =
          static_cast<float>(128.0 + 22.0 * diagonal + 11.0 * mid);
      image[1].Row(y)[x] =
          static_cast<float>(132.0 + 18.0 * diagonal + 9.0 * fine);
      image[2].Row(y)[x] =
          static_cast<float>(120.0 + 13.0 * mid + 7.0 * fine);
    }
  }
  return image;
}

float ClampPixel(double value) {
  return static_cast<float>(std::max(0.0, std::min(255.0, value)));
}

void AddPattern(std::vector<butteraugli::ImageF>* image, int type,
                double amplitude, double phase) {
  for (int y = 0; y < kHeight; ++y) {
    for (int x = 0; x < kWidth; ++x) {
      const double u = static_cast<double>(x) / kWidth;
      const double v = static_cast<double>(y) / kHeight;
      double delta = 0.0;
      switch (type) {
        case 0:
          delta = std::sin(2.0 * kPi * (1.0 * u + 0.5 * v + phase));
          break;
        case 1:
          delta = std::sin(2.0 * kPi * (5.0 * u - 4.0 * v + phase));
          break;
        case 2:
          delta =
              ((x + y) & 1) ? std::sin(phase) : -std::sin(phase);
          break;
        case 3:
          delta = std::sin(2.0 * kPi * (2.0 * u + phase)) *
                  std::cos(2.0 * kPi * (9.0 * v - phase));
          break;
        default:
          break;
      }
      for (int c = 0; c < 3; ++c) {
        const double channel_scale = c == 0 ? 1.0 : (c == 1 ? 0.7 : 1.3);
        (*image)[c].Row(y)[x] =
            ClampPixel((*image)[c].Row(y)[x] + channel_scale * amplitude *
                           delta);
      }
    }
  }
}

std::vector<double> AverageRanks(const std::vector<double>& values) {
  std::vector<int> order(values.size());
  for (size_t i = 0; i < order.size(); ++i) {
    order[i] = static_cast<int>(i);
  }
  std::sort(order.begin(), order.end(), [&](int a, int b) {
    if (values[a] != values[b]) return values[a] < values[b];
    return a < b;
  });
  std::vector<double> ranks(values.size());
  for (size_t begin = 0; begin < order.size();) {
    size_t end = begin + 1;
    while (end < order.size() &&
           values[order[end]] == values[order[begin]]) {
      ++end;
    }
    const double rank = 0.5 * (begin + end - 1);
    for (size_t i = begin; i < end; ++i) {
      ranks[order[i]] = rank;
    }
    begin = end;
  }
  return ranks;
}

double Pearson(const std::vector<double>& a,
               const std::vector<double>& b) {
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
  return covariance / std::sqrt(variance_a * variance_b);
}

double Spearman(const std::vector<double>& a,
                const std::vector<double>& b) {
  return Pearson(AverageRanks(a), AverageRanks(b));
}

bool SameBits(double a, double b) {
  return std::memcmp(&a, &b, sizeof(double)) == 0;
}

}  // namespace

int main() {
  struct Case {
    std::vector<butteraugli::ImageF> image0;
    std::vector<butteraugli::ImageF> image1;
    const char* kind;
    double amplitude;
  };
  std::vector<Case> cases;
  const double amplitudes[] = {1.0, 2.0, 4.0, 8.0, 16.0};
  const char* kinds[] = {"low", "mid", "high", "mixed"};
  for (int variant = 0; variant < 3; ++variant) {
    const std::vector<butteraugli::ImageF> base = MakeBase(variant);
    for (int type = 0; type < 4; ++type) {
      for (double amplitude : amplitudes) {
        Case sample;
        sample.image0 = butteraugli::CopyPlanes(base);
        sample.image1 = butteraugli::CopyPlanes(base);
        AddPattern(&sample.image1, type, amplitude, 0.17 + 0.09 * variant);
        sample.kind = kinds[type];
        sample.amplitude = amplitude;
        cases.push_back(std::move(sample));
      }
    }
  }

  std::vector<double> legacy_scores;
  std::vector<double> totals;
  std::vector<double> low_scores;
  std::vector<double> mid_scores;
  std::vector<double> high_scores;
  double max_sum_abs_error = 0.0;
  double max_sum_rel_error = 0.0;
  double max_total_rel_delta = 0.0;
  bool repeat_exact = true;

  std::printf("kind,amplitude,legacy,low,mid,high,total,share_low,share_mid,"
              "share_high\n");
  for (const Case& sample : cases) {
    butteraugli::ImageF diffmap;
    double legacy = 0.0;
    if (!butteraugli::ButteraugliInterface(
            sample.image0, sample.image1, 1.0, diffmap, legacy)) {
      std::fprintf(stderr, "ButteraugliInterface failed\n");
      return 1;
    }

    butteraugli::ButteraugliBandScores first;
    butteraugli::ButteraugliBandScores second;
    if (!butteraugli::ButteraugliBandInterface(
            sample.image0, sample.image1, 1.0, first) ||
        !butteraugli::ButteraugliBandInterface(
            sample.image0, sample.image1, 1.0, second)) {
      std::fprintf(stderr, "ButteraugliBandInterface failed\n");
      return 1;
    }
    repeat_exact =
        repeat_exact && SameBits(first.low, second.low) &&
        SameBits(first.mid, second.mid) &&
        SameBits(first.high, second.high) &&
        SameBits(first.total, second.total);

    const double sum = first.low + first.mid + first.high;
    const double sum_abs_error = std::fabs(sum - first.total);
    const double sum_rel_error =
        first.total > 0.0 ? sum_abs_error / first.total : sum_abs_error;
    const double total_rel_delta =
        legacy > 0.0 ? std::fabs(first.total - legacy) / legacy : 0.0;
    max_sum_abs_error = std::max(max_sum_abs_error, sum_abs_error);
    max_sum_rel_error = std::max(max_sum_rel_error, sum_rel_error);
    max_total_rel_delta = std::max(max_total_rel_delta, total_rel_delta);

    const double share_sum =
        first.low + first.mid + first.high > 0.0 ? sum : 1.0;
    std::printf("%s,%.0f,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,"
                "%.17g\n",
                sample.kind, sample.amplitude, legacy, first.low, first.mid,
                first.high, first.total, first.low / share_sum,
                first.mid / share_sum, first.high / share_sum);

    legacy_scores.push_back(legacy);
    totals.push_back(first.total);
    low_scores.push_back(first.low);
    mid_scores.push_back(first.mid);
    high_scores.push_back(first.high);
  }

  const double spearman = Spearman(legacy_scores, totals);
  const double pearson = Pearson(legacy_scores, totals);
  std::printf(
      "\nsummary cases=%zu spearman=%.17g pearson=%.17g "
      "max_sum_abs_error=%.17g max_sum_rel_error=%.17g "
      "max_total_rel_delta=%.17g repeat_exact=%s\n",
      cases.size(), spearman, pearson, max_sum_abs_error, max_sum_rel_error,
      max_total_rel_delta, repeat_exact ? "true" : "false");

  std::vector<butteraugli::ImageF> small0 =
      butteraugli::CreatePlanes<float>(5, 5, 3);
  std::vector<butteraugli::ImageF> small1 =
      butteraugli::CreatePlanes<float>(5, 5, 3);
  for (int c = 0; c < 3; ++c) {
    for (int y = 0; y < 5; ++y) {
      for (int x = 0; x < 5; ++x) {
        small0[c].Row(y)[x] = static_cast<float>(80 + 13 * x + 7 * y + c);
        small1[c].Row(y)[x] = static_cast<float>(
            small0[c].Row(y)[x] + ((x + y) & 1 ? 9 : -6));
      }
    }
  }
  butteraugli::ButteraugliBandScores small_first;
  butteraugli::ButteraugliBandScores small_second;
  if (!butteraugli::ButteraugliBandInterface(
          small0, small1, 1.0, small_first) ||
      !butteraugli::ButteraugliBandInterface(
          small0, small1, 1.0, small_second)) {
    std::fprintf(stderr, "small-image band interface failed\n");
    return 1;
  }
  const bool small_exact =
      SameBits(small_first.low, small_second.low) &&
      SameBits(small_first.mid, small_second.mid) &&
      SameBits(small_first.high, small_second.high) &&
      SameBits(small_first.total, small_second.total);
  const double small_sum =
      small_first.low + small_first.mid + small_first.high;
  const double small_error =
      small_first.total > 0.0
          ? std::fabs(small_sum - small_first.total) / small_first.total
          : std::fabs(small_sum - small_first.total);
  std::printf("small_image sum_rel_error=%.17g repeat_exact=%s\n",
              small_error, small_exact ? "true" : "false");

  if (max_sum_rel_error > 0.01 || spearman < 0.98 ||
      max_total_rel_delta > 0.01 || !repeat_exact ||
      small_error > 0.01 || !small_exact) {
    std::fprintf(stderr, "band score requirements not met\n");
    return 1;
  }
  return 0;
}
