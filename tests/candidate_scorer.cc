#include "scorer_api.h"

#include <cmath>
#include <cstring>
#include <vector>

#include "butteraugli/butteraugli.h"

namespace {

using butteraugli::ButteraugliInterface;
using butteraugli::ImageF;

std::vector<ImageF> PlanesFromInterleaved(const float *const plane[3],
                                          int width, int height) {
  std::vector<ImageF> planes;
  planes.reserve(3);
  for (int c = 0; c < 3; ++c) {
    planes.emplace_back(static_cast<size_t>(width),
                        static_cast<size_t>(height));
    for (int y = 0; y < height; ++y) {
      if (width > 0) {
        std::memcpy(planes.back().Row(y),
                    plane[c] + static_cast<size_t>(y) * width,
                    static_cast<size_t>(width) * sizeof(float));
      }
    }
  }
  return planes;
}

}  // namespace

extern "C" int ButteraugliCandidateScore(
    const float *r0, const float *g0, const float *b0,
    const float *r1, const float *g1, const float *b1,
    int width, int height, double *score) {
  if (width < 0 || height < 0 || score == nullptr) return 0;
  if (width > 0 && height > 0 &&
      (r0 == nullptr || g0 == nullptr || b0 == nullptr ||
       r1 == nullptr || g1 == nullptr || b1 == nullptr)) {
    return 0;
  }
  const float *plane0[3] = {r0, g0, b0};
  const float *plane1[3] = {r1, g1, b1};
  std::vector<ImageF> image0 = PlanesFromInterleaved(plane0, width, height);
  std::vector<ImageF> image1 = PlanesFromInterleaved(plane1, width, height);
  ImageF diffmap;
  double value = 0.0;
  const bool ok = ButteraugliInterface(image0, image1, 1.0, diffmap, value);
  if (!ok || !std::isfinite(value)) return 0;
  *score = value;
  return 1;
}
