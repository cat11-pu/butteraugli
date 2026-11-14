// Boundary / robustness scanner for butteraugli.
//
// Runs a single synthetic-image comparison per invocation so that a crash
// on one case does not abort the whole scan. The driver script
// (repro/run_scan.py) assembles the CSV rows.
//
// Usage:
//   boundary_scan --width W --height H --pattern_a PA --pattern_b PB
//
// Prints one line of machine-readable output:
//   ok,score,diffmap_width,diffmap_height,finite,nonnegative,max_diffmap
//
// Patterns (linear-light float values, nominal range [0, 255]):
//   0 black          all 0
//   1 white          all 255
//   2 checker        255 / 0 alternating, 255 at (0,0)
//   3 hot_pixel      single 255 pixel at (0,0), rest 0
//   4 ramp           horizontal gradient 0 .. 255
//   5 near_white     all 254
//   6 fltmax         all FLT_MAX
//   7 posinf         all +inf
//   8 nan            all NaN

#include <cfloat>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

#include "butteraugli/butteraugli.h"

namespace {

int ParseArg(const char* name, int argc, char** argv) {
  const std::string prefix = std::string("--") + name + "=";
  for (int i = 1; i < argc; ++i) {
    if (prefix.compare(0, prefix.size(), argv[i], prefix.size()) == 0) {
      return atoi(argv[i] + prefix.size());
    }
  }
  fprintf(stderr, "missing argument %s\n", name);
  exit(2);
}

float PatternValue(int pattern, size_t x, size_t y, size_t xsize,
                   size_t ysize) {
  switch (pattern) {
    case 0:  // black
      return 0.0f;
    case 1:  // white
      return 255.0f;
    case 2:  // checker, 255 at (0,0)
      return ((x + y) & 1) ? 0.0f : 255.0f;
    case 3:  // hot pixel at (0,0)
      return (x == 0 && y == 0) ? 255.0f : 0.0f;
    case 4:  // horizontal ramp 0 -> 255
      return xsize > 1 ? 255.0f * x / (xsize - 1) : 0.0f;
    case 5:  // near white
      return 254.0f;
    case 6:  // FLT_MAX
      return FLT_MAX;
    case 7:  // +inf
      return HUGE_VALF;
    case 8:  // NaN
      return nanf("");
    default:
      fprintf(stderr, "unknown pattern %d\n", pattern);
      exit(2);
  }
}

std::vector<butteraugli::ImageF> MakeImage(int pattern, size_t xsize,
                                           size_t ysize) {
  std::vector<butteraugli::ImageF> rgb =
      butteraugli::CreatePlanes<float>(xsize, ysize, 3);
  for (int c = 0; c < 3; ++c) {
    for (size_t y = 0; y < ysize; ++y) {
      float* const row = rgb[c].Row(y);
      for (size_t x = 0; x < xsize; ++x) {
        row[x] = PatternValue(pattern, x, y, xsize, ysize);
      }
    }
  }
  return rgb;
}

}  // namespace

int main(int argc, char** argv) {
  const int width = ParseArg("width", argc, argv);
  const int height = ParseArg("height", argc, argv);
  const int pattern_a = ParseArg("pattern_a", argc, argv);
  const int pattern_b = ParseArg("pattern_b", argc, argv);

  std::vector<butteraugli::ImageF> rgb0, rgb1;
  if (width > 0 && height > 0) {
    rgb0 = MakeImage(pattern_a, width, height);
    rgb1 = MakeImage(pattern_b, width, height);
  } else {
    // Zero-size case: exercise the interface with empty planes.
    rgb0 = butteraugli::CreatePlanes<float>(0, 0, 3);
    rgb1 = butteraugli::CreatePlanes<float>(0, 0, 3);
  }

  butteraugli::ImageF diffmap;
  double diffvalue = 0.0;
  const bool ok = butteraugli::ButteraugliInterface(rgb0, rgb1, 1.0f, diffmap,
                                                    diffvalue);
  if (!ok) {
    // Interface rejected the input (e.g. zero size); report a finite zero.
    printf("0,0,0,0,1,1,0\n");
    return 0;
  }

  bool finite = std::isfinite(diffvalue);
  bool nonnegative = diffvalue >= 0.0;
  float max_diffmap = 0.0f;
  for (size_t y = 0; y < diffmap.ysize(); ++y) {
    const float* const row = diffmap.Row(y);
    for (size_t x = 0; x < diffmap.xsize(); ++x) {
      if (!std::isfinite(row[x])) finite = false;
      if (row[x] < 0.0f) nonnegative = false;
      if (row[x] > max_diffmap) max_diffmap = row[x];
    }
  }
  printf("1,%.9g,%zu,%zu,%d,%d,%.9g\n", diffvalue, diffmap.xsize(),
         diffmap.ysize(), finite ? 1 : 0, nonnegative ? 1 : 0,
         static_cast<double>(max_diffmap));
  return 0;
}
