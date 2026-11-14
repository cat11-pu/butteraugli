// Calibration helper: finds which single-pixel position in a hot_pixel
// pattern reproduces a target butteraugli score.
#include <cstdio>
#include <cstdlib>
#include <vector>

#include "butteraugli/butteraugli.h"

int main(int argc, char** argv) {
  if (argc != 4) {
    fprintf(stderr, "usage: %s W H TARGET_SCORE\n", argv[0]);
    return 2;
  }
  const size_t xsize = atoi(argv[1]);
  const size_t ysize = atoi(argv[2]);
  const double target = atof(argv[3]);
  for (size_t py = 0; py < ysize; ++py) {
    for (size_t px = 0; px < xsize; ++px) {
      std::vector<butteraugli::ImageF> black =
          butteraugli::CreatePlanes<float>(xsize, ysize, 3);
      std::vector<butteraugli::ImageF> hot =
          butteraugli::CreatePlanes<float>(xsize, ysize, 3);
      for (int c = 0; c < 3; ++c) {
        for (size_t y = 0; y < ysize; ++y) {
          float* row0 = black[c].Row(y);
          float* row1 = hot[c].Row(y);
          for (size_t x = 0; x < xsize; ++x) {
            row0[x] = 0.0f;
            row1[x] = 0.0f;
          }
        }
        hot[c].Row(py)[px] = 255.0f;
      }
      butteraugli::ImageF diffmap;
      double diffvalue = 0.0;
      butteraugli::ButteraugliInterface(hot, black, 1.0f, diffmap, diffvalue);
      if (fabs(diffvalue - target) < 5e-4) {
        printf("MATCH pos=(%zu,%zu) score=%.7g\n", px, py, diffvalue);
      }
    }
  }
  return 0;
}
