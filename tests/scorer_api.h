#ifndef BUTTERAUGLI_TESTS_SCORER_API_H_
#define BUTTERAUGLI_TESTS_SCORER_API_H_

#ifdef __cplusplus
extern "C" {
#endif

// Scores two linear RGB images. Each plane is width * height floats in
// row-major order. Returns 1 when a finite score was produced, 0 otherwise.
int ButteraugliCandidateScore(const float *r0, const float *g0, const float *b0,
                              const float *r1, const float *g1, const float *b1,
                              int width, int height, double *score);

int ButteraugliReferenceScore(const float *r0, const float *g0, const float *b0,
                              const float *r1, const float *g1, const float *b1,
                              int width, int height, double *score);

#ifdef __cplusplus
}
#endif

#endif  // BUTTERAUGLI_TESTS_SCORER_API_H_
