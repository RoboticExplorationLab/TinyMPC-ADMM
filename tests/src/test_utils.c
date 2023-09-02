#include "test_utils.h"

// This will not work with sub-matrices
float SumOfSquaredError(const float* x, const float* y, const int len) {
  float err = 0;
  for (int i = 0; i < len; ++i) {
    float diff = x[i] - y[i];
    err += diff * diff;
  }
  return sqrt(err);
}

float SumOfSquaredErrorMatrices(const float* x, Matrix* Y, const int num) {
  float err = 0;
  int k = 0;
  for (int i = 0; i < num; ++i) {
    for (int j = 0; j < Y[i].cols*Y[i].rows; ++j) {
      float diff = x[k++] - Y[i].data[j];
      err += diff * diff;
    }
  }
  return sqrt(err);
}