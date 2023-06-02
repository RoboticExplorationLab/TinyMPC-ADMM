#ifndef TEST_UTILS_H
# define TEST_UTILS_H

# ifdef __cplusplus
extern "C" {
# endif // ifdef __cplusplus


#include <math.h>
#include "slap/slap.h"

float SumOfSquaredError(const float* x, const float* y, const int len);

float SumOfSquaredErrorMatrices(const float* x, Matrix* Y, const int num);


#ifdef __cplusplus
}
#endif

#endif // ifndef TEST_UTILS_H