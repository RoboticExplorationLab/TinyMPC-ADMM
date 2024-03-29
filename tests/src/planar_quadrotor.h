#ifndef PLANAR_QUADROTOR_H
# define PLANAR_QUADROTOR_H

# ifdef __cplusplus
extern "C" {
# endif // ifdef __cplusplus


#include <math.h>

#include "slap/slap.h"

//========================================
// Planar quadrotor model parameters
//========================================

//========================================
// Codes generated from julia/planar_quad_gen
// Discrete dynamics of planar quadrotor
//========================================
void tiny_PQuadNonlinearDynamics_Raw(float* xn, const float* x,
                                     const float* u);

void tiny_PQuadNonlinearDynamics(Matrix* xn, const Matrix x, const Matrix u);

//========================================
// Codes generated from julia/planar_quad_gen
// Jacobians of discrete dynamics of planar quadrotor
//========================================
void tiny_PQuadGetJacobianA_Raw(float* A, const float* x, const float* u);

void tiny_PQuadGetJacobianB_Raw(float* B, const float* x, const float* u);

void tiny_PQuadGetJacobians(Matrix* A, Matrix* B, const Matrix x,
                            const Matrix u);


# ifdef __cplusplus
}
# endif // ifdef __cplusplus

#endif // ifndef PLANAR_QUADROTOR_H