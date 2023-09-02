#ifndef BICYCLE_5D_H
# define BICYCLE_5D_H

#include <Eigen.h>

# ifdef __cplusplus
extern "C" {
# endif // ifdef __cplusplus


//========================================
// Bicycle model parameters
// X = [x; y; theta; v; delta] : x, y, yaw, linear vel, steering angle
// U = [a; delta_dot] : linear accel and steering rate
//========================================
// struct tiny_Model_Bicycle {
//   float drive_min[2];
//   float drive_max[2];
//   float u_min[2];
//   float u_max[2];
// } tiny_DefaultModel_Bicycle = {{-2, -0.5}, {2, 0.5}, {-4, -0.7}, {4, 0.7}};

//========================================
// Codes generated from julia/bicycle_tvlqr
// Discrete dynamics of bicycle model with predefined model params
//========================================
void tiny_Bicycle5dNonlinearDynamics_Raw(float* xn, const float* x,
                                         const float* u);

void tiny_Bicycle5dNonlinearDynamics(Eigen::MatrixNf* xn, Eigen::MatrixNf* x,
                                     Eigen::MatrixMf* u);

//========================================
// Codes generated from julia/bicycle_tvlqr
// Jacobians of discrete dynamics of bicycle model with predefined model params
//========================================
void tiny_Bicycle5dGetJacobianA_Raw(float* A, const float* x,
                                    const float* u);

void tiny_Bicycle5dGetJacobianB_Raw(float* B, const float* x,
                                    const float* u);

void tiny_Bicycle5dGetJacobians(Eigen::MatrixNf* A, Eigen::MatrixNMf* B, Eigen::MatrixNf* x,
                                Eigen::MatrixMf* u);


# ifdef __cplusplus
}
# endif // ifdef __cplusplus

#endif // ifndef BICYCLE_5D_H