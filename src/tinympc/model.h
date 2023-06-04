#ifndef MODEL_H
# define MODEL_H

#include "types.h"

# ifdef __cplusplus
extern "C" {
# endif // ifdef __cplusplus

enum tiny_ErrorCode tiny_InitModel(tiny_Model* model, const int nstates,
                                   const int ninputs, const int nhorizon,
                                   const int ltv, const int affine, 
                                   const float dt);

// User provides array of slap matrices
enum tiny_ErrorCode tiny_InitModelDataMatrix(tiny_Model* model, 
    Eigen::Matrix12f* A, Eigen::Matrix12x4f* B, Eigen::Vector12f* f);

// User provides matrix as column-major array
enum tiny_ErrorCode tiny_InitModelFromArray(tiny_Model* model, Eigen::Matrix12f* A, 
    Eigen::Matrix12x4f* B, Eigen::Vector12f* f, float* A_array, float* B_array, float* f_array);

// enum tiny_ErrorCode tiny_InitModelMemory(tiny_Model* model, Eigen::MatrixXf* mats,
//     float* data);

// Used after tiny_InitLtvModelMemory and before tiny_UpdateLtvModelJac
// enum tiny_ErrorCode tiny_FillModelMemory(tiny_Model* model, float* A_data, 
// float* B_data, float* f_data);

// enum tiny_ErrorCode tiny_SetModelJacFunc(
//     tiny_Model* model, 
//     void (*get_jacobians)(Eigen::VectorXf*, Eigen::VectorXf*, Eigen::VectorXf*, Eigen::VectorXf*));

// enum tiny_ErrorCode tiny_SetModelNonlFunc(
//     tiny_Model* model, 
//     void (*get_nonl_model)(Eigen::VectorXf*, Eigen::VectorXf*, Eigen::VectorXf*));

// k = 0 to use LTI model
enum tiny_ErrorCode tiny_EvalModel(Eigen::Vector12f* xn, Eigen::Vector12f* x, Eigen::Vector4f* u,
                                   tiny_Model* model, const int k);

enum tiny_ErrorCode tiny_RollOutClosedLoop(tiny_AdmmWorkspace* work);

enum tiny_ErrorCode tiny_RollOutOpenLoop(tiny_AdmmWorkspace* work);

// Update Model Jacobians based on get_jacobians, get_nonl_model, Xref, Uref
// enum tiny_ErrorCode tiny_UpdateModelJac(tiny_AdmmWorkspace* work);

// Update Model Jacobians based on get_jacobians, get_nonl_model, user-input X, U traj
// enum tiny_ErrorCode tiny_UpdateModelJacAbout(tiny_AdmmWorkspace* work,
//                                              Eigen::VectorXf* X, Eigen::VectorXf* U);


# ifdef __cplusplus
}
# endif // ifdef __cplusplus

#endif // ifndef MODEL_H