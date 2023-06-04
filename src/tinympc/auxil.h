#ifndef AUXIL_H
# define AUXIL_H

#include <Eigen.h>
#include "types.h"
#include "utils.h"

# ifdef __cplusplus
extern "C" {
# endif // ifdef __cplusplus

enum tiny_ErrorCode tiny_InitSettings(tiny_AdmmSettings* stgs);

enum tiny_ErrorCode tiny_SetUnconstrained(tiny_AdmmSettings* stgs);

// enum tiny_ErrorCode tiny_InitData(tiny_AdmmWorkspace* work);

enum tiny_ErrorCode tiny_InitDataQuadCostFromArray(tiny_AdmmWorkspace* work, 
float* Q_data, float* R_data);

enum tiny_ErrorCode tiny_InitDataLinearCostFromArray(tiny_AdmmWorkspace* work, 
Eigen::VectorNf* q, Eigen::VectorMf* r, Eigen::VectorMf* r_tilde, float* q_data, float* r_data, float* r_tilde_data);

enum tiny_ErrorCode tiny_InitSolnTrajFromArray(tiny_AdmmWorkspace* work,
Eigen::VectorNf* X, Eigen::VectorMf* U,
float* X_data, float* U_data);

enum tiny_ErrorCode tiny_InitSolnDualsFromArray(tiny_AdmmWorkspace* work,
Eigen::VectorNf* YX, Eigen::VectorMf* YU,
float* YX_data, float* YU_data, float* YG_data);

enum tiny_ErrorCode tiny_InitSolnGainsFromArray(tiny_AdmmWorkspace* work, 
Eigen::VectorMf* d, Eigen::VectorNf* p, float* d_data, float* p_data, 
float* Kinf_data, float* Pinf_data);

enum tiny_ErrorCode tiny_InitWorkspace(tiny_AdmmWorkspace* work,
                                       tiny_AdmmInfo* info,
                                       tiny_Model* model,
                                       tiny_AdmmData* data,
                                       tiny_AdmmSolution* soln,
                                       tiny_AdmmSettings* stgs);

enum tiny_ErrorCode tiny_InitWorkspaceTempData(tiny_AdmmWorkspace* work, 
Eigen::VectorMf* ZU, Eigen::VectorMf* ZU_new, Eigen::VectorNf* ZX, Eigen::VectorNf* ZX_new, float* temp_data);

// enum tiny_ErrorCode tiny_EvalPrimalCache(tiny_AdmmWorkspace* work);

enum tiny_ErrorCode tiny_InitPrimalCache(tiny_AdmmWorkspace* work, 
float* Quu_inv_data, float* AmBKt_data, float* coeff_d2p_data);

enum tiny_ErrorCode tiny_ResetInfo(tiny_AdmmWorkspace* work);

enum tiny_ErrorCode tiny_SetStateReference(tiny_AdmmWorkspace* work, Eigen::MatrixXf* Xref, 
float* Xref_data);

enum tiny_ErrorCode tiny_SetInputReference(tiny_AdmmWorkspace* work, Eigen::MatrixXf* Uref,
float* Uref_data);

enum tiny_ErrorCode tiny_SetReference(tiny_AdmmWorkspace* work, Eigen::MatrixXf* Xref, 
Eigen::MatrixXf* Uref, float* Xref_data, float* Uref_data);

enum tiny_ErrorCode tiny_SetGoalReference(tiny_AdmmWorkspace* work, Eigen::VectorNf* Xref,
Eigen::VectorMf* Uref, float* xg_data, float* ug_data);

enum tiny_ErrorCode tiny_SetInitialState(tiny_AdmmWorkspace* work, float* x0_data);

# ifdef __cplusplus
}
# endif // ifdef __cplusplus

#endif // ifndef AUXIL_H
