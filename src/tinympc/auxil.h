#ifndef AUXIL_H
# define AUXIL_H

# ifdef __cplusplus
extern "C" {
# endif // ifdef __cplusplus


#include "types.h"
#include "utils.h"

enum tiny_ErrorCode tiny_InitSettings(tiny_ADMMSettings* stgs);

enum tiny_ErrorCode tiny_SetUnconstrained(tiny_ADMMSettings* stgs);

// enum tiny_ErrorCode tiny_InitData(tiny_ADMMWorkspace* work);

enum tiny_ErrorCode tiny_InitDataQuadCostFromArray(tiny_ADMMWorkspace* work, 
sfloat* Q_data, sfloat* R_data);

enum tiny_ErrorCode tiny_InitDataLinearCostFromArray(tiny_ADMMWorkspace* work, 
Matrix* q, Matrix* r, sfloat* q_data, sfloat* r_data);

enum tiny_ErrorCode tiny_InitSolnTrajFromArray(tiny_ADMMWorkspace* work,
Matrix* X, Matrix* U,
sfloat* X_data, sfloat* U_data);

enum tiny_ErrorCode tiny_InitSolnDualsFromArray(tiny_ADMMWorkspace* work,
Matrix* YX, Matrix* YU,
sfloat* YX_data, sfloat* YU_data, sfloat* YG_data);

enum tiny_ErrorCode tiny_InitSolnGainsFromArray(tiny_ADMMWorkspace* work, 
Matrix* d, Matrix* p, sfloat* d_data, sfloat* p_data, 
sfloat* Kinf_data, sfloat* Pinf_data);

enum tiny_ErrorCode tiny_InitWorkspace(tiny_ADMMWorkspace* work,
                                       tiny_ADMMInfo* info,
                                       tiny_Model* model,
                                       tiny_ADMMData* data,
                                       tiny_ADMMSolution* soln,
                                       tiny_ADMMSettings* stgs);

enum tiny_ErrorCode tiny_InitWorkspaceTempData(tiny_ADMMWorkspace* work, 
Matrix* ZU, Matrix* ZU_new, Matrix* ZX, Matrix* ZX_new, sfloat* temp_data);

// enum tiny_ErrorCode tiny_EvalPrimalCache(tiny_ADMMWorkspace* work);

enum tiny_ErrorCode tiny_LoadPrimalCache(tiny_ADMMWorkspace* work, 
sfloat* Quu_inv_data, sfloat* AmBKt_data, sfloat* coeff_d2p_data);

enum tiny_ErrorCode tiny_ResetInfo(tiny_ADMMWorkspace* work);

# ifdef __cplusplus
}
# endif // ifdef __cplusplus

#endif // ifndef AUXIL_H
