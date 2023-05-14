#ifndef ADMM_H
# define ADMM_H

# ifdef __cplusplus
extern "C" {
# endif // ifdef __cplusplus


#include "types.h"
#include "utils.h"
#include "lqr.h"

enum tiny_ErrorCode tiny_SolveADMM(tiny_ADMMWorkspace* work);

enum tiny_ErrorCode UpdatePrimal(tiny_ADMMWorkspace* work);

enum tiny_ErrorCode UpdateSlackDual(tiny_ADMMWorkspace* work);

enum tiny_ErrorCode ComputePrimalResidual(tiny_ADMMWorkspace* work);

enum tiny_ErrorCode ComputeDualResidual(tiny_ADMMWorkspace* work);

enum tiny_ErrorCode ComputeObjectiveValue(tiny_ADMMWorkspace* work);

int CheckTermination(tiny_ADMMWorkspace* work);

enum tiny_ErrorCode UpdateInfo(tiny_ADMMWorkspace* work,
                                int                 iter,
                                int                 compute_objective);

enum tiny_ErrorCode tiny_WarmStartInput(tiny_ADMMWorkspace* work, sfloat* U_data);

# ifdef __cplusplus
}
# endif // ifdef __cplusplus

#endif // ifndef ADMM_H