#ifndef AL_LQR_H
# define AL_LQR_H

# ifdef __cplusplus
extern "C" {
# endif // ifdef __cplusplus


#include "types.h"
#include "model.h"
#include "cost_lqr.h"
#include "constraint_linear.h"
#include "auxil.h"
#include "lqr.h"

enum tiny_ErrorCode tiny_ConstrainedForwardPass(tiny_ADMMWorkspace* work);

enum tiny_ErrorCode tiny_ConstrainedBackwardPass(tiny_ADMMWorkspace* work);

enum tiny_ErrorCode tiny_SolveAlLqr(tiny_ADMMWorkspace* work);

int tiny_CheckRiccati(tiny_ADMMWorkspace* work);

int tiny_CheckAl(tiny_ADMMWorkspace* work);

enum tiny_ErrorCode tiny_ResetWorkspace(tiny_ADMMWorkspace* work);

enum tiny_ErrorCode tiny_WarmStartInput(tiny_ADMMWorkspace* work, sfloat* U_data);

enum tiny_ErrorCode tiny_UpdateDuals(tiny_ADMMWorkspace* work);

enum tiny_ErrorCode tiny_UpdatePenalty(tiny_ADMMWorkspace* work);

# ifdef __cplusplus
}
# endif // ifdef __cplusplus

#endif // ifndef AL_LQR_H