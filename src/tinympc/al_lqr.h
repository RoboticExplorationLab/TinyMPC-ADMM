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

enum tiny_ErrorCode tiny_ConstrainedForwardPass(tiny_AdmmWorkspace* work);

enum tiny_ErrorCode tiny_ConstrainedBackwardPass(tiny_AdmmWorkspace* work);

enum tiny_ErrorCode tiny_SolveAlLqr(tiny_AdmmWorkspace* work);

int tiny_CheckRiccati(tiny_AdmmWorkspace* work);

int tiny_CheckAl(tiny_AdmmWorkspace* work);

enum tiny_ErrorCode tiny_ResetWorkspace(tiny_AdmmWorkspace* work);

enum tiny_ErrorCode tiny_WarmStartInput(tiny_AdmmWorkspace* work, sfloat* U_data);

enum tiny_ErrorCode tiny_UpdateDuals(tiny_AdmmWorkspace* work);

enum tiny_ErrorCode tiny_UpdatePenalty(tiny_AdmmWorkspace* work);

# ifdef __cplusplus
}
# endif // ifdef __cplusplus

#endif // ifndef AL_LQR_H