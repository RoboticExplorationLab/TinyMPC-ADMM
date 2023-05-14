#ifndef LQR_H
# define LQR_H

# ifdef __cplusplus
extern "C" {
# endif // ifdef __cplusplus


#include "types.h"
#include "cost_lqr.h"
#include "model.h"
#include "auxil.h"

enum tiny_ErrorCode tiny_RollOutClosedLoopCost(tiny_ADMMWorkspace* work);

enum tiny_ErrorCode tiny_ForwardPass(tiny_ADMMWorkspace* work);

enum tiny_ErrorCode tiny_BackwardPassGrad(tiny_ADMMWorkspace* work);

enum tiny_ErrorCode tiny_SolveLqr(tiny_ADMMWorkspace* work);


# ifdef __cplusplus
}
# endif // ifdef __cplusplus

#endif // ifndef LQR_H