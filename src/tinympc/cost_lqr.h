#ifndef COST_LQR_H
# define COST_LQR_H

# ifdef __cplusplus
extern "C" {
# endif // ifdef __cplusplus


#include "types.h"

enum tiny_ErrorCode tiny_AddStageCost(tiny_ADMMWorkspace* work, const int k);

enum tiny_ErrorCode tiny_AddTerminalCost(tiny_ADMMWorkspace* work);

enum tiny_ErrorCode tiny_UpdateLinearCost(tiny_ADMMWorkspace* work);

enum tiny_ErrorCode tiny_UpdateConstrainedLinearCost(tiny_ADMMWorkspace* work);

# ifdef __cplusplus
}
# endif // ifdef __cplusplus

#endif // ifndef COST_LQR_H