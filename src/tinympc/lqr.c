#include "lqr.h"

enum tiny_ErrorCode tiny_RollOutClosedLoopCost(tiny_ADMMWorkspace* work) {
  tiny_Model* model = work->data->model;
  int N = model[0].nhorizon;
  int adaptive_horizon = work->stgs->adaptive_horizon;
  work->info->obj_val = 0.0;
  
  if (model[0].ltv) {
    return TINY_NOT_SUPPORTED;
  }
  else {
    for (int k = 0; k < N - 1; ++k) {
      // Control input: u = - d - K*x
      slap_Copy(work->soln->U[k], work->soln->d[k]); // u[k] = -d[k]
      slap_MatMulAdd(work->soln->U[k], work->soln->Kinf, work->soln->X[k], -1, -1);  // u[k] -= Kinf * x[k]
      // Next state: x = A*x + B*u + f
      if (adaptive_horizon && k > adaptive_horizon - 1) {        
        tiny_EvalModel(&(work->soln->X[k + 1]), work->soln->X[k], work->soln->U[k], &model[1], 0);
      }
      else {
        tiny_EvalModel(&(work->soln->X[k + 1]), work->soln->X[k], work->soln->U[k], &model[0], 0);
      }
      tiny_AddStageCost(work, k);
      // printf("%f\n", work->info->obj_val);
    } 
    // printf("%f\n", work->info->obj_val);
    tiny_AddTerminalCost(work);       
  }
  return TINY_NO_ERROR;
}

enum tiny_ErrorCode tiny_ForwardPass(tiny_ADMMWorkspace* work) {
  tiny_RollOutClosedLoop(work);
  // tiny_RollOutClosedLoopCost(work);
  return TINY_NO_ERROR;
}

enum tiny_ErrorCode tiny_BackwardPassGrad(tiny_ADMMWorkspace* work) {
  tiny_Model* model = work->data->model;
  int N = model[0].nhorizon;
  // int n = model[0].nstates;
  // int m = model[0].ninputs;

  if (model[0].ltv && model[0].affine) {
    return TINY_NOT_SUPPORTED;
  }
  if (model[0].ltv && !model[0].affine) {
    return TINY_NOT_SUPPORTED;
  }
  if (!model[0].ltv && model[0].affine) {
    return TINY_NOT_SUPPORTED;
  }
  // LTI model
  if (!model[0].ltv && !model[0].affine) {
    for (int k = N - 2; k >= 0; --k) {
      // Qu = B'*p[k+1] + r[k]
      slap_MatMulAtB(work->Qu, model[0].B[0], work->soln->p[k+1]);
      slap_MatrixAddition(work->Qu, work->Qu, work->data->r[k], 1);
      // d = Quu\Qu
      slap_MatMulAB(work->soln->d[k], work->Quu_inv, work->Qu);
      // p[k] .= q[k] + AmBKt*p[k+1] - Kinf'*r[k] + coeff_d2p*d[k]
      slap_Copy(work->soln->p[k], work->data->q[k]);
      slap_MatMulAdd(work->soln->p[k], work->AmBKt, work->soln->p[k+1], 1, 1);
      slap_MatMulAdd(work->soln->p[k], slap_Transpose(work->soln->Kinf), work->data->r[k], -1, 1);
      slap_MatMulAdd(work->soln->p[k], work->coeff_d2p, work->soln->d[k], 1, 1);
    }
  }
  return TINY_NO_ERROR;
}

enum tiny_ErrorCode tiny_SolveLqr(tiny_ADMMWorkspace* work) {

  slap_Copy(work->soln->X[0], work->data->x0);

  tiny_BackwardPassGrad(work);
  tiny_ForwardPass(work);

  // work->info->status_val = TINY_SOLVED;
  
  return TINY_NO_ERROR;
}