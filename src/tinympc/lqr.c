#include "lqr.h"

enum tiny_ErrorCode tiny_ForwardPass(tiny_AdmmWorkspace* work) {
  tiny_RollOutClosedLoop(work);
  return TINY_NO_ERROR;
}

enum tiny_ErrorCode tiny_BackwardPassGrad(tiny_AdmmWorkspace* work) {
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
      // slap_Copy(work->Qu, work->data->r[k]);
      // slap_MatMulAdd(work->Qu, slap_Transpose(model[0].B[0]), work->soln->p[k+1], 1, 1);
      slap_MatMulAtB(work->Qu, model[0].B[0], work->soln->p[k+1]);
      slap_MatrixAddition(work->Qu, work->Qu, work->data->r[k], 1);
      // d = Quu\Qu
      slap_MatMulAB(work->soln->d[k], work->Quu_inv, work->Qu);
      // slap_MatMulAdd(work->soln->d[k], work->Quu_inv, work->Qu, 1, 0);
      // p[k] .= q[k] + AmBKt*p[k+1] - Kinf'*r[k] + coeff_d2p*d[k]
      slap_Copy(work->soln->p[k], work->data->q[k]);
      slap_MatMulAdd(work->soln->p[k], work->AmBKt, work->soln->p[k+1], 1, 1);
      slap_MatMulAdd(work->soln->p[k], slap_Transpose(work->soln->Kinf), work->data->r[k], -1, 1);
      slap_MatMulAdd(work->soln->p[k], work->coeff_d2p, work->soln->d[k], 1, 1);
    }
  }
  return TINY_NO_ERROR;
}

enum tiny_ErrorCode tiny_SolveLqr(tiny_AdmmWorkspace* work) {

  slap_Copy(work->soln->X[0], work->data->x0);

  tiny_BackwardPassGrad(work);
  
  tiny_ForwardPass(work);

  // work->info->status_val = TINY_SOLVED;
  
  return TINY_NO_ERROR;
}