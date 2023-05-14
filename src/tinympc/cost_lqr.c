#include "cost_lqr.h"

enum tiny_ErrorCode tiny_AddStageCost(tiny_ADMMWorkspace* work, const int k) {
  int n = work->data->model[0].nstates;
  int m = work->data->model[0].ninputs;
  sfloat dx_data[n];
  Matrix dx = slap_MatrixFromArray(n, 1, dx_data);
  slap_MatrixAddition(dx, work->soln->X[k], work->data->X_ref[k], -1);
  work->info->obj_val += 0.5 * slap_QuadraticForm(dx, work->data->Q, dx);
  // PrintMatrixT(dx);
  // printf("%f\n", work->info->obj_val);
  Matrix du = slap_MatrixFromArray(m, 1, dx_data);
  slap_MatrixAddition(du, work->soln->U[k], work->data->U_ref[k], -1);
  work->info->obj_val += 0.5 * slap_QuadraticForm(du, work->data->R, du);

  // PrintMatrixT(du);
  // PrintMatrix(work->data->R);
  // printf("%f\n", work->info->obj_val);
  return TINY_NO_ERROR;
}

enum tiny_ErrorCode tiny_AddTerminalCost(tiny_ADMMWorkspace* work) {
  int n = work->data->model[0].nstates;
  int N = work->data->model[0].nhorizon;
  sfloat dx_data[n];
  Matrix dx = slap_MatrixFromArray(n, 1, dx_data);
  slap_MatrixAddition(dx, work->soln->X[N - 1], work->data->X_ref[N - 1], -1);
  work->info->obj_val += 0.5 * slap_QuadraticForm(dx, work->soln->Pinf, dx);
  return TINY_NO_ERROR;
}

enum tiny_ErrorCode tiny_UpdateLinearCost(tiny_ADMMWorkspace* work) {
  int N = work->data->model[0].nhorizon;
  for (int k = 0; k < N - 1; ++k) {
    slap_MatMulAB(work->data->q[k], work->data->Q, work->data->X_ref[k]);
    slap_ScaleByConst(work->data->q[k], -1);
  
    // r[k] .= - R*Uref[k]
    slap_MatMulAB(work->data->r[k], work->data->R, work->data->U_ref[k]);
    slap_ScaleByConst(work->data->r[k], -1);
  }
    slap_MatMulAB(work->soln->p[N-1], work->soln->Pinf, work->data->X_ref[N-1]);
    slap_ScaleByConst(work->soln->p[N-1], -1);
  return TINY_NO_ERROR;
}

enum tiny_ErrorCode tiny_UpdateConstrainedLinearCost(tiny_ADMMWorkspace* work) {
  int N = work->data->model[0].nhorizon;
  for (int k = 0; k < N - 1; ++k) {
    slap_MatMulAB(work->data->q[k], work->data->Q, work->data->X_ref[k]);
    slap_ScaleByConst(work->data->q[k], -1);
  
    // r[k] .= -ρ*(z[k]-y[k]) - R*Uref[k]
    slap_MatMulAB(work->data->r[k], work->data->R, work->data->U_ref[k]);
    slap_ScaleByConst(work->data->r[k], -1);
    slap_MatrixAddition(work->data->r[k], work->data->r[k], work->ZU_new[k], -work->rho);
    slap_MatrixAddition(work->data->r[k], work->data->r[k], work->soln->YU[k], work->rho);
  }
    slap_MatMulAB(work->soln->p[N-1], work->soln->Pinf, work->data->X_ref[N-1]);
    slap_ScaleByConst(work->soln->p[N-1], -1);
  return TINY_NO_ERROR;
}