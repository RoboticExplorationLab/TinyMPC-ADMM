#include "cost_lqr.h"

enum tiny_ErrorCode tiny_AddStageCost(tiny_AdmmWorkspace* work, const int k) {
  int n = work->data->model[0].nstates;
  int m = work->data->model[0].ninputs;
  sfloat dx_data[n];
  Matrix dx = slap_MatrixFromArray(n, 1, dx_data);
  // PrintMatrixInfo(dx);
  // PrintMatrixInfo(work->data->Uref[k]);
  MatAdd(dx, work->soln->X[k], work->data->Xref[k], -1);
  work->info->obj_val += 0.5 * slap_QuadraticForm(dx, work->data->Q, dx);
  // PrintMatrixT(dx);
  // printf("%f\n", work->info->obj_val);
  Matrix du = slap_MatrixFromArray(m, 1, dx_data);
  MatAdd(du, work->soln->U[k], work->data->Uref[k], -1);
  work->info->obj_val += 0.5 * slap_QuadraticForm(du, work->data->R, du);
  return TINY_NO_ERROR;
}

enum tiny_ErrorCode tiny_AddTerminalCost(tiny_AdmmWorkspace* work) {
  int n = work->data->model[0].nstates;
  int N = work->data->model[0].nhorizon;
  sfloat dx_data[n];
  Matrix dx = slap_MatrixFromArray(n, 1, dx_data);
  MatAdd(dx, work->soln->X[N - 1], work->data->Xref[N - 1], -1);
  work->info->obj_val += 0.5 * slap_QuadraticForm(dx, work->soln->Pinf, dx);
  return TINY_NO_ERROR;
}

enum tiny_ErrorCode tiny_UpdateLinearCost(tiny_AdmmWorkspace* work) {
  int N = work->data->model[0].nhorizon;
  for (int k = 0; k < N - 1; ++k) {
    // q[k] .= - Q*Xref[k]
    slap_MatMulAB(work->data->q[k], work->data->Q, work->data->Xref[k]);
    MatScale(work->data->q[k], -1);
  
    // r[k] .= - R*Uref[k]
    slap_MatMulAB(work->data->r[k], work->data->R, work->data->Uref[k]);
    MatScale(work->data->r[k], -1);
    // slap_MatrixAddition(work->data->r[k], work->data->r[k], work->data->Uref[k], 
    // work->rho);
  }
    slap_MatMulAB(work->soln->p[N-1], work->soln->Pinf, work->data->Xref[N-1]);
    MatScale(work->soln->p[N-1], -1);
  return TINY_NO_ERROR;
}

enum tiny_ErrorCode tiny_UpdateConstrainedLinearCost(tiny_AdmmWorkspace* work) {
  int N = work->data->model[0].nhorizon;
  for (int k = 0; k < N - 1; ++k) {  
    // slap_MatMulAB(work->data->q[k], work->data->Q, work->data->Xref[k]);
    // slap_MatMulAdd(work->data->q[k], work->data->Q, work->data->Xref[k], 1, 0);
    // MatScale(work->data->q[k], -1);
    MatMulAdd(work->data->q[k], work->data->Q, work->data->Xref[k], -1, 0);

    // r[k] .= -ρ*(z[k]-y[k]) - R*Uref[k]
    // slap_MatMulAB(work->data->r[k], work->data->R, work->data->Uref[k]);
    // MatScale(work->data->r[k], -1);
    MatMulAdd(work->data->r[k], work->data->R, work->data->Uref[k], -1, 0);
    MatAdd(work->data->r[k], work->data->r[k], work->ZU_new[k], -work->rho);
    MatAdd(work->data->r[k], work->data->r[k], work->soln->YU[k], work->rho);
  }
  slap_MatMulAB(work->soln->p[N-1], work->soln->Pinf, work->data->Xref[N-1]);
  MatScale(work->soln->p[N-1], -1);
  return TINY_NO_ERROR;
}