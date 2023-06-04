#include "cost_lqr.h"

# ifdef __cplusplus
extern "C" {
# endif // ifdef __cplusplus

enum tiny_ErrorCode tiny_AddStageCost(tiny_AdmmWorkspace* work, const int k) {
  // int n = work->data->model[0].nstates;
  // int m = work->data->model[0].ninputs;
  // float dx_data[n];
  // Matrix dx = slap_MatrixFromArray(n, 1, dx_data);
  // MatAdd(dx, work->soln->X[k], work->data->Xref[k], -1);
  // work->info->obj_val += 0.5 * slap_QuadraticForm(dx, work->data->Q, dx);
  // Matrix du = slap_MatrixFromArray(m, 1, dx_data);
  // MatAdd(du, work->soln->U[k], work->data->Uref[k], -1);
  // work->info->obj_val += 0.5 * slap_QuadraticForm(du, work->data->R, du);
  work->info->obj_val = (0.5 * (work->soln->X[k] - work->data->Xref[k]).transpose() * 
                         work->data->Q * (work->soln->X[k] - work->data->Xref[k]) +
                         0.5 * (work->soln->U[k] - work->data->Uref[k]).transpose() * 
                         work->data->R * (work->soln->U[k] - work->data->Uref[k]))(0);
  return TINY_NO_ERROR;
}

enum tiny_ErrorCode tiny_AddTerminalCost(tiny_AdmmWorkspace* work) {
  // int n = work->data->model[0].nstates;
  int N = work->data->model[0].nhorizon;
  // float dx_data[n];
  // Matrix dx = slap_MatrixFromArray(n, 1, dx_data);
  // MatAdd(dx, work->soln->X[N - 1], work->data->Xref[N - 1], -1);
  work->info->obj_val += 0.5 * (work->soln->X[N-1] - work->data->Xref[N-1]).transpose() * 
                         work->soln->Pinf * (work->soln->X[N-1] - work->data->Xref[N-1]) ;
  return TINY_NO_ERROR;
}

enum tiny_ErrorCode tiny_UpdateLinearCost(tiny_AdmmWorkspace* work) {
  int N = work->data->model[0].nhorizon;
  for (int k = 0; k < N - 1; ++k) {
    /* Compute q[k] = -Q*Xref[k] */  
    // MatMulAdd(work->data->q[k], work->data->Q, work->data->Xref[k], -1, 0);
    work->data->q[k] = -work->data->Q * work->data->Xref[k];

    /* Compute r[k] = -R*Uref[k] */ 
    // MatMulAdd(work->data->r[k], work->data->R, work->data->Uref[k], -1, 0);
    work->data->r[k] = -work->data->R * work->data->Uref[k];
  }
  /* Compute q[N-1] = -Pinf*Xref[N-1] */ 
  // MatMulAdd(work->soln->p[N-1], work->soln->Pinf, work->data->Xref[N-1], -1, 0);
  work->soln->p[N-1] = -work->soln->Pinf * work->data->Xref[N-1];
  return TINY_NO_ERROR;
}

enum tiny_ErrorCode tiny_UpdateConstrainedLinearCost(tiny_AdmmWorkspace* work) {
  int N = work->data->model[0].nhorizon;
  for (int k = 0; k < N - 1; ++k) {
    /* Compute r_tilde[k] = r[k] - ρ*(z[k]-y[k]) */ 
    // MatAdd(work->data->r_tilde[k], work->ZU_new[k], work->soln->YU[k], -1);
    // MatAdd(work->data->r_tilde[k], work->data->r[k], work->data->r_tilde[k], -work->rho);
    work->data->r_tilde[k] = work->data->r[k] - work->rho * (work->ZU_new[k] - work->soln->YU[k]);
  }
  return TINY_NO_ERROR;
}

# ifdef __cplusplus
}
# endif // ifdef __cplusplus