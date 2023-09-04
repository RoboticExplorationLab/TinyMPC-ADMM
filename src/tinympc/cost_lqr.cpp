#include "cost_lqr.h"

# ifdef __cplusplus
extern "C" {
# endif // ifdef __cplusplus

enum tiny_ErrorCode tiny_AddStageCost(tiny_AdmmWorkspace* work, const int k) {
  work->info->obj_val += (0.5 * (work->soln->X[k] - work->data->Xref[k]).transpose() * 
                         (*(work->data->Q)) * (work->soln->X[k] - work->data->Xref[k]) +
                         0.5 * (work->soln->U[k] - work->data->Uref[k]).transpose() * 
                         (*(work->data->R)) * (work->soln->U[k] - work->data->Uref[k])).value();
  return TINY_NO_ERROR;
}

enum tiny_ErrorCode tiny_AddTerminalCost(tiny_AdmmWorkspace* work) {
  int N = work->data->model[0].nhorizon;
  work->info->obj_val += 0.5 * (work->soln->X[N-1] - work->data->Xref[N-1]).transpose() * 
                         (*(work->soln->Pinf)) * (work->soln->X[N-1] - work->data->Xref[N-1]) ;
  return TINY_NO_ERROR;
}

enum tiny_ErrorCode tiny_UpdateLinearCost(tiny_AdmmWorkspace* work) {
  int N = work->data->model[0].nhorizon;
  for (int k = 0; k < N - 1; ++k) {
    /* Compute q[k] = -Q*Xref[k] */  
    (work->data->q[k]).noalias() = -(*(work->data->Q)).lazyProduct(work->data->Xref[k]);

    /* Compute r[k] = -R*Uref[k] */ 
    (work->data->r[k]).noalias() = -(*(work->data->R)).lazyProduct(work->data->Uref[k]);
  }
  /* Compute q[N-1] = -Pinf*Xref[N-1] */ 
  (work->data->q[N-1]).noalias() = -(*(work->soln->Pinf)).lazyProduct(work->data->Xref[N-1]);
  // (work->soln->p[N-1]) = work->data->q[N-1];
  // PrintMatrixT(work->data->q[N-1]);
  return TINY_NO_ERROR;
}

enum tiny_ErrorCode tiny_UpdateConstrainedLinearCost(tiny_AdmmWorkspace* work) {
  int N = work->data->model[0].nhorizon;

  // if (work->stgs->en_cstr_inputs) {
  for (int k = 0; k < N - 1; ++k) {
    /* Compute r_tilde[k] = r[k] - ρ*(z[k]-y[k]) */ 
    work->data->r_tilde[k] = work->data->r[k] - work->rho * (work->ZU_new[k] - work->soln->YU[k]);
  }
  // }

  // if (work->stgs->en_cstr_states) {
  for (int k = 0; k < N - 1; ++k) {
    /* Compute q_tilde[k] = q[k] - ρ*(z[k]-y[k]) */ 
    work->data->q_tilde[k] = work->data->q[k] - work->rho * (work->ZX_new[k] - work->soln->YX[k]);   
    }  
  work->soln->p[N-1] = work->data->q[N-1] - work->rho * (work->ZX_new[N-1] - work->soln->YX[N-1]);  
    // work->data->q_tilde[N-1] = work->soln->p[N-1];
  // }
  return TINY_NO_ERROR;
}

# ifdef __cplusplus
}
# endif // ifdef __cplusplus