#include "auxil.h"

#ifdef __cplusplus
extern "C" {
#endif

enum tiny_ErrorCode tiny_InitSettings(tiny_AdmmSettings* stgs) {
  stgs->reg_min       = (float)REG_MIN;
  stgs->reg_max       = (float)REG_MAX;
  stgs->reg_mul       = (float)REG_MUL;
  stgs->en_reg_update = EN_REG_UPDATE;

  stgs->rho_init  = (float)RHO_INIT;
  stgs->rho_max   = (float)RHO_MAX;
  stgs->rho_mul   = (float)RHO_MUL;

  stgs->alpha_mul = (float)ALPHA_MUL;

  stgs->max_iter          = MAX_ITER;
  stgs->max_iter_riccati  = MAX_ITER_RICCATI;
  stgs->max_iter_ls       = MAX_ITER_LS;

  stgs->tol_abs_prim    = (float)TOL_ABS_PRIM;
  stgs->tol_abs_dual    = (float)TOL_ABS_DUAL;

  stgs->en_cstr_states  = EN_CSTR_STATES;
  stgs->en_cstr_inputs  = EN_CSTR_INPUTS;
  stgs->en_cstr_goal    = EN_CSTR_GOAL;

  stgs->verbose           = VERBOSE;
  stgs->adaptive_horizon  = ADAPTIVE_HORIZON;
  stgs->check_riccati     = CHECK_RICCATI;
  stgs->check_termination = CHECK_TERMINATION;
  stgs->warm_start        = WARM_START;
  stgs->time_limit        = TIME_LIMIT;

  return TINY_NO_ERROR;
}

enum tiny_ErrorCode tiny_SetUnconstrained(tiny_AdmmSettings* stgs) {
  stgs->en_cstr_states  = 0;
  stgs->en_cstr_inputs  = 0;
  stgs->en_cstr_goal    = 0;
  stgs->check_termination = 0;

  return TINY_NO_ERROR;
}

enum tiny_ErrorCode tiny_InitSolnTrajFromArray(tiny_AdmmWorkspace* work,
Eigen::Vector12f* X, Eigen::Vector4f* U,
float* X_data, float* U_data) {

  int N = work->data->model->nhorizon;
  int n = work->data->model[0].nstates;
  int m = work->data->model->ninputs;

  work->soln->X = X;  
  work->soln->U = U;

  for (int i = 0; i < N - 1; ++i) {
    U[i] = Eigen::Map<Eigen::Vector4f>(&U_data[i * m]);
  }

  for (int i = 0; i < N; ++i) {
    X[i] = Eigen::Map<Eigen::Vector12f>(&X_data[i * n]);
  }

  return TINY_NO_ERROR;  
}

enum tiny_ErrorCode tiny_InitSolnDualsFromArray(tiny_AdmmWorkspace* work,
Eigen::Vector12f* YX, Eigen::Vector4f* YU,
float* YX_data, float* YU_data, float* YG_data) {

  int N = work->data->model->nhorizon;
  int n = work->data->model[0].nstates;
  int m = work->data->model->ninputs;

  work->soln->YX = YX;  
  work->soln->YU = YU;
  if (YG_data) {
    work->soln->YG = Eigen::Map<Eigen::Vector12f>(YG_data);
  }

  if (YU_data) {
    for (int i = 0; i < N - 1; ++i) {
      YU[i] = Eigen::Map<Eigen::Vector4f>(&YU_data[i * m]);
    }
  }

  if (YX_data) {
    for (int i = 0; i < N; ++i) {
      YX[i] = Eigen::Map<Eigen::Vector12f>(&YX_data[i * n]);
    }
  }

  return TINY_NO_ERROR;  
}

enum tiny_ErrorCode tiny_InitSolnGainsFromArray(tiny_AdmmWorkspace* work, 
Eigen::Vector4f* d, Eigen::Vector12f* p, float* d_data, float* p_data, 
float* Kinf_data, float* Pinf_data) {

  int N = work->data->model->nhorizon;
  int n = work->data->model[0].nstates;
  int m = work->data->model->ninputs;

  work->soln->Kinf = Eigen::Map<Eigen::Matrix4x12f>(Kinf_data);
  work->soln->d = d;
  work->soln->Pinf = Eigen::Map<Eigen::Matrix12f>(Pinf_data);
  work->soln->p = p;

  for (int i = 0; i < N - 1; ++i) {
    d[i] = Eigen::Map<Eigen::Vector4f>(&d_data[i * m]);
  }
  for (int i = 0; i < N; ++i) {
    p[i] = Eigen::Map<Eigen::Vector12f>(&p_data[i * n]);
  }

  return TINY_NO_ERROR;  
}

enum tiny_ErrorCode tiny_InitDataQuadCostFromArray(tiny_AdmmWorkspace* work, 
float* Q_data, 
float* R_data) {

  int n = work->data->model[0].nstates;
  int m = work->data->model->ninputs;
  work->data->Q = Eigen::Map<Eigen::Matrix12f>(Q_data);
  work->data->R = Eigen::Map<Eigen::Matrix4f>(R_data);
  return TINY_NO_ERROR;
}

enum tiny_ErrorCode tiny_InitDataLinearCostFromArray(tiny_AdmmWorkspace* work, 
Eigen::Vector12f* q, Eigen::Vector4f* r, Eigen::Vector4f* r_tilde, float* q_data, float* r_data, float* r_tilde_data) {

  int N = work->data->model[0].nhorizon;
  int n = work->data->model[0].nstates;
  int m = work->data->model[0].ninputs;
  work->data->q = q;
  work->data->r = r;
  work->data->r_tilde = r_tilde;
  
  for (int i = 0; i < N - 1; ++i) {
    q[i] = Eigen::Map<Eigen::Vector12f>(&q_data[i * n]);
    r[i] = Eigen::Map<Eigen::Vector4f>(&r_data[i * m]);
    r_tilde[i] = Eigen::Map<Eigen::Vector4f>(&r_tilde_data[i * m]);
  }

  return TINY_NO_ERROR;
}

enum tiny_ErrorCode tiny_InitWorkspace(tiny_AdmmWorkspace* work,
                                       tiny_AdmmInfo* info,
                                       tiny_Model* model,
                                       tiny_AdmmData* data,
                                       tiny_AdmmSolution* soln,
                                       tiny_AdmmSettings* stgs) {
  work->data = data;
  work->info = info;
  work->soln = soln;
  work->stgs = stgs;
  work->data->model = model;

  // tiny_InitSolution(work);
  // tiny_InitData(work);

  // int n = model->nstates;
  int m = model->ninputs;
  int N = model->nhorizon;

  work->reg = (float)REG_MIN;
  work->alpha = (float)ALPHA;
  work->rho = work->stgs->rho_init;

  work->data_size = m + 2*m*(N - 1);  // only input constraints
  work->first_run = 1;
  return TINY_NO_ERROR;
}

enum tiny_ErrorCode tiny_InitWorkspaceTempData(tiny_AdmmWorkspace* work, 
Eigen::Vector4f* ZU, Eigen::Vector4f* ZU_new, Eigen::Vector12f* ZX, Eigen::Vector12f* ZX_new, float* temp_data) {
  int N = work->data->model[0].nhorizon;
  int n = work->data->model[0].nstates;
  int m = work->data->model[0].ninputs;

  float* ptr = temp_data;
  // work->Quu_inv = slap_MatrixFromArray(m, m, ptr); 
  // ptr += m*m;
  work->Qu = Eigen::Map<Eigen::Vector4f>(ptr);   
  ptr += m;
  // work->AmBKt = slap_MatrixFromArray(n, n, ptr);
  // ptr += n*n; 
  // work->coeff_d2p = slap_MatrixFromArray(n, m, ptr); 
  // ptr += n*m; 

  work->ZU     = ZU;
  work->ZU_new = ZU_new;
  work->ZX     = ZX;
  work->ZX_new = ZX_new;

  if (ZU) {
    for (int i = 0; i < N - 1; ++i) {
      ZU[i] = Eigen::Map<Eigen::Vector4f>(ptr);   
      ptr += m;
      ZU_new[i] = Eigen::Map<Eigen::Vector4f>(ptr);   
      ptr += m;
    }
  }

  if (ZX) {
    for (int i = 0; i < N; ++i) {
      ZX[i] = Eigen::Map<Eigen::Vector12f>(ptr);   
      ptr += n;
      ZX_new[i] = Eigen::Map<Eigen::Vector12f>(ptr);   
      ptr += n;
    }
  }

  return TINY_NO_ERROR;
}

enum tiny_ErrorCode tiny_InitPrimalCache(tiny_AdmmWorkspace* work, 
float* Quu_inv_data, float* AmBKt_data, float* coeff_d2p_data) {

  int n = work->data->model[0].nstates;
  int m = work->data->model->ninputs;

  // work->Quu_inv   = slap_MatrixFromArray(m, m, Quu_inv_data); 
  // work->AmBKt     = slap_MatrixFromArray(n, n, AmBKt_data); 
  // work->coeff_d2p = slap_MatrixFromArray(n, m, coeff_d2p_data); 

  work->Quu_inv   = Eigen::Map<Eigen::Matrix4f>(Quu_inv_data);  
  work->AmBKt     = Eigen::Map<Eigen::Matrix12f>(AmBKt_data); 
  work->coeff_d2p = Eigen::Map<Eigen::Matrix12x4f>(coeff_d2p_data); 

  return TINY_NO_ERROR;
}

enum tiny_ErrorCode tiny_ResetInfo(tiny_AdmmWorkspace* work) {
  work->info->iter = 0;
  work->info->iter_riccati = 0;
  work->info->status_val = TINY_UNSOLVED;

  work->info->obj_val = 0.0;
  work->info->pri_res = 0.0;
  work->info->dua_res = 0.0;
  return TINY_NO_ERROR;
}

// enum tiny_ErrorCode tiny_SetStateReference(tiny_AdmmWorkspace* work, Matrix* Xref, 
// float* Xref_data) {
//   int n = work->data->model[0].nstates;
//   int N = work->data->model[0].nhorizon;
//   work->data->Xref = Xref;
//   for (int i = 0; i < N; ++i) {
//     Xref[i] = slap_MatrixFromArray(n, 1, &Xref_data[i * n]);
//   }
//   return TINY_NO_ERROR;
// }

// enum tiny_ErrorCode tiny_SetInputReference(tiny_AdmmWorkspace* work, Matrix* Uref,
// float* Uref_data) {
//   int m = work->data->model[0].ninputs;
//   int N = work->data->model[0].nhorizon;
//   work->data->Uref = Uref;
//   for (int i = 0; i < N - 1; ++i) {
//     Uref[i] = slap_MatrixFromArray(m, 1, &Uref_data[i * m]);
//   }
//   return TINY_NO_ERROR;
// }

// enum tiny_ErrorCode tiny_SetReference(tiny_AdmmWorkspace* work, Matrix* Xref, 
// Matrix* Uref, float* Xref_data, float* Uref_data) {

//   tiny_SetStateReference(work, Xref, Xref_data);
//   tiny_SetInputReference(work, Uref, Uref_data);
//   return TINY_NO_ERROR;
// }

enum tiny_ErrorCode tiny_SetGoalReference(tiny_AdmmWorkspace* work, Eigen::Vector12f* Xref,
Eigen::Vector4f* Uref, float* xg_data, float* ug_data) {
  int n = work->data->model[0].nstates;
  int m = work->data->model[0].ninputs;
  int N = work->data->model[0].nhorizon;
  work->data->Xref = Xref;
  work->data->Uref = Uref;
  for (int i = 0; i < N; ++i) {
    Xref[i] =  Eigen::Map<Eigen::Vector12f>(xg_data);
  }
  for (int i = 0; i < N - 1; ++i) {
    Uref[i] =  Eigen::Map<Eigen::Vector4f>(ug_data);
  }
  return TINY_NO_ERROR;
}

enum tiny_ErrorCode tiny_SetInitialState(tiny_AdmmWorkspace* work, float* x0_data) {
  int n = work->data->model[0].nstates;
  // work->data->x0 = slap_MatrixFromArray(n, 1, x0_data);
  work->data->x0 = Eigen::Map<Eigen::Vector12f>(x0_data);
  return TINY_NO_ERROR;
}

#ifdef __cplusplus
}
#endif