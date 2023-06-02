// Task: Test AL-LQR on double integrator with input/state box constraints and
// goal constraint. Scenerio: drive from initial state to goal state.

#include <time.h>

#include "simpletest.h"
#include "slap/slap.h"
#include "test_utils.h"
#include "tinympc/admm.h"
#include "tinympc/auxil.h"
#include "tinympc/constraint_linear.h"

#define NSTATES 4
#define NINPUTS 2
#define NHORIZON 51

void MpcLtiTest() {
  float A_data[NSTATES * NSTATES] = {1,   0, 0, 0, 0, 1,   0, 0,
                                      0.1, 0, 1, 0, 0, 0.1, 0, 1};
  float B_data[NSTATES * NINPUTS] = {0.005, 0, 0.1, 0, 0, 0.005, 0, 0.1};
  float f_data[NSTATES] = {0};
  // float x0_data[NSTATES] = {5, 7, 2, -1.4};
  float x0_data[NSTATES] = {1, 0, 0, 1.0};  
  // float xg_data[NSTATES] = {2, 5, 0, 0};
  float xg_data[NSTATES] = {0};
  float ug_data[NINPUTS] = {0};
  float X_data[NSTATES * NHORIZON] = {0};
  float U_data[NINPUTS * (NHORIZON - 1)] = {0};
  float Kinf_data[NINPUTS*NSTATES] = {
    1.255647f,0.000000f,
    0.000000f,1.255647f,
    2.021867f,0.000000f,
    0.000000f,2.021867f,
  };
  float Pinf_data[NSTATES*NSTATES] = {
    161.021861f,0.000000f,71.589105f,0.000000f,
    0.000000f,161.021861f,0.000000f,71.589105f,
    71.589105f,0.000000f,116.694654f,0.000000f,
    0.000000f,71.589105f,0.000000f,116.694654f,
  };
  float Quu_inv_data[NINPUTS*NINPUTS] = {
    0.157665f,0.000000f,
    0.000000f,0.157665f,
  };
  float AmBKt_data[NSTATES*NSTATES] = {
    0.993722f,0.000000f,0.089891f,0.000000f,
    0.000000f,0.993722f,0.000000f,0.089891f,
    -0.125565f,0.000000f,0.797813f,0.000000f,
    0.000000f,-0.125565f,0.000000f,0.797813f,
  };
  float coeff_d2p_data[NSTATES*NINPUTS] = {
    0.000000f,0.000000f,-0.000000f,0.000000f,
    0.000000f,0.000000f,0.000000f,-0.000000f,
  };
  float d_data[NINPUTS * (NHORIZON - 1)] = {0};
  float p_data[NSTATES * NHORIZON] = {0};
  float Q_data[NSTATES * NSTATES] = {0};
  float R_data[NINPUTS * NINPUTS] = {0};
  float q_data[NSTATES*(NHORIZON-1)] = {0};
  float r_data[NINPUTS*(NHORIZON-1)] = {0};
  float r_tilde_data[NINPUTS*(NHORIZON-1)] = {0};

  float umin_data[NINPUTS] = {-1, -1};
  float umax_data[NINPUTS] = {1, 1};
  float Acu_data[NINPUTS * NINPUTS] = {0};  
  float YU_data[NINPUTS * (NHORIZON - 1)] = {0};

  Matrix X[NHORIZON];
  Matrix U[NHORIZON - 1];
  Matrix Xref[NHORIZON];
  Matrix Uref[NHORIZON - 1];
  Matrix d[NHORIZON - 1];
  Matrix p[NHORIZON];
  Matrix YU[NHORIZON - 1];
  // Matrix YX[NHORIZON];
  Matrix ZU[NHORIZON - 1];
  Matrix ZU_new[NHORIZON - 1];
  Matrix q[NHORIZON-1];
  Matrix r[NHORIZON-1];
  Matrix r_tilde[NHORIZON-1];
  Matrix A;
  Matrix B;
  Matrix f;

  tiny_Model model;
  tiny_InitModel(&model, NSTATES, NINPUTS, NHORIZON, 0, 0, 0.1);
  tiny_AdmmSettings stgs;
  tiny_InitSettings(&stgs);  //if switch on/off during run, initialize all
  stgs.rho_init = 5.0;
  tiny_AdmmData data;
  tiny_AdmmInfo info;
  tiny_AdmmSolution soln;
  tiny_AdmmWorkspace work;  
  tiny_InitWorkspace(&work, &info, &model, &data, &soln, &stgs);
  
  float temp_data[work.data_size];
  T_INIT_ZEROS(temp_data);

  tiny_InitWorkspaceTempData(&work, ZU, ZU_new, 0, 0, temp_data);
  tiny_InitPrimalCache(&work, Quu_inv_data, AmBKt_data, coeff_d2p_data);
  
  tiny_InitModelFromArray(&model, &A, &B, &f, A_data, B_data, f_data);
  tiny_InitSolnTrajFromArray(&work, X, U, X_data, U_data);
  tiny_InitSolnGainsFromArray(&work, d, p, d_data, p_data, Kinf_data, Pinf_data);
  tiny_InitSolnDualsFromArray(&work, 0, YU, 0, YU_data, 0);

  tiny_SetInitialState(&work, x0_data);  
  tiny_SetGoalReference(&work, Xref, Uref, xg_data, ug_data);

  tiny_InitDataQuadCostFromArray(&work, Q_data, R_data);
  slap_SetIdentity(data.Q, 10);  
  slap_SetIdentity(data.R, 0.1);
  slap_AddIdentity(data.R, work.rho); // \tilde{R}
  tiny_InitDataLinearCostFromArray(&work, q, r, r_tilde, q_data, r_data, r_tilde_data);

  tiny_SetInputBound(&work, Acu_data, umin_data, umax_data);
  tiny_UpdateLinearCost(&work);

  if (0) {
    printf("\nProblem Info: \n");
    PrintMatrix(work.data->model->A[0]);
    PrintMatrix(work.data->model->B[0]);
    // PrintMatrix(work.data->model->f[0]);
    PrintMatrix(work.data->Q);
    PrintMatrix(work.data->R);
    PrintMatrixT(work.data->x0);
    PrintMatrixT(work.data->Xref[NHORIZON-5]);
    PrintMatrixT(work.data->Uref[NHORIZON-5]);
    PrintMatrixT(work.data->q[NHORIZON-5]);
    PrintMatrixT(work.data->r[NHORIZON-5]);
    PrintMatrixT(work.soln->Kinf);
  }

  stgs.en_cstr_goal = 0;
  stgs.en_cstr_inputs = 1;
  stgs.en_cstr_states = 0;
  stgs.max_iter = 200;
  stgs.verbose = 0;
  stgs.check_termination = 10;
  
  clock_t start, end;
  double cpu_time_used;
  start = clock();
  tiny_SolveAdmm(&work);
  end = clock();
  cpu_time_used = ((double)(end - start)) * 1000 / CLOCKS_PER_SEC;
  // printf("time: %f ms\n", cpu_time_used);

  if (0) {
    for (int k = 0; k < NHORIZON - 1; ++k) {
      printf("\n=>k = %d\n", k);
      // PrintMatrix(p[k]);
      // PrintMatrixT(Xref[k]);
      PrintMatrixT(U[k]);
      PrintMatrixT(X[k]);
    }
    PrintMatrixT(X[NHORIZON - 1]);
  }  

  // ========== Test ==========
  for (int k = 0; k < NHORIZON - 1; ++k) {
    // for (int i = 0; i < NSTATES; ++i) {
    //   TEST(X[k].data[i] < xmax_data[i] + stgs.tol_abs_dual);
    //   TEST(X[k].data[i] > xmin_data[i] - stgs.tol_abs_dual);
    // }
    for (int i = 0; i < NINPUTS; ++i) {
      TEST(U[k].data[i] > umin_data[i] - stgs.tol_abs_dual);
      TEST(U[k].data[i] < umax_data[i] + stgs.tol_abs_dual);
    }
  }
  TEST(SumOfSquaredError(X[NHORIZON - 1].data, xg_data, NSTATES) < 1e-1);
}

int main() {
  printf("=== AL LQR LTI Test ===\n");
  MpcLtiTest();
  PrintTestResult();
  return TestResult();
}
