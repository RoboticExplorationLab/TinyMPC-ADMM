// MPC
// Scenerio: drive the Crazyflie quadrotor to track a reference
//

#include <stdlib.h>

#include "data/quadrotor_fig8.h"
#include "quadrotor.h"
#include "simpletest.h"
#include "slap/slap.h"
#include "time.h"
#include "tinympc/tinympc.h"

// Macro variables
#define H 0.02       // dt
#define NSTATES 12   // no. of states (error state)
#define NINPUTS 4    // no. of controls
#define NHORIZON 3  // horizon steps (NHORIZON states and NHORIZON-1 controls)
#define NSIM 500     // simulation steps (fixed with reference data)

int main() {
  // ===== Start MPC setup =====
  // ===== Created data =====
  sfloat x0_data[NSTATES] = {0, 1, 0, 0.1, 0, 0,
                             0, 0, 0, 0,   0, 0};  // initial state
  sfloat ug_data[NINPUTS] = {0., 0., 0., 0.};      // goal input if needed
  sfloat Xhrz_data[NSTATES * NHORIZON] = {0};      // save X for one horizon
  sfloat X_data[NSTATES * NSIM] = {0};             // save X for the whole run
  sfloat Uhrz_data[NINPUTS * (NHORIZON - 1)] = {0};
  sfloat A_data[NSTATES * NSTATES] = {
      1.000000f, 0.000000f,  0.000000f, 0.000000f, 0.000000f, 0.000000f,
      0.000000f, 0.000000f,  0.000000f, 0.000000f, 0.000000f, 0.000000f,
      0.000000f, 1.000000f,  0.000000f, 0.000000f, 0.000000f, 0.000000f,
      0.000000f, 0.000000f,  0.000000f, 0.000000f, 0.000000f, 0.000000f,
      0.000000f, 0.000000f,  1.000000f, 0.000000f, 0.000000f, 0.000000f,
      0.000000f, 0.000000f,  0.000000f, 0.000000f, 0.000000f, 0.000000f,
      0.000000f, -0.003924f, 0.000000f, 1.000000f, 0.000000f, 0.000000f,
      0.000000f, -0.392400f, 0.000000f, 0.000000f, 0.000000f, 0.000000f,
      0.003924f, 0.000000f,  0.000000f, 0.000000f, 1.000000f, 0.000000f,
      0.392400f, 0.000000f,  0.000000f, 0.000000f, 0.000000f, 0.000000f,
      0.000000f, 0.000000f,  0.000000f, 0.000000f, 0.000000f, 1.000000f,
      0.000000f, 0.000000f,  0.000000f, 0.000000f, 0.000000f, 0.000000f,
      0.020000f, 0.000000f,  0.000000f, 0.000000f, 0.000000f, 0.000000f,
      1.000000f, 0.000000f,  0.000000f, 0.000000f, 0.000000f, 0.000000f,
      0.000000f, 0.020000f,  0.000000f, 0.000000f, 0.000000f, 0.000000f,
      0.000000f, 1.000000f,  0.000000f, 0.000000f, 0.000000f, 0.000000f,
      0.000000f, 0.000000f,  0.020000f, 0.000000f, 0.000000f, 0.000000f,
      0.000000f, 0.000000f,  1.000000f, 0.000000f, 0.000000f, 0.000000f,
      0.000000f, -0.000013f, 0.000000f, 0.010000f, 0.000000f, 0.000000f,
      0.000000f, -0.001962f, 0.000000f, 1.000000f, 0.000000f, 0.000000f,
      0.000013f, 0.000000f,  0.000000f, 0.000000f, 0.010000f, 0.000000f,
      0.001962f, 0.000000f,  0.000000f, 0.000000f, 1.000000f, 0.000000f,
      0.000000f, 0.000000f,  0.000000f, 0.000000f, 0.000000f, 0.010000f,
      0.000000f, 0.000000f,  0.000000f, 0.000000f, 0.000000f, 1.000000f,
  };
  sfloat B_data[NSTATES * NINPUTS] = {
      -0.000019f, -0.000001f, 0.000981f, 0.001264f,  -0.029414f, 0.004771f,
      -0.003847f, -0.000165f, 0.098100f, 0.252748f,  -5.882783f, 0.954290f,
      -0.000001f, -0.000019f, 0.000981f, 0.029044f,  -0.001057f, -0.003644f,
      -0.000138f, -0.003799f, 0.098100f, 5.808852f,  -0.211410f, -0.728857f,
      0.000019f,  0.000001f,  0.000981f, -0.001493f, 0.028771f,  0.001265f,
      0.003763f,  0.000195f,  0.098100f, -0.298680f, 5.754175f,  0.252942f,
      0.000001f,  0.000019f,  0.000981f, -0.028815f, 0.001700f,  -0.002392f,
      0.000222f,  0.003769f,  0.098100f, -5.762921f, 0.340018f,  -0.478376f,
  };
  sfloat f_data[NSTATES] = {0};
  sfloat Kinf_data[NINPUTS*NSTATES] = {
    -0.204719f,-0.035010f,0.298543f,-0.058814f,
    -0.006771f,-0.233680f,-0.030598f,0.271050f,
    1.394464f,1.394464f,1.394464f,1.394464f,
    0.033636f,1.065312f,0.142948f,-1.241896f,
    -0.928855f,-0.166981f,1.372039f,-0.276202f,
    0.266463f,-0.269349f,0.275893f,-0.273007f,
    -0.153553f,-0.026797f,0.225082f,-0.044732f,
    -0.005271f,-0.175610f,-0.023223f,0.204104f,
    0.580173f,0.580173f,0.580173f,0.580173f,
    0.002878f,0.083979f,0.011550f,-0.098408f,
    -0.072860f,-0.013767f,0.109059f,-0.022432f,
    0.280709f,-0.283762f,0.290688f,-0.287635f,
  };
  sfloat Pinf_data[NSTATES*NSTATES] = {
    376.164204f,-0.020264f,-0.000000f,0.089198f,338.786043f,0.080492f,112.737963f,-0.015070f,-0.000000f,0.006606f,1.669946f,0.083527f,
    -0.020264f,376.151433f,0.000000f,-338.729395f,-0.089198f,-0.032255f,-0.015070f,112.728425f,0.000000f,-1.665752f,-0.006606f,-0.033471f,
    -0.000000f,0.000000f,208.027074f,-0.000000f,-0.000000f,0.000000f,-0.000000f,0.000000f,16.194993f,0.000000f,-0.000000f,0.000000f,
    0.089198f,-338.729395f,-0.000000f,1449.708366f,0.402843f,0.151939f,0.066983f,-248.002042f,-0.000000f,7.286598f,0.030505f,0.157838f,
    338.786043f,-0.089198f,0.000000f,0.402843f,1449.968369f,0.379094f,248.044903f,-0.066983f,0.000000f,0.030505f,7.306424f,0.393812f,
    0.080492f,-0.032255f,0.000000f,0.151939f,0.379094f,97.767107f,0.061311f,-0.024570f,0.000000f,0.012225f,0.030496f,0.936579f,
    112.737963f,-0.015070f,-0.000000f,0.066983f,248.044903f,0.061311f,66.921344f,-0.011252f,-0.000000f,0.004991f,1.234055f,0.063631f,
    -0.015070f,112.728425f,0.000000f,-248.002042f,-0.066983f,-0.024570f,-0.011252f,66.914183f,0.000000f,-1.230856f,-0.004991f,-0.025500f,
    -0.000000f,0.000000f,16.194993f,-0.000000f,-0.000000f,0.000000f,-0.000000f,0.000000f,7.076044f,-0.000000f,0.000000f,0.000000f,
    0.006606f,-1.665752f,0.000000f,7.286598f,0.030505f,0.012225f,0.004991f,-1.230856f,0.000000f,1.051781f,0.002432f,0.012870f,
    1.669946f,-0.006606f,0.000000f,0.030505f,7.306424f,0.030496f,1.234055f,-0.004991f,0.000000f,0.002432f,1.053397f,0.032106f,
    0.083527f,-0.033471f,0.000000f,0.157838f,0.393812f,0.936579f,0.063631f,-0.025500f,0.000000f,0.012870f,0.032106f,1.481670f,
  };
  sfloat Quu_inv_data[NINPUTS*NINPUTS] = {
    0.258936f,0.106184f,0.251217f,0.105341f,
    0.106184f,0.261898f,0.102398f,0.251199f,
    0.251217f,0.102398f,0.268936f,0.099128f,
    0.105341f,0.251199f,0.099128f,0.266011f,
  };
  sfloat AmBKt_data[NSTATES*NSTATES] = {
    0.999990f,-0.000000f,0.000000f,0.000000f,0.003881f,0.000000f,0.019993f,-0.000000f,-0.000000f,0.000000f,0.000010f,0.000000f,
    -0.000000f,0.999990f,0.000000f,-0.003880f,-0.000000f,-0.000000f,-0.000000f,0.019993f,0.000000f,-0.000010f,-0.000000f,-0.000000f,
    -0.000000f,-0.000000f,0.994528f,0.000000f,0.000000f,0.000000f,-0.000000f,-0.000000f,0.017723f,0.000000f,-0.000000f,0.000000f,
    0.000027f,0.014560f,-0.000000f,0.933445f,0.000114f,0.000032f,0.000020f,0.010954f,-0.000000f,0.004739f,0.000008f,0.000033f,
    -0.014548f,-0.000027f,0.000000f,0.000114f,0.933497f,0.000079f,-0.010945f,-0.000020f,0.000000f,0.000008f,0.004743f,0.000082f,
    0.000331f,-0.000132f,0.000000f,0.000571f,0.001428f,0.996745f,0.000243f,-0.000097f,0.000000f,0.000042f,0.000106f,0.006571f,
    -0.001903f,-0.000003f,0.000000f,0.000015f,0.383701f,0.000010f,0.998568f,-0.000003f,0.000000f,0.000001f,0.001274f,0.000011f,
    -0.000003f,-0.001904f,0.000000f,-0.383695f,-0.000015f,-0.000004f,-0.000003f,0.998567f,0.000000f,-0.001274f,-0.000001f,-0.000004f,
    -0.000000f,-0.000000f,-0.547188f,0.000000f,0.000000f,0.000000f,-0.000000f,-0.000000f,0.772340f,0.000000f,-0.000000f,0.000000f,
    0.005335f,2.912025f,0.000000f,-13.310992f,0.022809f,0.006349f,0.003906f,2.190721f,-0.000000f,-0.052217f,0.001686f,0.006588f,
    -2.909586f,-0.005335f,0.000000f,0.022808f,-13.300594f,0.015893f,-2.188937f,-0.003906f,0.000000f,0.001686f,-0.051449f,0.016491f,
    0.066195f,-0.026455f,0.000000f,0.114111f,0.285517f,-0.650985f,0.048672f,-0.019452f,0.000000f,0.008465f,0.021180f,0.314175f,
  };
  sfloat coeff_d2p_data[NSTATES*NINPUTS] = {
    -0.000001f,0.000000f,-0.000000f,-0.000001f,-0.000003f,-0.000879f,-0.000001f,0.000000f,-0.000000f,-0.000000f,-0.000000f,-0.000009f,
    0.000001f,-0.000000f,-0.000000f,0.000001f,0.000003f,0.000728f,0.000000f,-0.000000f,-0.000000f,0.000000f,0.000000f,0.000007f,
    -0.000000f,0.000000f,-0.000000f,-0.000001f,-0.000002f,-0.000409f,-0.000000f,0.000000f,-0.000000f,-0.000000f,-0.000000f,-0.000004f,
    0.000000f,0.000000f,-0.000000f,0.000002f,0.000002f,0.000560f,0.000000f,-0.000000f,-0.000000f,0.000000f,0.000000f,0.000005f,
  };
  sfloat d_data[NINPUTS * (NHORIZON - 1)] = {0};
  sfloat p_data[NSTATES * NHORIZON] = {0};
  sfloat Q_data[NSTATES * NSTATES] = {0};
  sfloat R_data[NINPUTS * NINPUTS] = {0};
  sfloat q_data[NSTATES*(NHORIZON-1)] = {0};
  sfloat r_data[NINPUTS*(NHORIZON-1)] = {0};
  sfloat r_tilde_data[NINPUTS*(NHORIZON-1)] = {0};

  sfloat umin_data[NINPUTS] = {-0.5, -0.5, -0.5, -0.5};
  sfloat umax_data[NINPUTS] = {0.5, 0.5, 0.5, 0.5};
  // Put constraints on u, x
  sfloat Acu_data[NINPUTS * NINPUTS] = {0};  
  sfloat YU_data[NINPUTS * (NHORIZON - 1)] = {0};

  // ===== Created matrices =====
  Matrix X[NSIM];
  Matrix Xref[NSIM];
  Matrix Uref[NSIM - 1];
  Matrix Xhrz[NHORIZON];
  Matrix Uhrz[NHORIZON - 1];
  Matrix d[NHORIZON - 1];
  Matrix p[NHORIZON];
  Matrix YU[NHORIZON - 1];
  Matrix ZU[NHORIZON - 1];
  Matrix ZU_new[NHORIZON - 1];
  Matrix q[NHORIZON-1];
  Matrix r[NHORIZON-1];
  Matrix r_tilde[NHORIZON-1];
  Matrix A;
  Matrix B;
  Matrix f;

  for (int i = 0; i < NSIM; ++i) {
    if (i < NSIM - 1) {
      Uref[i] = slap_MatrixFromArray(NINPUTS, 1, ug_data);
      // PrintMatrix(Uref[i]);
    }
    X[i] = slap_MatrixFromArray(NSTATES, 1, &X_data[i * NSTATES]);
    Xref[i] = slap_MatrixFromArray(NSTATES, 1, &Xref_data[i * NSTATES]);
    // PrintMatrix(Xref[i]);
  }

  // ===== Created tinyMPC struct =====
  tiny_Model model;
  tiny_InitModel(&model, NSTATES, NINPUTS, NHORIZON, 0, 0, 0.1);
  // tiny_InitModel(&model, NSTATES, NINPUTS, NHORIZON, 0, 1, 0.1);
  tiny_AdmmSettings stgs;
  tiny_InitSettings(&stgs);  //if switch on/off during run, initialize all
  stgs.rho_init = 1e0;

  tiny_AdmmData data;
  tiny_AdmmInfo info;
  tiny_AdmmSolution soln;
  tiny_AdmmWorkspace work;
  tiny_InitWorkspace(&work, &info, &model, &data, &soln, &stgs);
  
  // ===== Fill in the remaining struct =====
  sfloat temp_data[work.data_size];
  T_INIT_ZEROS(temp_data);
  tiny_InitWorkspaceTempData(&work, ZU, ZU_new, 0, 0, temp_data);
  tiny_InitPrimalCache(&work, Quu_inv_data, AmBKt_data, coeff_d2p_data);

  tiny_InitModelFromArray(&model, &A, &B, &f, A_data, B_data, f_data);
  tiny_InitSolnTrajFromArray(&work, Xhrz, Uhrz, Xhrz_data, Uhrz_data);
  tiny_InitSolnGainsFromArray(&work, d, p, d_data, p_data, Kinf_data, Pinf_data);
  tiny_InitSolnDualsFromArray(&work, 0, YU, 0, YU_data, 0);

  tiny_SetInitialState(&work, x0_data);  
  data.Xref = Xref;
  data.Uref = Uref;
  tiny_InitDataQuadCostFromArray(&work, Q_data, R_data);
  // slap_SetIdentity(prob.Q, 1000e-1);
  sfloat Qdiag[NSTATES] = {10, 10, 10, 1, 1, 1, 1, 1, 1, 0.1, 0.1, 0.1};
  slap_SetDiagonal(data.Q, Qdiag, NSTATES);
  slap_SetIdentity(data.R, 1);
  slap_AddIdentity(data.R, work.rho); // \tilde{R}
  tiny_InitDataLinearCostFromArray(&work, q, r, r_tilde, q_data, r_data, r_tilde_data);

  // Set up constraints
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
  }

  stgs.en_cstr_goal = 0;
  stgs.en_cstr_inputs = 1;
  stgs.en_cstr_states = 0;
  stgs.max_iter = 8;
  stgs.verbose = 0;
  stgs.check_termination = 2;
  stgs.tol_abs_dual = 1e-2;
  stgs.tol_abs_prim = 1e-2;

  // Absolute formulation:
  // Warm-starting since horizon data is reused
  // At each time step (stop earlier as horizon exceeds the end)
  slap_Copy(X[0], work.data->x0);  
  srand(1);  // random seed
  // ----- End MPC setup -----

  // ===== Start MPC loop =====
  for (int k = 0; k < NSIM - NHORIZON - 1; ++k) {
    // printf("\n=> k = %d\n", k);
    Matrix pose = slap_CreateSubMatrix(X[k], 0, 0, 6, 1);
    Matrix pose_ref = slap_CreateSubMatrix(Xref[k], 0, 0, 6, 1);
    // printf("ex[%d] = %.4f\n", k, slap_NormedDifference(X[k], Xref[k]));
    printf("ex[%d] = %.4f\n", k, slap_NormedDifference(pose, pose_ref));
    // printf("%.4f\n", slap_NormedDifference(pose, pose_ref));

    // === 1. Setup and solve MPC ===
    for (int j = 0; j < NSTATES; ++j) {
      X[k].data[j] += X[k].data[j] * T_NOISE(5);
    }
    slap_Copy(work.data->x0, X[k]);  // update current measurement

    // Update reference
    data.Xref = &Xref[k];
    data.Uref = &Uref[k];
    tiny_UpdateLinearCost(&work);

    clock_t start, end;
    double cpu_time_used;
    start = clock();

    // Solve optimization problem using Augmented Lagrangian TVLQR
    tiny_SolveAdmm(&work);

    end = clock();
    cpu_time_used = ((double)(end - start)) * 1000 / CLOCKS_PER_SEC;  // ms
    printf("solve time:       %f\n", cpu_time_used);
    // printf("%f\n", cpu_time_used);

    // if(work.info->status_val != TINY_SOLVED) {
    //   printf("!!! STOP AS SOLVER FAILED !!!\n");
    //   return 0;
    // }

    // Test control constraints here (since we didn't save U)
    for (int i = 0; i < NINPUTS; ++i) {
      TEST(Uhrz[0].data[i] < umax_data[i] + stgs.tol_abs_dual);
      TEST(Uhrz[0].data[i] > umin_data[i] - stgs.tol_abs_dual);
    }
    // PrintMatrixT(Uhrz[0]);

    // Matrix pos = slap_CreateSubMatrix(X[k], 0, 0, 3, 1);
    // PrintMatrixT(pos);

    // === 2. Simulate dynamics using the first control solution ===
    // tiny_QuadNonlinearDynamics(&X[k + 1], X[k], Uref[k]);
    tiny_QuadNonlinearDynamics(&X[k + 1], X[k], Uhrz[0]);
    // tiny_DynamicsLti(&X[k + 1], X[k], Uref[k], model);
    tiny_ShiftFill(Uhrz, T_ARRAY_SIZE(Uhrz));
  }

  // ========== Test ==========
  // Test state constraints
  // for (int k = 0; k < NSIM - NHORIZON - 1; ++k) {
  //   for (int i = 0; i < NSTATES; ++i) {
  //     TEST(X[k].data[i] < bcx_data[i] + stgs.tol_abs_dual);
  //     TEST(X[k].data[i] > -bcx_data[i] - stgs.tol_abs_dual);
  //   }
  // }
  // Test tracking performance
  for (int k = NSIM - NHORIZON - 5; k < NSIM - NHORIZON; ++k) {
    TEST(slap_NormedDifference(X[k], Xref[k]) < 0.5);
  }
  // --------------------------

  PrintTestResult();
  return TestResult();
}
