// MPC
// Scenerio: make Crazyflie hovering
//

#include "quadrotor.h"
#include "slap/slap.h"
#include "time.h"
#include "tinympc/tinympc.h"

// Macro variables
#define DT 0.01       // dt
#define NSTATES 12   // no. of states (error state)
#define NINPUTS 4    // no. of controls
#define NHORIZON 3  // horizon steps (NHORIZON states and NHORIZON-1 controls)
#define NSIM 200     // simulation steps (fixed with reference data)

int main() {
  /* Start MPC initialization*/

  // Create data array 
  float x0_data[NSTATES] = {0, 0, 0, 1, 0, 0,
                             0, 0, 0, 0, 0, 0};  // initial state
  float xg_data[NSTATES] = {0};  
  float ug_data[NINPUTS] = {0};      // goal input if needed
  float Xhrz_data[NSTATES * NHORIZON] = {0};      // save X for one horizon
  float X_data[NSTATES * NSIM] = {0};             // save X for the whole run
  float Uhrz_data[NINPUTS * (NHORIZON - 1)] = {0};
  float A_data[NSTATES*NSTATES] = {
    1.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,
    0.000000f,1.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,
    0.000000f,0.000000f,1.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,
    0.000000f,-0.000962f,0.000000f,1.000000f,0.000000f,0.000000f,0.000000f,-0.196200f,0.000000f,0.000000f,0.000000f,0.000000f,
    0.000962f,0.000000f,0.000000f,0.000000f,1.000000f,0.000000f,0.196200f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,
    0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,1.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,
    0.010000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,1.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,
    0.000000f,0.010000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,1.000000f,0.000000f,0.000000f,0.000000f,0.000000f,
    0.000000f,0.000000f,0.010000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,1.000000f,0.000000f,0.000000f,0.000000f,
    0.000000f,-0.000002f,0.000000f,0.005000f,0.000000f,0.000000f,0.000000f,-0.000500f,0.000000f,1.000000f,0.000000f,0.000000f,
    0.000002f,0.000000f,0.000000f,0.000000f,0.005000f,0.000000f,0.000500f,0.000000f,0.000000f,0.000000f,1.000000f,0.000000f,
    0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.005000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,1.000000f,
  };

  float B_data[NSTATES*NINPUTS] = {
    -0.000001f,0.000001f,0.000241f,-0.006729f,-0.006718f,-0.000157f,-0.000456f,0.000457f,0.048118f,-2.691475f,-2.687185f,-0.062638f,
    0.000001f,0.000001f,0.000241f,-0.007461f,0.007394f,0.000464f,0.000502f,0.000507f,0.048118f,-2.984480f,2.957661f,0.185499f,
    0.000001f,-0.000001f,0.000241f,0.006785f,0.006876f,-0.001324f,0.000467f,-0.000461f,0.048118f,2.714004f,2.750267f,-0.529510f,
    -0.000001f,-0.000001f,0.000241f,0.007405f,-0.007552f,0.001017f,-0.000513f,-0.000503f,0.048118f,2.961950f,-3.020743f,0.406649f,
  };
  float f_data[NSTATES] = {0};

  float Kinf_data[NINPUTS*NSTATES] = {
    -0.141668f,0.140640f,0.084249f,-0.083221f,
    0.126596f,0.103198f,-0.104505f,-0.125288f,
    1.542556f,1.542556f,1.542556f,1.542556f,
    -0.658536f,-0.524130f,0.534755f,0.647910f,
    -0.744417f,0.735134f,0.424451f,-0.415168f,
    -0.530601f,0.594336f,-0.754518f,0.690783f,
    -0.125704f,0.124566f,0.073728f,-0.072590f,
    0.111952f,0.090537f,-0.091916f,-0.110573f,
    0.505438f,0.505438f,0.505438f,0.505438f,
    -0.050258f,-0.039023f,0.040098f,0.049183f,
    -0.057398f,0.056417f,0.031262f,-0.030281f,
    -0.225917f,0.251652f,-0.316332f,0.290597f,
  };
  float Pinf_data[NSTATES*NSTATES] = {
    8789.322258f,-12.626764f,0.000000f,63.952679f,7521.021129f,240.248353f,2577.368956f,-11.116774f,0.000000f,4.499760f,61.868897f,97.516668f,
    -12.626764f,8777.858717f,0.000000f,-7462.140483f,-63.963994f,-92.720769f,-11.117516f,2567.216299f,0.000000f,-57.724868f,-4.500876f,-37.637455f,
    -0.000000f,0.000000f,327662.866955f,-0.000000f,-0.000000f,-0.000000f,-0.000000f,0.000000f,32043.163369f,-0.000000f,-0.000000f,-0.000000f,
    63.952679f,-7462.140483f,-0.000000f,35804.776438f,338.665751f,513.456568f,57.119328f,-6422.597373f,-0.000000f,292.033792f,24.843787f,210.172462f,
    7521.021129f,-63.963994f,0.000000f,338.665751f,36124.828480f,1323.784196f,6475.711343f,-57.125516f,0.000000f,24.842231f,315.794800f,541.690966f,
    240.248353f,-92.720769f,-0.000000f,513.456568f,1323.784196f,204805.908985f,216.891848f,-83.855656f,0.000000f,40.268235f,103.455293f,2579.674924f,
    2577.368956f,-11.117516f,0.000000f,57.119328f,6475.711343f,216.891848f,1887.863355f,-9.835203f,0.000000f,4.061864f,54.362070f,88.176819f,
    -11.116774f,2567.216299f,0.000000f,-6422.597373f,-57.125516f,-83.855656f,-9.835203f,1878.815746f,0.000000f,-50.571678f,-4.062576f,-34.095133f,
    -0.000000f,-0.000000f,32043.163369f,0.000000f,-0.000000f,0.000000f,-0.000000f,-0.000000f,10539.139141f,0.000000f,-0.000000f,-0.000000f,
    4.499760f,-57.724868f,-0.000000f,292.033792f,24.842231f,40.268235f,4.061864f,-50.571678f,-0.000000f,27.238629f,2.000591f,17.193802f,
    61.868897f,-4.500876f,0.000000f,24.843787f,315.794800f,103.455293f,54.362070f,-4.062576f,0.000000f,2.000591f,29.219562f,44.108880f,
    97.516668f,-37.637455f,-0.000000f,210.172462f,541.690966f,2579.674924f,88.176819f,-34.095133f,-0.000000f,17.193802f,44.108880f,1281.585482f,
  };
  float Quu_inv_data[NINPUTS*NINPUTS] = {
    0.000718f,0.000039f,0.000149f,0.000046f,
    0.000039f,0.000690f,0.000064f,0.000159f,
    0.000149f,0.000064f,0.000654f,0.000085f,
    0.000046f,0.000159f,0.000085f,0.000662f,
  };
  float AmBKt_data[NSTATES*NSTATES] = {
    0.999999f,-0.000000f,0.000000f,0.000000f,0.000960f,0.000000f,0.010000f,-0.000000f,0.000000f,0.000000f,0.000001f,0.000000f,
    -0.000000f,0.999999f,0.000000f,-0.000960f,-0.000000f,-0.000000f,-0.000000f,0.010000f,0.000000f,-0.000001f,-0.000000f,-0.000000f,
    -0.000000f,-0.000000f,0.998516f,-0.000000f,-0.000000f,0.000000f,0.000000f,-0.000000f,0.009514f,-0.000000f,0.000000f,-0.000000f,
    0.000141f,0.003259f,-0.000000f,0.983232f,0.000670f,0.000868f,0.000121f,0.002871f,-0.000000f,0.003734f,0.000047f,0.000352f,
    -0.003199f,-0.000140f,0.000000f,0.000668f,0.983510f,0.002445f,-0.002821f,-0.000120f,-0.000000f,0.000047f,0.003754f,0.000991f,
    0.000109f,-0.000039f,0.000000f,0.000189f,0.000526f,0.997940f,0.000094f,-0.000034f,0.000000f,0.000013f,0.000037f,0.004134f,
    -0.000217f,-0.000010f,0.000000f,0.000045f,0.195081f,0.000166f,0.999809f,-0.000008f,0.000000f,0.000003f,0.000415f,0.000067f,
    -0.000010f,-0.000221f,0.000000f,-0.195062f,-0.000046f,-0.000059f,-0.000008f,0.999805f,0.000000f,-0.000414f,-0.000003f,-0.000024f,
    -0.000000f,-0.000000f,-0.296899f,-0.000000f,-0.000000f,0.000000f,0.000000f,0.000000f,0.902717f,0.000000f,0.000000f,-0.000000f,
    0.056286f,1.303447f,-0.000000f,-6.707093f,0.268158f,0.347386f,0.048345f,1.148494f,-0.000000f,0.493768f,0.018734f,0.140793f,
    -1.279752f,-0.056083f,0.000000f,0.267042f,-6.596131f,0.978138f,-1.128262f,-0.048161f,-0.000000f,0.018652f,0.501446f,0.396432f,
    0.043490f,-0.015602f,-0.000000f,0.075663f,0.210584f,-0.823915f,0.037578f,-0.013489f,-0.000000f,0.005323f,0.014807f,0.653496f,
  };
  float coeff_d2p_data[NSTATES*NINPUTS] = {
    0.032891f,-0.030505f,-0.000003f,0.043732f,0.048005f,0.219590f,0.013422f,-0.012356f,-0.000001f,0.000412f,0.000532f,0.003710f,
    -0.033762f,-0.028101f,-0.000003f,0.038800f,-0.048985f,-0.279279f,-0.013738f,-0.011207f,-0.000001f,0.000227f,-0.000523f,-0.004084f,
    -0.024565f,0.027276f,-0.000003f,-0.037878f,-0.032999f,0.439999f,-0.009701f,0.010909f,-0.000001f,-0.000235f,-0.000087f,0.005100f,
    0.025436f,0.031329f,-0.000003f,-0.044655f,0.033979f,-0.380310f,0.010017f,0.012653f,-0.000001f,-0.000403f,0.000078f,-0.004726f,
  };
  float Q_data[NSTATES*NSTATES] = {
    100.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,
    0.000000f,100.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,
    0.000000f,0.000000f,10000.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,
    0.000000f,0.000000f,0.000000f,25.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,
    0.000000f,0.000000f,0.000000f,0.000000f,25.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,
    0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,2500.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,
    0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,25.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,
    0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,25.000000f,0.000000f,0.000000f,0.000000f,0.000000f,
    0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,400.000000f,0.000000f,0.000000f,0.000000f,
    0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,11.111111f,0.000000f,0.000000f,
    0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,11.111111f,0.000000f,
    0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,400.000000f,
  };
  float R_data[NINPUTS*NINPUTS] = {
    900.000000f,0.000000f,0.000000f,0.000000f,
    0.000000f,900.000000f,0.000000f,0.000000f,
    0.000000f,0.000000f,900.000000f,0.000000f,
    0.000000f,0.000000f,0.000000f,900.000000f,
  };

  float d_data[NINPUTS * (NHORIZON - 1)] = {0};
  float p_data[NSTATES * NHORIZON] = {0};
  float q_data[NSTATES*(NHORIZON-1)] = {0};
  float r_data[NINPUTS*(NHORIZON-1)] = {0};
  float r_tilde_data[NINPUTS*(NHORIZON-1)] = {0};

  float umin_data[NINPUTS] = {0};
  float umax_data[NINPUTS] = {0};
  // Put constraints on u, x
  float Acu_data[NINPUTS * NINPUTS] = {0};  
  float YU_data[NINPUTS * (NHORIZON - 1)] = {0};

  // Created matrices
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
    X[i] = slap_MatrixFromArray(NSTATES, 1, &X_data[i * NSTATES]);
  }

  /* Create TinyMPC struct and problem data*/
  tiny_Model model;
  tiny_InitModel(&model, NSTATES, NINPUTS, NHORIZON, 0, 0, DT);
  tiny_AdmmSettings stgs;
  tiny_InitSettings(&stgs);
  stgs.rho_init = 50.0;  // Important (select offline, associated with precomp.)

  tiny_AdmmData data;
  tiny_AdmmInfo info;
  tiny_AdmmSolution soln;
  tiny_AdmmWorkspace work;
  tiny_InitWorkspace(&work, &info, &model, &data, &soln, &stgs);
  
  // Fill in the remaining struct 
  float temp_data[work.data_size];
  T_INIT_ZEROS(temp_data);
  tiny_InitWorkspaceTempData(&work, ZU, ZU_new, 0, 0, temp_data);
  tiny_InitPrimalCache(&work, Quu_inv_data, AmBKt_data, coeff_d2p_data);

  tiny_InitModelFromArray(&model, &A, &B, &f, A_data, B_data, f_data);
  tiny_InitSolnTrajFromArray(&work, Xhrz, Uhrz, Xhrz_data, Uhrz_data);
  tiny_InitSolnGainsFromArray(&work, d, p, d_data, p_data, Kinf_data, Pinf_data);
  tiny_InitSolnDualsFromArray(&work, 0, YU, 0, YU_data, 0);

  tiny_SetInitialState(&work, x0_data);  
  tiny_SetGoalReference(&work, Xref, Uref, xg_data, ug_data);

  /* Set up LQR cost */
  tiny_InitDataQuadCostFromArray(&work, Q_data, R_data);
  tiny_InitDataLinearCostFromArray(&work, q, r, r_tilde, q_data, r_data, r_tilde_data);

  /* Set up constraints */
  tiny_SetInputBound(&work, Acu_data, umin_data, umax_data);
  slap_SetConst(data.ucu, 0.5);
  slap_SetConst(data.lcu, -0.5);

  tiny_UpdateLinearCost(&work);

  if (1) {
    printf("\nProblem Info: \n");
    PrintMatrix(work.data->model->A[0]);
    PrintMatrix(work.data->model->B[0]);
    PrintMatrix(work.data->Q);
    PrintMatrix(work.data->R);
    PrintMatrixT(work.data->x0);
    PrintMatrixT(work.data->Xref[NHORIZON-5]);
    PrintMatrixT(work.data->Uref[NHORIZON-5]);
    PrintMatrixT(work.data->q[NHORIZON-5]);
    PrintMatrixT(work.data->r[NHORIZON-5]);
  }

  /* Solver settings */
  stgs.en_cstr_goal = 0;
  stgs.en_cstr_inputs = 1;
  stgs.en_cstr_states = 0;
  stgs.max_iter = 100;           // limit this if needed
  stgs.verbose = 0;
  stgs.check_termination = 1;
  stgs.tol_abs_dual = 5e-2;
  stgs.tol_abs_prim = 5e-2;

  // Absolute formulation:
  // Warm-starting since horizon data is reused
  // Stop earlier as horizon exceeds the end
  MatCpy(X[0], work.data->x0);  
  srand(1);  // random seed

  /* End of MPC initialization*/

  /* Start MPC loop */

  for (int k = 0; k < NSIM - NHORIZON - 1; ++k) {
    Matrix pose = slap_CreateSubMatrix(X[k], 0, 0, 6, 1);
    Matrix pose_ref = slap_CreateSubMatrix(Xref[0], 0, 0, 6, 1);
    // printf("ex[%d] = %.4f\n", k, slap_NormedDifference(X[k], Xref[0]));
    printf("ex[%d] =  %.4f\n", k, slap_NormedDifference(pose, pose_ref));
    // printf("%.4f\n", slap_NormedDifference(pose, pose_ref));

    // Inject noise into measurement
    for (int j = 0; j < NSTATES; ++j) {
      X[k].data[j] += X[k].data[j] * T_NOISE(0);
    }

    clock_t start, end;
    double cpu_time_used;
    start = clock();

    MatCpy(work.data->x0, X[k]);  // update current measurement

    // Warm-start by previous solution
    tiny_ShiftFill(Uhrz, T_ARRAY_SIZE(Uhrz));

    // Solve optimization problem using Augmented Lagrangian TVLQR
    tiny_SolveAdmm(&work);

    end = clock();
    cpu_time_used = ((double)(end - start)) * 1000 / CLOCKS_PER_SEC;  // ms
    printf("solve time:        %f\n", cpu_time_used);
    // printf("%f\n", cpu_time_used);

    // if(work.info->status_val != TINY_SOLVED) {
    //   printf("!!! STOP AS SOLVER FAILED !!!\n");
    //   return 0;
    // }

    // PrintMatrixT(Uhrz[0]);

    // Matrix pos = slap_CreateSubMatrix(X[k], 0, 0, 3, 1);
    // PrintMatrixT(pos);

    // === 2. Simulate dynamics using the first control solution ===
    // tiny_QuadNonlinearDynamics(&X[k + 1], X[k], Uref[k]);
    // tiny_Clamp(ZU_new[0].data, umin_data[0], umax_data[0], NINPUTS);
    // If no constraints, use Uhrz[0]
    tiny_QuadNonlinearDynamics(&X[k + 1], X[k], ZU_new[0]);
    // tiny_DynamicsLti(&X[k + 1], X[k], Uref[k], model);
  }

  return 0;
}