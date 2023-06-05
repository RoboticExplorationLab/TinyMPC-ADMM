// MPC
// Scenerio: make Crazyflie hovering
//

#include "time.h"
#include "tinympc/tinympc.h"
#include "Eigen.h"

// Macro variables
#define DT 0.01       // dt
// #define NSTATES 12   // no. of states (error state)
// #define NINPUTS 4    // no. of controls
// These are already defined in `constants.h`
#define NHORIZON 3  // horizon steps (NHORIZON states and NHORIZON-1 controls)
#define NSIM 10     // simulation steps (fixed with reference data)

using namespace Eigen;

// Create TinyMPC struct
static tiny_Model model;
static tiny_AdmmSettings stgs;
static tiny_AdmmData data;
static tiny_AdmmInfo info;
static tiny_AdmmSolution soln;
static tiny_AdmmWorkspace work;

void init_mpc() {
  /* Start MPC initialization*/

  /* Create TinyMPC struct and problem data*/
  tiny_InitModel(&model, NSTATES, NINPUTS, NHORIZON, 0, 0, DT);
  tiny_InitSettings(&stgs);
  stgs.rho_init = 50.0;  // Important (select offline, associated with precomp.)
  tiny_InitWorkspace(&work, &info, &model, &data, &soln, &stgs);
  
  // Fill in the remaining struct 
  tiny_InitWorkspaceTempData(&work, ZU, ZU_new, 0, 0);
  tiny_InitPrimalCache(&work, Quu_inv_data, AmBKt_data, coeff_d2p_data);

  tiny_InitModel(&model, &A, &B, &f, A_data, B_data, f_data);
  tiny_InitSolnTraj(&work, Xhrz, Uhrz, Xhrz_data, Uhrz_data);
  tiny_InitSolnGains(&work, d, p, d_data, p_data, Kinf_data, Pinf_data);
  tiny_InitSolnDuals(&work, 0, YU, 0, YU_data, 0);

  tiny_SetInitialState(&work, x0_data);  
  tiny_SetGoalReference(&work, Xref, Uref, xg_data, ug_data);

  /* Set up LQR cost */
  tiny_InitDataQuadCost(&work, Q_data, R_data);
  tiny_InitDataLinearCost(&work, q, r, r_tilde, q_data, r_data, r_tilde_data);

  // /* Set up constraints */
  tiny_SetInputBound(&work, Acu_data, umin_data, umax_data);
  (data.ucu).fill(0.5);
  (data.lcu).fill(-0.5);

  tiny_UpdateLinearCost(&work);

  if (0) {
    printf("\nProblem Info: \n");
    PrintMatrix(work.data->model->A[0]);
    PrintMatrix(work.data->model->B[0]);
    PrintMatrix(work.data->Q);
    PrintMatrix(work.data->R);
    PrintMatrixT(work.data->x0);
    PrintMatrixT(work.data->Xref[0]);
    PrintMatrixT(work.data->Uref[0]);
    PrintMatrixT(work.data->q[0]);
    PrintMatrixT(work.data->r[0]);
  }
  // tiny_AddStageCost(&work, 0);
  // printf("%f\n",info.obj_val);
  /* Solver settings */
  stgs.en_cstr_goal = 0;
  stgs.en_cstr_inputs = 1;
  stgs.en_cstr_states = 0;
  stgs.max_iter = 1;           // limit this if needed
  stgs.verbose = 0;
  stgs.check_termination = 1;
  stgs.tol_abs_dual = 5e-2;
  stgs.tol_abs_prim = 5e-2;

  // Absolute formulation:
  // Warm-starting since horizon data is reused
  // Stop earlier as horizon exceeds the end
  X[0] = work.data->x0;  
  srand(1);  // random seed

  /* End of MPC initialization*/
}

int main() {
  init_mpc();

  /* Start MPC loop */

  for (int k = 0; k < NSIM - NHORIZON - 1; ++k) {
    MatrixXf pose = X[k](seq(0,5));
    MatrixXf pose_ref = Xref[0](seq(0,5));
    printf("ex[%d] =  %.4f\n", k, (pose - pose_ref).norm());

    // Inject noise into measurement
    for (int j = 0; j < NSTATES; ++j) {
      X[k](j) += X[k](j) * T_NOISE(0);
    }

    clock_t start, end;
    double cpu_time_used;
    start = clock();

    work.data->x0 = X[k]; // update current measurement
    // Warm-start by previous solution
    // tiny_ShiftFill(Uhrz, T_ARRAY_SIZE(Uhrz));

    // Solve optimization problem using Augmented Lagrangian TVLQR
    tiny_SolveAdmm(&work);

    end = clock();
    cpu_time_used = ((double)(end - start)) * 1000 / CLOCKS_PER_SEC;  // ms
    // printf("solve time:        %f\n", cpu_time_used);
    printf("%f\n", cpu_time_used);

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
    // tiny_QuadNonlinearDynamics(&X[k + 1], X[k], ZU_new[0]);
    tiny_EvalModel(&X[k + 1], &X[k], &ZU_new[0], &model, 0);
  }

  return 0;
}