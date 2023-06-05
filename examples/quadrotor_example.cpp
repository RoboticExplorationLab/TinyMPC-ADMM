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
#define NSIM 200     // simulation steps (fixed with reference data)

using namespace Eigen;

// Create data, can assign values here
static MatrixNf A = (Eigen::MatrixNf() << 
1.000000f,0.000000f,0.000000f,0.000000f,0.000841f,0.000000f,0.010000f,0.000000f,0.000000f,0.000000f,0.000001f,0.000000f,
0.000000f,1.000000f,0.000000f,-0.000841f,0.000000f,0.000000f,0.000000f,0.010000f,0.000000f,-0.000001f,0.000000f,0.000000f,
0.000000f,0.000000f,1.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.010000f,0.000000f,0.000000f,0.000000f,
0.000000f,0.000000f,0.000000f,1.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.005000f,0.000000f,0.000000f,
0.000000f,0.000000f,0.000000f,0.000000f,1.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.005000f,0.000000f,
0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,1.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.005000f,
0.000000f,0.000000f,0.000000f,0.000000f,0.168171f,0.000000f,1.000000f,0.000000f,0.000000f,0.000000f,0.000420f,0.000000f,
0.000000f,0.000000f,0.000000f,-0.168171f,0.000000f,0.000000f,0.000000f,1.000000f,0.000000f,-0.000420f,0.000000f,0.000000f,
0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,1.000000f,0.000000f,0.000000f,0.000000f,
0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,1.000000f,0.000000f,0.000000f,
0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,1.000000f,0.000000f,
0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,1.000000f).finished();

static MatrixNMf B = (Eigen::MatrixNMf() << 
-0.000001f,0.000001f,0.000001f,-0.000001f,
0.000001f,0.000001f,-0.000001f,-0.000001f,
0.000210f,0.000210f,0.000210f,0.000210f,
-0.006884f,-0.007581f,0.006892f,0.007573f,
-0.006918f,0.007607f,0.006939f,-0.007629f,
0.000494f,-0.000181f,-0.000696f,0.000383f,
-0.000388f,0.000426f,0.000389f,-0.000428f,
0.000386f,0.000425f,-0.000386f,-0.000425f,
0.042043f,0.042043f,0.042043f,0.042043f,
-2.753546f,-3.032340f,2.756626f,3.029260f,
-2.767070f,3.042784f,2.775695f,-3.051409f,
0.197477f,-0.072236f,-0.278438f,0.153197f).finished();

static MatrixMNf Kinf = (Eigen::MatrixMNf() << 
-0.314445f,0.306110f,0.428867f,-0.809923f,-0.884296f,-0.818635f,-0.190177f,0.181338f,0.241446f,-0.058038f,-0.069714f,-0.306873f,
0.306995f,0.283354f,0.428867f,-0.680051f,0.856478f,0.838550f,0.184898f,0.162637f,0.241446f,-0.040300f,0.067488f,0.311041f,
0.277270f,-0.291312f,0.428867f,0.708923f,0.631356f,-0.886749f,0.157200f,-0.168189f,0.241446f,0.042585f,0.031055f,-0.321113f,
-0.269820f,-0.298153f,0.428867f,0.781051f,-0.603538f,0.866834f,-0.151921f,-0.175786f,0.241446f,0.055753f,-0.028829f,0.316945f).finished();

static MatrixNf Pinf = (Eigen::MatrixNf() << 
5859.874747f,-17.678697f,0.000000f,44.138684f,3273.627695f,260.971826f,1488.601283f,-10.589332f,0.000000f,2.413299f,24.851934f,62.230614f,
-17.678697f,5835.614041f,-0.000000f,-3212.351189f,-44.250231f,-103.947495f,-10.601268f,1473.806239f,-0.000000f,-21.712125f,-2.423730f,-24.796499f,
0.000000f,-0.000000f,5623.057189f,-0.000000f,-0.000000f,-0.000000f,0.000000f,0.000000f,1353.708470f,-0.000000f,-0.000000f,-0.000000f,
44.138684f,-3212.351189f,-0.000000f,7560.040376f,159.191176f,596.889441f,30.404704f,-1816.632359f,-0.000000f,56.273483f,12.762304f,166.420091f,
3273.627695f,-44.250231f,-0.000000f,159.191176f,7814.149825f,1494.826895f,1861.974084f,-30.440053f,-0.000000f,12.753600f,77.976971f,416.613040f,
260.971826f,-103.947495f,-0.000000f,596.889441f,1494.826895f,29390.506140f,218.160793f,-87.012755f,-0.000000f,68.774961f,172.093004f,2831.348002f,
1488.601283f,-10.601268f,-0.000000f,30.404704f,1861.974084f,218.160793f,673.739858f,-6.706953f,0.000000f,1.890693f,15.557089f,54.737563f,
-10.589332f,1473.806239f,0.000000f,-1816.632359f,-30.440053f,-87.012755f,-6.706953f,664.003392f,-0.000000f,-12.732963f,-1.894930f,-21.840459f,
0.000000f,-0.000000f,1353.708470f,-0.000000f,-0.000000f,-0.000000f,0.000000f,-0.000000f,757.339155f,-0.000000f,-0.000000f,-0.000000f,
2.413299f,-21.712125f,-0.000000f,56.273483f,12.753600f,68.774961f,1.890693f,-12.732963f,-0.000000f,5.013491f,1.729803f,25.273475f,
24.851934f,-2.423730f,-0.000000f,12.762304f,77.976971f,172.093004f,15.557089f,-1.894930f,-0.000000f,1.729803f,8.299841f,63.216602f,
62.230614f,-24.796499f,-0.000000f,166.420091f,416.613040f,2831.348002f,54.737563f,-21.840459f,-0.000000f,25.273475f,63.216602f,1048.018795f).finished();

static MatrixMf Quu_inv = (Eigen::MatrixMf() << 
0.005454f,0.000027f,0.001900f,0.000003f,
0.000027f,0.005356f,0.000007f,0.001993f,
0.001900f,0.000007f,0.005436f,0.000039f,
0.000003f,0.001993f,0.000039f,0.005348f).finished();

static MatrixNf AmBKt = (Eigen::MatrixNf() << 
0.999999f,-0.000000f,-0.000000f,0.000295f,-0.008493f,0.000507f,-0.000476f,-0.000017f,-0.000000f,0.118098f,-3.397160f,0.202810f,
-0.000000f,0.999999f,-0.000000f,0.008521f,-0.000291f,-0.000189f,-0.000016f,-0.000478f,-0.000000f,3.408336f,-0.116350f,-0.075417f,
0.000000f,0.000000f,0.999639f,-0.000000f,0.000000f,-0.000000f,0.000000f,0.000000f,-0.072123f,-0.000000f,0.000000f,-0.000000f,
0.000000f,-0.000838f,-0.000000f,0.978469f,0.000609f,0.000471f,0.000034f,-0.166964f,-0.000000f,-8.612550f,0.243687f,0.188553f,
0.000838f,-0.000000f,-0.000000f,0.000625f,0.978382f,0.001262f,0.166960f,-0.000035f,-0.000000f,0.250042f,-8.647076f,0.504750f,
0.000000f,-0.000000f,-0.000000f,0.000268f,0.000724f,0.999606f,0.000041f,-0.000015f,-0.000000f,0.107190f,0.289663f,-0.157465f,
0.009999f,-0.000000f,-0.000000f,0.000160f,-0.004972f,0.000295f,0.999721f,-0.000009f,-0.000000f,0.063878f,-1.988747f,0.117956f,
-0.000000f,0.009999f,0.000000f,0.004972f,-0.000157f,-0.000110f,-0.000009f,0.999721f,0.000000f,1.988630f,-0.062645f,-0.043962f,
-0.000000f,0.000000f,0.009797f,0.000000f,-0.000000f,0.000000f,-0.000000f,-0.000000f,0.959396f,-0.000000f,-0.000000f,0.000000f,
0.000000f,-0.000001f,-0.000000f,0.003579f,0.000035f,0.000030f,0.000002f,-0.000341f,-0.000000f,0.431706f,0.013950f,0.011866f,
0.000001f,-0.000000f,-0.000000f,0.000036f,0.003569f,0.000079f,0.000340f,-0.000002f,-0.000000f,0.014408f,0.427576f,0.031706f,
0.000000f,-0.000000f,-0.000000f,0.000058f,0.000157f,0.004863f,0.000009f,-0.000003f,-0.000000f,0.023274f,0.062874f,0.945104f).finished();

static MatrixNMf coeff_d2p = (Eigen::MatrixNMf() << 
0.004006f,-0.003878f,-0.003854f,0.003726f,
-0.003726f,-0.003522f,0.003666f,0.003582f,
-0.005144f,-0.005144f,-0.005144f,-0.005144f,
0.006674f,0.006236f,-0.006519f,-0.006391f,
0.007212f,-0.006957f,-0.006819f,0.006563f,
0.005378f,-0.006494f,0.009382f,-0.008265f,
0.002474f,-0.002381f,-0.002409f,0.002316f,
-0.002295f,-0.002167f,0.002270f,0.002193f,
-0.002049f,-0.002049f,-0.002049f,-0.002049f,
0.000071f,0.000013f,-0.000010f,-0.000074f,
0.000112f,-0.000114f,0.000041f,-0.000039f,
0.001054f,-0.001176f,0.001492f,-0.001370f).finished();

static MatrixNf Q = (Eigen::MatrixNf() << 
100.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,
0.000000f,100.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,
0.000000f,0.000000f,100.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,
0.000000f,0.000000f,0.000000f,4.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,
0.000000f,0.000000f,0.000000f,0.000000f,4.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,
0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,400.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,
0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,4.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,
0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,4.000000f,0.000000f,0.000000f,0.000000f,0.000000f,
0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,4.000000f,0.000000f,0.000000f,0.000000f,
0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,2.040816f,0.000000f,0.000000f,
0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,2.040816f,0.000000f,
0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,25.000000f).finished();

static MatrixMf R = (Eigen::MatrixMf() << 
100.000000f,0.000000f,0.000000f,0.000000f,
0.000000f,100.000000f,0.000000f,0.000000f,
0.000000f,0.000000f,100.000000f,0.000000f,
0.000000f,0.000000f,0.000000f,100.000000f).finished();

static VectorNf  Xhrz[NHORIZON];
static VectorMf  Uhrz[NHORIZON-1]; 
static VectorMf  d[NHORIZON-1];
static VectorNf  p[NHORIZON];
static VectorMf  YU[NHORIZON];

static VectorNf x0 = (Eigen::VectorNf() << 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0).finished();

static VectorNf q[NHORIZON-1];
static VectorMf r[NHORIZON-1];
static VectorMf r_tilde[NHORIZON-1];

static VectorNf Xref[NSIM];
static VectorMf Uref[NSIM-1];

static MatrixMf Acu;
static VectorMf ucu;
static VectorMf lcu;

static VectorMf  Qu;
static VectorMf  ZU[NHORIZON-1]; 
static VectorMf  ZU_new[NHORIZON-1];

static VectorNf xg = (Eigen::VectorNf() << 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0).finished();
static VectorMf ug = (Eigen::VectorMf() << 0, 0, 0, 0).finished();;

static VectorNf X[NSIM];

// Create TinyMPC struct
static tiny_Model model;
static tiny_AdmmSettings stgs;
static tiny_AdmmData data;
static tiny_AdmmInfo info;
static tiny_AdmmSolution soln;
static tiny_AdmmWorkspace work;

void init_mpc() {
  /* Start MPC initialization*/

  tiny_InitModel(&model, NSTATES, NINPUTS, NHORIZON, 0, 0, DT, &A, &B, 0);
  tiny_InitSettings(&stgs);
  stgs.rho_init = 30.0;  // Important (select offline, associated with precomp.)
  tiny_InitWorkspace(&work, &info, &model, &data, &soln, &stgs);
  
  // Fill in the remaining struct 
  tiny_InitWorkspaceTemp(&work, &Qu, ZU, ZU_new, 0, 0);
  tiny_InitPrimalCache(&work, &Quu_inv, &AmBKt, &coeff_d2p);
  tiny_InitSolution(&work, Xhrz, Uhrz, 0, YU, 0, &Kinf, d, &Pinf, p);

  tiny_SetInitialState(&work, &x0);  
  tiny_SetGoalReference(&work, Xref, Uref, &xg, &ug);

  /* Set up LQR cost */
  tiny_InitDataCost(&work, &Q, q, &R, r, r_tilde);

  // /* Set up constraints */
  tiny_SetInputBound(&work, &Acu, &lcu, &ucu);
  ucu.fill(0.5);
  lcu.fill(-0.5);

  tiny_UpdateLinearCost(&work);

  if (1) {
    printf("\nProblem Info: \n");
    PrintMatrix(work.data->model->A[0]);
    PrintMatrix(work.data->model->B[0]);
    PrintMatrix((*(work.data->Q)));
    PrintMatrix((*(work.data->R)));
    PrintMatrixT((*(work.data->x0)));
    PrintMatrixT(work.data->Xref[0]);
    PrintMatrixT(work.data->Uref[0]);
    PrintMatrixT(work.data->q[0]);
    PrintMatrixT(work.data->r[0]);
  }

  /* Solver settings */
  stgs.max_iter = 100;           // limit this if needed
  stgs.verbose = 0;
  stgs.check_termination = 1;
  stgs.tol_abs_dual = 5e-2;
  stgs.tol_abs_prim = 5e-2;

  // Absolute formulation:
  // Warm-starting since horizon data is reused
  // Stop earlier as horizon exceeds the end
  X[0] = x0;  
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

    work.data->x0 = &(X[k]); // update current measurement
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
    // tiny_Clamp(ZU_new[0].data, umin[0], umax[0], NINPUTS);
    // If no constraints, use Uhrz[0]
    // tiny_QuadNonlinearDynamics(&X[k + 1], X[k], ZU_new[0]);
    tiny_EvalModel(&X[k + 1], &X[k], &ZU_new[0], &model, 0);
  }

  return 0;
}