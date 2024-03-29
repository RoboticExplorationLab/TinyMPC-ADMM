// MPC
// Scenerio: make Crazyflie hovering
//

#include "time.h"
#include "tinympc/tinympc.h"
#include "Eigen.h"

# ifdef __cplusplus
extern "C" {
# endif // ifdef __cplusplus

// Macro variables
#define DT 0.01       // dt
// #define NSTATES 12   // no. of states (error state)
// #define NINPUTS 4    // no. of controls
// These are already defined in `constants.h`
#define NHORIZON 3  // horizon steps (NHORIZON states and NHORIZON-1 controls)
#define NSIM 200     // simulation steps (fixed with reference data)

using namespace Eigen;

// Precomputed data and cache
static MatrixNf A = (Eigen::MatrixNf() << 
1.000000f,0.000000f,0.000000f,0.000000f,0.000981f,0.000000f,0.010000f,0.000000f,0.000000f,0.000000f,0.000002f,0.000000f,
0.000000f,1.000000f,0.000000f,-0.000981f,0.000000f,0.000000f,0.000000f,0.010000f,0.000000f,-0.000002f,0.000000f,0.000000f,
0.000000f,0.000000f,1.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.010000f,0.000000f,0.000000f,0.000000f,
0.000000f,0.000000f,0.000000f,1.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.005000f,0.000000f,0.000000f,
0.000000f,0.000000f,0.000000f,0.000000f,1.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.005000f,0.000000f,
0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,1.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.005000f,
0.000000f,0.000000f,0.000000f,0.000000f,0.196200f,0.000000f,1.000000f,0.000000f,0.000000f,0.000000f,0.000491f,0.000000f,
0.000000f,0.000000f,0.000000f,-0.196200f,0.000000f,0.000000f,0.000000f,1.000000f,0.000000f,-0.000491f,0.000000f,0.000000f,
0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,1.000000f,0.000000f,0.000000f,0.000000f,
0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,1.000000f,0.000000f,0.000000f,
0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,1.000000f,0.000000f,
0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,1.000000f).finished();

static MatrixNMf B = (MatrixNMf() << 
-0.000001f,0.000001f,0.000001f,-0.000001f,
0.000001f,0.000001f,-0.000001f,-0.000001f,
0.000210f,0.000210f,0.000210f,0.000210f,
-0.006884f,-0.007581f,0.006892f,0.007573f,
-0.006918f,0.007607f,0.006939f,-0.007629f,
0.000494f,-0.000181f,-0.000696f,0.000383f,
-0.000452f,0.000497f,0.000454f,-0.000499f,
0.000450f,0.000496f,-0.000451f,-0.000495f,
0.042043f,0.042043f,0.042043f,0.042043f,
-2.753546f,-3.032340f,2.756626f,3.029260f,
-2.767070f,3.042784f,2.775695f,-3.051409f,
0.197477f,-0.072236f,-0.278438f,0.153197f).finished();

static MatrixMNf Kinf = (MatrixMNf() << 
-0.312260f,0.304883f,0.428867f,-0.855749f,-0.931259f,-0.821987f,-0.180993f,0.173195f,0.241446f,-0.058366f,-0.070072f,-0.307631f,
0.304775f,0.283439f,0.428867f,-0.722432f,0.901351f,0.840037f,0.175893f,0.156169f,0.241446f,-0.040561f,0.067789f,0.311377f,
0.278452f,-0.291426f,0.428867f,0.753441f,0.675057f,-0.883715f,0.151735f,-0.161529f,0.241446f,0.042902f,0.031386f,-0.320428f,
-0.270967f,-0.296895f,0.428867f,0.824740f,-0.645149f,0.865664f,-0.146635f,-0.167835f,0.241446f,0.056025f,-0.029104f,0.316681f).finished();

static MatrixNf Pinf = (MatrixNf() << 
5621.717598f,-16.999819f,0.000000f,44.642509f,3293.481980f,238.946972f,1352.895591f,-9.736375f,0.000000f,2.280192f,24.531210f,56.725447f,
-16.999819f,5598.859190f,-0.000000f,-3232.957759f,-44.763451f,-95.176553f,-9.747921f,1339.573926f,0.000000f,-21.672346f,-2.290968f,-22.603260f,
0.000000f,-0.000000f,5623.057189f,0.000000f,-0.000000f,-0.000000f,0.000000f,-0.000000f,1353.708470f,-0.000000f,-0.000000f,-0.000000f,
44.642509f,-3232.957759f,-0.000000f,8044.073369f,170.914880f,606.771943f,29.442129f,-1751.840465f,-0.000000f,59.494802f,13.068636f,168.651920f,
3293.481980f,-44.763451f,-0.000000f,170.914880f,8312.966125f,1519.441300f,1794.885702f,-29.479318f,-0.000000f,13.058986f,81.515019f,422.170302f,
238.946972f,-95.176553f,-0.000000f,606.771943f,1519.441300f,29392.836577f,194.635185f,-77.634157f,-0.000000f,68.858907f,172.288849f,2831.879754f,
1352.895591f,-9.747921f,0.000000f,29.442129f,1794.885702f,194.635185f,589.167197f,-5.894599f,0.000000f,1.720063f,14.691560f,48.626542f,
-9.736375f,1339.573926f,0.000000f,-1751.840465f,-29.479318f,-77.634157f,-5.894599f,580.781214f,0.000000f,-12.180653f,-1.724297f,-19.403118f,
0.000000f,-0.000000f,1353.708470f,0.000000f,-0.000000f,-0.000000f,0.000000f,-0.000000f,757.339155f,-0.000000f,-0.000000f,-0.000000f,
2.280192f,-21.672346f,-0.000000f,59.494802f,13.058986f,68.858907f,1.720063f,-12.180653f,-0.000000f,5.035126f,1.733579f,25.292495f,
24.531210f,-2.290968f,-0.000000f,13.068636f,81.515019f,172.288849f,14.691560f,-1.724297f,-0.000000f,1.733579f,8.325373f,63.260972f,
56.725447f,-22.603260f,-0.000000f,168.651920f,422.170302f,2831.879754f,48.626542f,-19.403118f,-0.000000f,25.292495f,63.260972f,1048.140886f).finished();

static MatrixMf Quu_inv = (MatrixMf() << 
0.005445f,0.000027f,0.001908f,0.000003f,
0.000027f,0.005348f,0.000007f,0.002001f,
0.001908f,0.000007f,0.005428f,0.000039f,
0.000003f,0.002001f,0.000039f,0.005340f).finished();

static MatrixNf AmBKt = (MatrixNf() << 
0.999999f,-0.000000f,0.000000f,0.000294f,-0.008478f,0.000507f,-0.000554f,-0.000019f,0.000000f,0.117602f,-3.391143f,0.202723f,
-0.000000f,0.999999f,0.000000f,0.008504f,-0.000290f,-0.000188f,-0.000019f,-0.000556f,0.000000f,3.401720f,-0.115848f,-0.075393f,
-0.000000f,0.000000f,0.999639f,0.000000f,0.000000f,-0.000000f,0.000000f,0.000000f,-0.072123f,0.000000f,0.000000f,-0.000000f,
0.000000f,-0.000977f,-0.000000f,0.977194f,0.000639f,0.000501f,0.000042f,-0.194709f,-0.000000f,-9.122312f,0.255584f,0.200244f,
0.000977f,-0.000000f,0.000000f,0.000656f,0.977095f,0.001340f,0.194702f,-0.000043f,0.000000f,0.262382f,-9.161845f,0.535809f,
0.000000f,-0.000000f,-0.000000f,0.000244f,0.000660f,0.999611f,0.000043f,-0.000016f,-0.000000f,0.097650f,0.263871f,-0.155672f,
0.009999f,-0.000000f,0.000000f,0.000152f,-0.004762f,0.000283f,0.999689f,-0.000010f,0.000000f,0.060915f,-1.904636f,0.113160f,
-0.000000f,0.009999f,-0.000000f,0.004760f,-0.000149f,-0.000105f,-0.000010f,0.999689f,-0.000000f,1.904151f,-0.059726f,-0.042185f,
0.000000f,0.000000f,0.009797f,0.000000f,0.000000f,-0.000000f,0.000000f,0.000000f,0.959396f,-0.000000f,0.000000f,-0.000000f,
0.000000f,-0.000001f,0.000000f,0.003571f,0.000034f,0.000030f,0.000002f,-0.000397f,0.000000f,0.428313f,0.013789f,0.011959f,
0.000001f,-0.000000f,-0.000000f,0.000036f,0.003560f,0.000080f,0.000396f,-0.000002f,-0.000000f,0.014255f,0.423915f,0.031932f,
0.000000f,-0.000000f,-0.000000f,0.000053f,0.000143f,0.004864f,0.000009f,-0.000003f,-0.000000f,0.021116f,0.057044f,0.945509f).finished();

static MatrixNMf coeff_d2p = (MatrixNMf() << 
0.001658f,-0.001578f,-0.001371f,0.001291f,
-0.001457f,-0.001255f,0.001342f,0.001369f,
-0.005144f,-0.005144f,-0.005144f,-0.005144f,
0.004056f,0.003558f,-0.003805f,-0.003809f,
0.004565f,-0.004337f,-0.003934f,0.003707f,
0.006094f,-0.006809f,0.008656f,-0.007940f,
0.001188f,-0.001123f,-0.001051f,0.000987f,
-0.001055f,-0.000931f,0.001001f,0.000985f,
-0.002049f,-0.002049f,-0.002049f,-0.002049f,
0.000056f,-0.000005f,0.000006f,-0.000056f,
0.000100f,-0.000100f,0.000055f,-0.000054f,
0.001152f,-0.001220f,0.001394f,-0.001326f).finished();

static MatrixNf Q = (MatrixNf() << 
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

static MatrixMf R = (MatrixMf() << 
100.000000f,0.000000f,0.000000f,0.000000f,
0.000000f,100.000000f,0.000000f,0.000000f,
0.000000f,0.000000f,100.000000f,0.000000f,
0.000000f,0.000000f,0.000000f,100.000000f).finished();

static VectorNf Xhrz[NHORIZON];
static VectorMf Uhrz[NHORIZON-1]; 
static VectorMf d[NHORIZON-1];
static VectorNf p[NHORIZON];
static VectorMf YU[NHORIZON];

static VectorNf q[NHORIZON-1];
static VectorMf r[NHORIZON-1];
static VectorMf r_tilde[NHORIZON-1];

static VectorNf Xref[NSIM];
static VectorMf Uref[NSIM-1];

static MatrixMf Acu;
static VectorMf ucu;
static VectorMf lcu;

static VectorMf Qu;
static VectorMf ZU[NHORIZON-1]; 
static VectorMf ZU_new[NHORIZON-1];

static VectorNf x0 = (VectorNf() << 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0).finished();
static VectorNf xg = (VectorNf() << 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0).finished();
static VectorMf ug = (VectorMf() << 0, 0, 0, 0).finished();;

static VectorNf X[NSIM];

// Create TinyMPC struct
static tiny_Model model;
static tiny_AdmmSettings stgs;
static tiny_AdmmData data;
static tiny_AdmmInfo info;
static tiny_AdmmSolution soln;
static tiny_AdmmWorkspace work;

void InitMpc() {
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
  tiny_SetGoalState(&work, Xref, &xg);
  tiny_SetGoalInput(&work, Uref, &ug);

  /* Set up LQR cost */
  tiny_InitDataCost(&work, &Q, q, &R, r, r_tilde);

  // /* Set up constraints */
  tiny_SetInputBound(&work, &Acu, &lcu, &ucu);
  ucu.fill(0.5);
  lcu.fill(-0.5);
  printf("Is constrained? %d\n", IsConstrained(&work));
  tiny_UpdateLinearCost(&work);

  /* Solver settings */
  stgs.max_iter = 1;           // limit this if needed
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

  InitMpc();

  /* Start MPC loop */

  if (0) {
    printf("\nHOVERING QUADROTOR\n");
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

  printf("\n*** Start MPC loop ***\n");
  for (int k = 0; k < NSIM - NHORIZON - 1; ++k) {
    MatrixXf pose = X[k](seq(0,5));
    MatrixXf pose_ref = Xref[0](seq(0,5));
    printf("ex[%d] =  %.4f\n", k, (pose - pose_ref).norm());

    // Inject noise into measurement
    for (int j = 0; j < NSTATES; ++j) {
      X[k](j) += X[k](j) * T_NOISE(1);
    }

    clock_t start, end;
    double cpu_time_used;
    start = clock();

    work.data->x0 = &(X[k]); // update current measurement
    // PrintMatrixT(*(work.data->x0));
    // Warm-start by previous solution
    // tiny_ShiftFill(Uhrz, T_ARRAY_SIZE(Uhrz));

    // Solve optimization problem using Augmented Lagrangian TVLQR
    tiny_SolveAdmm(&work);

    end = clock();
    cpu_time_used = ((double)(end - start)) * 1000 / CLOCKS_PER_SEC;  // ms
    // printf("solve time:        %f\n", cpu_time_used);
    // printf("%f\n", cpu_time_used);

    // if(work.info->status_val != TINY_SOLVED) {
    //   printf("!!! STOP AS SOLVER FAILED !!!\n");
    //   return 0;
    // }

    // PrintMatrixT(Uhrz[0]);

    // PrintMatrixT(pose);
    // PrintMatrixT(work.data->Xref[0]);

    // === 2. Simulate dynamics using the first control solution ===
    // tiny_Clamp(ZU_new[0].data, umin[0], umax[0], NINPUTS);
    // If no constraints, use Uhrz[0]
    tiny_EvalModel(&X[k + 1], &X[k], ZU_new, &model, 0);
  }

  return 0;
}

# ifdef __cplusplus
}
# endif // ifdef __cplusplus