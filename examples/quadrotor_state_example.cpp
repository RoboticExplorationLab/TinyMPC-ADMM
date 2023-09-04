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
#define NHORIZON 10  // horizon steps (NHORIZON states and NHORIZON-1 controls)
#define NSIM 200     // simulation steps (fixed with reference data)

using namespace Eigen;

// Precomputed data and cache
static MatrixNf A;
static MatrixNMf B;
static MatrixMNf Kinf;
static MatrixNf Pinf;
static MatrixMf Quu_inv;
static MatrixNf AmBKt;
static MatrixNMf coeff_d2p;
static MatrixNf Q;
static MatrixMf R;

static VectorNf Xhrz[NHORIZON];
static VectorMf Uhrz[NHORIZON-1]; 
static VectorMf d[NHORIZON-1];
static VectorNf p[NHORIZON];
static VectorMf YU[NHORIZON-1];
static VectorNf YX[NHORIZON];

static VectorNf q[NHORIZON-1];
static VectorMf r[NHORIZON-1];
static VectorNf q_tilde[NHORIZON];
static VectorMf r_tilde[NHORIZON-1];

static VectorNf Xref[NSIM];
static VectorMf Uref[NSIM-1];

static MatrixMf Acu;
static VectorMf ucu;
static VectorMf lcu;
static VectorNf Acx;
static VectorNf ucx;
static VectorNf lcx;

static VectorMf Qu;
static VectorMf ZU[NHORIZON-1]; 
static VectorMf ZU_new[NHORIZON-1];
static VectorNf ZX[NHORIZON]; 
static VectorNf ZX_new[NHORIZON];

static VectorNf x0 = (VectorNf() << 0.5, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0).finished();
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
  A << 
  1.000000f,0.000000f,0.000000f,0.000000f,0.003924f,0.000000f,0.020000f,0.000000f,0.000000f,0.000000f,0.000013f,0.000000f,
  0.000000f,1.000000f,0.000000f,-0.003924f,0.000000f,0.000000f,0.000000f,0.020000f,0.000000f,-0.000013f,0.000000f,0.000000f,
  0.000000f,0.000000f,1.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.020000f,0.000000f,0.000000f,0.000000f,
  0.000000f,0.000000f,0.000000f,1.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.010000f,0.000000f,0.000000f,
  0.000000f,0.000000f,0.000000f,0.000000f,1.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.010000f,0.000000f,
  0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,1.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.010000f,
  0.000000f,0.000000f,0.000000f,0.000000f,0.392400f,0.000000f,1.000000f,0.000000f,0.000000f,0.000000f,0.001962f,0.000000f,
  0.000000f,0.000000f,0.000000f,-0.392400f,0.000000f,0.000000f,0.000000f,1.000000f,0.000000f,-0.001962f,0.000000f,0.000000f,
  0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,1.000000f,0.000000f,0.000000f,0.000000f,
  0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,1.000000f,0.000000f,0.000000f,
  0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,1.000000f,0.000000f,
  0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,1.000000f;

  B << 
  -0.000018f,0.000020f,0.000018f,-0.000020f,
  0.000018f,0.000020f,-0.000018f,-0.000020f,
  0.000841f,0.000841f,0.000841f,0.000841f,
  -0.027535f,-0.030323f,0.027566f,0.030293f,
  -0.027671f,0.030428f,0.027757f,-0.030514f,
  0.001975f,-0.000722f,-0.002784f,0.001532f,
  -0.003619f,0.003980f,0.003631f,-0.003991f,
  0.003602f,0.003966f,-0.003606f,-0.003962f,
  0.084086f,0.084086f,0.084086f,0.084086f,
  -5.507092f,-6.064681f,5.513253f,6.058520f,
  -5.534140f,6.085568f,5.551390f,-6.102818f,
  0.394954f,-0.144473f,-0.556875f,0.306394f;

  Kinf << 
  -0.073195f,0.070245f,0.276622f,-0.536677f,-0.586364f,-0.476395f,-0.078963f,0.074782f,0.245855f,-0.051780f,-0.061284f,-0.260023f,
  0.067812f,0.060312f,0.276622f,-0.426069f,0.545932f,0.477777f,0.073251f,0.062875f,0.245855f,-0.035092f,0.057581f,0.260461f,
  0.061833f,-0.065711f,0.276622f,0.466578f,0.410927f,-0.480819f,0.063484f,-0.068601f,0.245855f,0.038800f,0.028817f,-0.261396f,
  -0.056450f,-0.064847f,0.276622f,0.496168f,-0.370495f,0.479437f,-0.057772f,-0.069056f,0.245855f,0.048072f,-0.025114f,0.260959f;

  Pinf << 
  7910.584235f,-2.220603f,0.000000f,16.570079f,9455.855819f,82.095002f,2742.912303f,-2.379319f,0.000000f,1.326068f,52.534196f,32.782798f,
  -2.220603f,7907.685455f,-0.000000f,-9433.627861f,-16.574446f,-32.814779f,-2.379543f,2739.758491f,-0.000000f,-50.818248f,-1.326645f,-13.104224f,
  -0.000000f,-0.000000f,6665.817548f,0.000000f,0.000000f,-0.000000f,-0.000000f,-0.000000f,1545.551178f,0.000000f,-0.000000f,-0.000000f,
  16.570079f,-9433.627861f,0.000000f,65547.924008f,156.466643f,459.095923f,19.126554f,-9747.746749f,0.000000f,363.591593f,16.286952f,200.378978f,
  9455.855819f,-16.574446f,0.000000f,156.466643f,65785.164874f,1148.127202f,9774.701909f,-19.129518f,0.000000f,16.285494f,389.885404f,501.094015f,
  82.095002f,-32.814779f,-0.000000f,459.095923f,1148.127202f,24504.557690f,108.809637f,-43.500064f,-0.000000f,68.070635f,170.206741f,2800.530483f,
  2742.912303f,-2.379543f,0.000000f,19.126554f,9774.701909f,108.809637f,2411.129065f,-2.611207f,0.000000f,1.648608f,55.349765f,44.474085f,
  -2.379319f,2739.758491f,-0.000000f,-9747.746749f,-19.129518f,-43.500064f,-2.611207f,2407.574458f,-0.000000f,-53.015893f,-1.649075f,-17.780539f,
  -0.000000f,-0.000000f,1545.551178f,0.000000f,0.000000f,-0.000000f,-0.000000f,0.000000f,1385.193982f,0.000000f,-0.000000f,-0.000000f,
  1.326068f,-50.818248f,-0.000000f,363.591593f,16.285494f,68.070635f,1.648608f,-53.015893f,-0.000000f,58.114161f,2.627363f,36.848883f,
  52.534196f,-1.326645f,0.000000f,16.286952f,389.885404f,170.206741f,55.349765f,-1.649075f,-0.000000f,2.627363f,62.917362f,92.132822f,
  32.782798f,-13.104224f,-0.000000f,200.378978f,501.094015f,2800.530483f,44.474085f,-17.780539f,-0.000000f,36.848883f,92.132822f,1549.011876f;

  Quu_inv << 
  0.001083f,-0.000001f,0.000958f,-0.000000f,
  -0.000001f,0.001072f,0.000001f,0.000968f,
  0.000958f,0.000001f,0.001078f,0.000003f,
  -0.000000f,0.000968f,0.000003f,0.001070f;

  AmBKt << 
  0.999995f,-0.000000f,-0.000000f,0.000046f,-0.007528f,0.000452f,-0.000985f,-0.000006f,-0.000000f,0.009263f,-1.505513f,0.090435f,
  -0.000000f,0.999995f,0.000000f,0.007539f,-0.000046f,-0.000179f,-0.000006f,-0.000986f,0.000000f,1.507781f,-0.009248f,-0.035754f,
  0.000000f,0.000000f,0.999070f,-0.000000f,0.000000f,-0.000000f,0.000000f,0.000000f,-0.093040f,-0.000000f,0.000000f,-0.000000f,
  0.000000f,-0.003888f,0.000000f,0.944410f,0.000303f,0.001291f,0.000040f,-0.385129f,0.000000f,-11.117913f,0.060693f,0.258210f,
  0.003888f,-0.000000f,0.000000f,0.000304f,0.944452f,0.003264f,0.385134f,-0.000040f,0.000000f,0.060844f,-11.109604f,0.652812f,
  0.000000f,-0.000000f,-0.000000f,0.000101f,0.000256f,0.999213f,0.000033f,-0.000013f,-0.000000f,0.020214f,0.051155f,-0.157473f,
  0.019995f,-0.000000f,0.000000f,0.000047f,-0.007939f,0.000474f,0.998962f,-0.000006f,0.000000f,0.009398f,-1.587772f,0.094824f,
  -0.000000f,0.019995f,-0.000000f,0.007949f,-0.000047f,-0.000187f,-0.000006f,0.998960f,-0.000000f,1.589739f,-0.009380f,-0.037496f,
  0.000000f,0.000000f,0.019173f,-0.000000f,0.000000f,-0.000000f,0.000000f,0.000000f,0.917308f,-0.000000f,0.000000f,-0.000000f,
  0.000000f,-0.000010f,0.000000f,0.004984f,0.000025f,0.000111f,0.000003f,-0.001306f,0.000000f,-0.003140f,0.004975f,0.022259f,
  0.000010f,-0.000000f,0.000000f,0.000025f,0.004986f,0.000281f,0.001306f,-0.000003f,0.000000f,0.004990f,-0.002804f,0.056265f,
  0.000000f,-0.000000f,-0.000000f,0.000039f,0.000098f,0.009574f,0.000013f,-0.000005f,-0.000000f,0.007761f,0.019640f,0.914805f;

  coeff_d2p << 
  0.000597f,-0.000646f,-0.000579f,0.000629f,
  -0.000587f,-0.000629f,0.000579f,0.000636f,
  -0.000014f,-0.000014f,-0.000014f,-0.000014f,
  0.001267f,0.001296f,-0.001194f,-0.001369f,
  0.001332f,-0.001435f,-0.001150f,0.001253f,
  0.001612f,-0.001672f,0.001838f,-0.001778f,
  0.000273f,-0.000294f,-0.000252f,0.000274f,
  -0.000264f,-0.000278f,0.000256f,0.000286f,
  -0.000007f,-0.000007f,-0.000007f,-0.000007f,
  0.000014f,0.000000f,0.000001f,-0.000014f,
  0.000024f,-0.000025f,0.000012f,-0.000011f,
  0.000276f,-0.000283f,0.000300f,-0.000294f;

  Q << 
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
  0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,25.000000f;

  R << 
  400.000000f,0.000000f,0.000000f,0.000000f,
  0.000000f,400.000000f,0.000000f,0.000000f,
  0.000000f,0.000000f,400.000000f,0.000000f,
  0.000000f,0.000000f,0.000000f,400.000000f;

  tiny_InitModel(&model, NSTATES, NINPUTS, NHORIZON, 0, 0, DT, &A, &B, 0);
  tiny_InitSettings(&stgs);
  stgs.rho_init = 50.0;  // Important (select offline, associated with precomp.)
  tiny_InitWorkspace(&work, &info, &model, &data, &soln, &stgs);
  
  // Fill in the remaining struct 
  tiny_InitWorkspaceTemp(&work, &Qu, ZU, ZU_new, ZX, ZX_new);
  tiny_InitPrimalCache(&work, &Quu_inv, &AmBKt, &coeff_d2p);
  tiny_InitSolution(&work, Xhrz, Uhrz, YX, YU, 0, &Kinf, d, &Pinf, p);

  tiny_SetInitialState(&work, &x0);  
  tiny_SetGoalState(&work, Xref, &xg);
  tiny_SetGoalInput(&work, Uref, &ug);

  /* Set up LQR cost */
  tiny_InitDataCost(&work, &Q, q, &R, r, q_tilde, r_tilde);

  /* Set up constraints */
  tiny_SetInputBound(&work, &Acu, &lcu, &ucu);
  ucu.fill(0.5);
  lcu.fill(-0.5);
  tiny_SetStateConstraint(&work, &Acx, &lcx, &ucx);
  Acx(1) = -1.0f;
  ucx(0) = 0.3f;
  PrintMatrixT(Acx);
  PrintMatrixT(ucx);
  printf("Is constrained? %d\n", IsConstrained(&work));
  tiny_UpdateLinearCost(&work);

  /* Solver settings */
  stgs.en_cstr_states = 1;
  stgs.max_iter = 100;           // limit this if needed
  stgs.verbose = 0;
  stgs.check_termination = 0;
  stgs.tol_abs_dual = 1e-2;
  stgs.tol_abs_prim = 1e-2;

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
  // for (int k = 0; k < 2; ++k) {
    MatrixXf pose = X[k].head(6);
    MatrixXf pose_ref = Xref[0].head(6);
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
    printf("Start to solve\n");
    tiny_SolveAdmm(&work);

    end = clock();
    cpu_time_used = ((double)(end - start)) * 1000 / CLOCKS_PER_SEC;  // ms
    // printf("solve time:        %f\n", cpu_time_used);
    // printf("%f\n", cpu_time_used);

    // if(work.info->status_val != TINY_SOLVED) {
    //   printf("!!! STOP AS SOLVER FAILED !!!\n");
    //   return 0;
    // }

    PrintMatrixT(Uhrz[0]);

    PrintMatrixT(pose);
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