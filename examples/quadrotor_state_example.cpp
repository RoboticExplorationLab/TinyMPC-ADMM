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
#define NHORIZON 15  // horizon steps (NHORIZON states and NHORIZON-1 controls)
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
  -0.057232f,0.050739f,0.379355f,-0.553259f,-0.653429f,-0.459974f,-0.080869f,0.070832f,0.407553f,-0.057094f,-0.071759f,-0.394156f,
  0.053220f,0.037661f,0.379355f,-0.373387f,0.611098f,0.460681f,0.075315f,0.051340f,0.407553f,-0.033037f,0.067617f,0.394649f,
  0.034565f,-0.041680f,0.379355f,0.415788f,0.309550f,-0.462159f,0.046026f,-0.056905f,0.407553f,0.037185f,0.021970f,-0.395660f,
  -0.030554f,-0.046719f,0.379355f,0.510857f,-0.267219f,0.461452f,-0.040471f,-0.065267f,0.407553f,0.052945f,-0.017827f,0.395167f;

  Pinf << 
  41420.035434f,-4.091413f,0.000000f,44.279630f,58036.704784f,167.018016f,15601.959249f,-5.740568f,0.000000f,3.992407f,294.777299f,128.525720f,
  -4.091413f,41413.045874f,-0.000000f,-57960.011987f,-44.281669f,-66.788621f,-5.740649f,15592.102523f,-0.000000f,-287.896741f,-3.992693f,-51.396030f,
  0.000000f,-0.000000f,32207.127893f,0.000000f,0.000000f,0.000000f,0.000000f,-0.000000f,4376.443528f,0.000000f,0.000000f,0.000000f,
  44.279630f,-57960.011987f,0.000000f,576434.924391f,547.611934f,995.271282f,64.471042f,-78921.519453f,0.000000f,2919.623337f,56.915221f,778.327314f,
  58036.704784f,-44.281669f,0.000000f,547.611934f,577416.579908f,2488.613706f,79034.457834f,-64.473019f,0.000000f,56.914230f,3023.378749f,1946.151255f,
  167.018016f,-66.788621f,0.000000f,995.271282f,2488.613706f,75299.251217f,255.890258f,-102.331134f,0.000000f,131.747777f,329.405832f,5536.068343f,
  15601.959249f,-5.740649f,0.000000f,64.471042f,79034.457834f,255.890258f,18699.995920f,-8.138366f,0.000000f,6.018898f,404.800977f,197.453821f,
  -5.740568f,15592.102523f,-0.000000f,-78921.519453f,-64.473019f,-102.331134f,-8.138366f,18685.925369f,-0.000000f,-394.239464f,-6.019215f,-78.962379f,
  0.000000f,-0.000000f,4376.443528f,0.000000f,0.000000f,0.000000f,0.000000f,-0.000000f,4909.957070f,0.000000f,0.000000f,0.000000f,
  3.992407f,-287.896741f,0.000000f,2919.623337f,56.914230f,131.747777f,6.018898f,-394.239464f,0.000000f,526.171216f,7.805173f,112.487991f,
  294.777299f,-3.992693f,0.000000f,56.915221f,3023.378749f,329.405832f,404.800977f,-6.019215f,0.000000f,7.805173f,541.030271f,281.246887f,
  128.525720f,-51.396030f,0.000000f,778.327314f,1946.151255f,5536.068343f,197.453821f,-78.962379f,0.000000f,112.487991f,281.246887f,4975.759484f;

  Quu_inv << 
  0.000491f,-0.000004f,0.000476f,-0.000003f,
  -0.000004f,0.000489f,-0.000002f,0.000477f,
  0.000476f,-0.000002f,0.000489f,-0.000002f,
  -0.000003f,0.000477f,-0.000002f,0.000488f;

  AmBKt << 
  0.999997f,-0.000000f,0.000000f,0.000011f,-0.005095f,0.000295f,-0.000666f,-0.000001f,0.000000f,0.002128f,-1.018954f,0.058903f,
  -0.000000f,0.999997f,-0.000000f,0.005103f,-0.000011f,-0.000117f,-0.000001f,-0.000668f,-0.000000f,1.020662f,-0.002126f,-0.023495f,
  0.000000f,0.000000f,0.998724f,0.000000f,0.000000f,-0.000000f,0.000000f,-0.000000f,-0.127593f,0.000000f,0.000000f,-0.000000f,
  0.000000f,-0.003889f,-0.000000f,0.946506f,0.000100f,0.001198f,0.000013f,-0.385403f,-0.000000f,-10.698703f,0.019929f,0.239586f,
  0.003889f,-0.000000f,0.000000f,0.000100f,0.946579f,0.003003f,0.385412f,-0.000013f,0.000000f,0.019946f,-10.684270f,0.600617f,
  0.000000f,-0.000000f,-0.000000f,0.000065f,0.000164f,0.999247f,0.000021f,-0.000009f,-0.000000f,0.013047f,0.032716f,-0.150526f,
  0.019995f,-0.000000f,-0.000000f,0.000014f,-0.007042f,0.000404f,0.999079f,-0.000002f,-0.000000f,0.002851f,-1.408367f,0.080851f,
  -0.000000f,0.019995f,-0.000000f,0.007053f,-0.000014f,-0.000161f,-0.000002f,0.999077f,-0.000000f,1.410594f,-0.002849f,-0.032250f,
  0.000000f,0.000000f,0.018629f,0.000000f,0.000000f,-0.000000f,0.000000f,-0.000000f,0.862922f,0.000000f,0.000000f,-0.000000f,
  0.000000f,-0.000010f,-0.000000f,0.004797f,0.000009f,0.000111f,0.000001f,-0.001281f,-0.000000f,-0.040563f,0.001772f,0.022262f,
  0.000010f,-0.000000f,-0.000000f,0.000009f,0.004803f,0.000279f,0.001282f,-0.000001f,-0.000000f,0.001774f,-0.039367f,0.055807f,
  0.000000f,-0.000000f,-0.000000f,0.000050f,0.000126f,0.009356f,0.000016f,-0.000007f,-0.000000f,0.010015f,0.025112f,0.871280f;

  coeff_d2p << 
  0.122765f,-0.134452f,-0.121471f,0.133158f,
  -0.121569f,-0.132709f,0.121036f,0.133243f,
  -0.007238f,-0.007238f,-0.007238f,-0.007238f,
  0.254058f,0.274010f,-0.249267f,-0.278801f,
  0.258748f,-0.283460f,-0.246850f,0.271562f,
  0.111646f,-0.137566f,0.209995f,-0.184076f,
  0.055577f,-0.060867f,-0.054211f,0.059501f,
  -0.054845f,-0.059580f,0.054292f,0.060134f,
  -0.001171f,-0.001171f,-0.001171f,-0.001171f,
  0.001501f,0.001092f,-0.000868f,-0.001724f,
  0.001873f,-0.002085f,-0.000292f,0.000503f,
  0.010028f,-0.011700f,0.016371f,-0.014700f;

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
  stgs.rho_init = 500.0;  // Important (select offline, associated with precomp.)
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
  ucx(0) = 0.35f;  // make it more aggressive
  PrintMatrixT(Acx);
  PrintMatrixT(ucx);
  printf("Is constrained? %d\n", IsConstrained(&work));
  tiny_UpdateLinearCost(&work);

  /* Solver settings */
  stgs.en_cstr_states = 1;
  stgs.max_iter = 50;           // limit this if needed
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
  // for (int k = 0; k < NSIM - NHORIZON - 1; ++k) {
  for (int k = 0; k < 30; ++k) {
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