// MPC
// Scenerio: drive the Crazyflie quadrotor to track fig-8 reference (easy)
//

#include "data/traj_pos_fig8.h"
#include "time.h"
#include "tinympc/tinympc.h"
#include "Eigen.h"

// Macro variables
#define H 0.01       // dt
// #define NSTATES 12   // no. of states (error state)
// #define NINPUTS 4    // no. of controls
#define NHORIZON 10  // horizon steps (NHORIZON states and NHORIZON-1 controls)
#define NSIM 500     // simulation steps (fixed with reference data)

using namespace Eigen;

int main() {
  /* Start MPC initialization*/

  // Create data array 
  float x0_data[NSTATES] = {0, 1, 0, 0.1, 0, 0,
                             0, 0, 0, 0,   0, 0};  // initial state
  float ug_data[NINPUTS] = {0., 0., 0., 0.};      // goal input if needed
  float Xhrz_data[NSTATES * NHORIZON] = {0};      // save X for one horizon
  float X_data[NSTATES * NSIM] = {0};             // save X for the whole run
  float Uhrz_data[NINPUTS * (NHORIZON - 1)] = {0};
  float Xref_data[NSTATES * NHORIZON] = {0};

  float A_data[NSTATES*NSTATES] = {
    1.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,
    0.000000f,1.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,
    0.000000f,0.000000f,1.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,
    0.000000f,-0.000981f,0.000000f,1.000000f,0.000000f,0.000000f,0.000000f,-0.196200f,0.000000f,0.000000f,0.000000f,0.000000f,
    0.000981f,0.000000f,0.000000f,0.000000f,1.000000f,0.000000f,0.196200f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,
    0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,1.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,
    0.010000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,1.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,
    0.000000f,0.010000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,1.000000f,0.000000f,0.000000f,0.000000f,0.000000f,
    0.000000f,0.000000f,0.010000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,1.000000f,0.000000f,0.000000f,0.000000f,
    0.000000f,-0.000002f,0.000000f,0.005000f,0.000000f,0.000000f,0.000000f,-0.000491f,0.000000f,1.000000f,0.000000f,0.000000f,
    0.000002f,0.000000f,0.000000f,0.000000f,0.005000f,0.000000f,0.000491f,0.000000f,0.000000f,0.000000f,1.000000f,0.000000f,
    0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.005000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,1.000000f,
  };

  float B_data[NSTATES*NINPUTS] = {
    -0.000001f,0.000001f,0.000245f,-0.006884f,-0.006918f,0.000494f,-0.000452f,0.000450f,0.049050f,-2.753546f,-2.767070f,0.197477f,
    0.000001f,0.000001f,0.000245f,-0.007581f,0.007607f,-0.000181f,0.000497f,0.000496f,0.049050f,-3.032340f,3.042784f,-0.072236f,
    0.000001f,-0.000001f,0.000245f,0.006892f,0.006939f,-0.000696f,0.000454f,-0.000451f,0.049050f,2.756626f,2.775695f,-0.278438f,
    -0.000001f,-0.000001f,0.000245f,0.007573f,-0.007629f,0.000383f,-0.000499f,-0.000495f,0.049050f,3.029260f,-3.051409f,0.153197f,
  };

  float f_data[NSTATES] = {0};

  static float Kinf_data[NINPUTS*NSTATES] = {
    -0.445949f,0.430873f,0.444882f,-0.429807f,
    0.446807f,0.431253f,-0.446393f,-0.431667f,
    0.662483f,0.662483f,0.662483f,0.662483f,
    -1.079486f,-0.997766f,1.039091f,1.038161f,
    -1.107855f,1.066641f,1.006780f,-0.965566f,
    -0.664215f,0.664103f,-0.663748f,0.663860f,
    -0.226285f,0.218172f,0.220855f,-0.212741f,
    0.225168f,0.214860f,-0.223002f,-0.217025f,
    0.268173f,0.268173f,0.268173f,0.268173f,
    -0.059505f,-0.046573f,0.049043f,0.057035f,
    -0.067296f,0.064829f,0.041135f,-0.038669f,
    -0.214773f,0.214086f,-0.212402f,0.213090f,
  };
  static float Pinf_data[NSTATES*NSTATES] = {
    5012.158703f,-6.325300f,0.000000f,14.273477f,2301.168460f,36.019690f,1181.204851f,-3.165742f,0.000000f,0.603986f,13.813608f,1.702194f,
    -6.325300f,5010.429420f,-0.000000f,-2297.620821f,-14.286189f,-14.399647f,-3.167017f,1180.288583f,-0.000000f,-13.789302f,-0.605188f,-0.686372f,
    0.000000f,-0.000000f,4047.756786f,0.000000f,-0.000000f,-0.000000f,0.000000f,-0.000000f,749.077879f,-0.000000f,-0.000000f,-0.000000f,
    14.273477f,-2297.620821f,-0.000000f,5106.055209f,43.400204f,154.241836f,7.890619f,-1112.350622f,-0.000000f,32.472343f,3.134461f,28.417408f,
    2301.168460f,-14.286189f,0.000000f,43.400204f,5134.931751f,385.641796f,1115.639010f,-7.894015f,0.000000f,3.133074f,36.013817f,71.010338f,
    36.019690f,-14.399647f,-0.000000f,154.241836f,385.641796f,6387.584060f,37.505567f,-14.998374f,-0.000000f,22.314037f,55.785526f,908.134818f,
    1181.204851f,-3.167017f,0.000000f,7.890619f,1115.639010f,37.505567f,469.573259f,-1.645398f,0.000000f,0.390670f,7.052074f,4.192225f,
    -3.165742f,1180.288583f,-0.000000f,-1112.350622f,-7.894015f,-14.998374f,-1.645398f,468.979777f,-0.000000f,-6.864903f,-0.391130f,-1.679634f,
    0.000000f,-0.000000f,749.077879f,0.000000f,-0.000000f,-0.000000f,0.000000f,-0.000000f,299.982271f,-0.000000f,-0.000000f,-0.000000f,
    0.603986f,-13.789302f,-0.000000f,32.472343f,3.133074f,22.314037f,0.390670f,-6.864903f,-0.000000f,2.247327f,0.512500f,7.132266f,
    13.813608f,-0.605188f,-0.000000f,3.134461f,36.013817f,55.785526f,7.052074f,-0.391130f,-0.000000f,0.512500f,3.162347f,17.828607f,
    1.702194f,-0.686372f,-0.000000f,28.417408f,71.010338f,908.134818f,4.192225f,-1.679634f,-0.000000f,7.132266f,17.828607f,290.342050f,
  };
  static float Quu_inv_data[NINPUTS*NINPUTS] = {
    0.012850f,-0.000080f,0.004881f,-0.000096f,
    -0.000080f,0.012580f,-0.000090f,0.005146f,
    0.004881f,-0.000090f,0.012832f,-0.000066f,
    -0.000096f,0.005146f,-0.000066f,0.012572f,
  };
  static float AmBKt_data[NSTATES*NSTATES] = {
    0.999998f,-0.000000f,-0.000000f,0.000000f,0.000976f,0.000000f,0.009999f,-0.000000f,-0.000000f,0.000000f,0.000001f,-0.000000f,
    -0.000000f,0.999998f,-0.000000f,-0.000976f,-0.000000f,-0.000000f,-0.000000f,0.009999f,0.000000f,-0.000001f,-0.000000f,0.000000f,
    -0.000000f,0.000000f,0.999350f,0.000000f,-0.000000f,-0.000000f,-0.000000f,0.000000f,0.009737f,-0.000000f,-0.000000f,-0.000000f,
    0.000386f,0.012690f,0.000000f,0.969982f,0.000834f,0.000009f,0.000185f,0.006359f,0.000000f,0.003467f,0.000038f,-0.000005f,
    -0.012728f,-0.000385f,-0.000000f,0.000832f,0.969870f,0.000024f,-0.006380f,-0.000185f,-0.000000f,0.000037f,0.003461f,-0.000015f,
    0.000772f,-0.000288f,-0.000000f,0.000678f,0.001810f,0.999732f,0.000386f,-0.000144f,0.000000f,0.000033f,0.000088f,0.004915f,
    -0.000832f,-0.000025f,0.000000f,0.000054f,0.194230f,0.000002f,0.999583f,-0.000012f,0.000000f,0.000002f,0.000390f,-0.000001f,
    -0.000025f,-0.000830f,-0.000000f,-0.194237f,-0.000055f,-0.000001f,-0.000012f,0.999584f,0.000000f,-0.000390f,-0.000002f,0.000000f,
    -0.000000f,-0.000000f,-0.129979f,0.000000f,-0.000000f,-0.000000f,-0.000000f,0.000000f,0.947384f,-0.000000f,-0.000000f,-0.000000f,
    0.154237f,5.076178f,0.000000f,-12.007224f,0.333520f,0.003540f,0.074118f,2.543699f,-0.000000f,0.386957f,0.015023f,-0.002197f,
    -5.091399f,-0.154006f,-0.000000f,0.332625f,-12.051921f,0.009416f,-2.552182f,-0.073966f,-0.000000f,0.014964f,0.384352f,-0.005920f,
    0.308906f,-0.115244f,-0.000000f,0.271378f,0.724073f,-0.107374f,0.154532f,-0.057789f,0.000000f,0.013305f,0.035350f,0.966092f,
  };
  static float coeff_d2p_data[NSTATES*NINPUTS] = {
    0.000348f,-0.000320f,0.000033f,0.000410f,0.000664f,0.008163f,0.000106f,-0.000113f,-0.000036f,0.000037f,0.000089f,0.001435f,
    -0.000358f,-0.000297f,0.000033f,0.000073f,-0.000685f,-0.008484f,-0.000109f,-0.000127f,-0.000036f,-0.000035f,-0.000092f,-0.001476f,
    -0.000263f,0.000286f,0.000033f,-0.000051f,0.000236f,0.009318f,-0.000135f,0.000125f,-0.000036f,0.000037f,0.000095f,0.001585f,
    0.000273f,0.000330f,0.000033f,-0.000432f,-0.000215f,-0.008997f,0.000138f,0.000116f,-0.000036f,-0.000039f,-0.000093f,-0.001543f,
  };
  static float Q_data[NSTATES*NSTATES] = {
    100.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,
    0.000000f,100.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,
    0.000000f,0.000000f,100.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,
    0.000000f,0.000000f,0.000000f,100.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,
    0.000000f,0.000000f,0.000000f,0.000000f,100.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,
    0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,100.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,
    0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,1.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,
    0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,1.000000f,0.000000f,0.000000f,0.000000f,0.000000f,
    0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,1.000000f,0.000000f,0.000000f,0.000000f,
    0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,1.000000f,0.000000f,0.000000f,
    0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,1.000000f,0.000000f,
    0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,0.000000f,1.000000f,
  };
  static float R_data[NINPUTS*NINPUTS] = {
    4.000000f,0.000000f,0.000000f,0.000000f,
    0.000000f,4.000000f,0.000000f,0.000000f,
    0.000000f,0.000000f,4.000000f,0.000000f,
    0.000000f,0.000000f,0.000000f,4.000000f,
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
  VectorNf X[NSIM];
  VectorNf Xref[NHORIZON];
  VectorMf Uref[NHORIZON - 1];
  VectorNf Xhrz[NHORIZON];
  VectorMf Uhrz[NHORIZON - 1];
  VectorMf d[NHORIZON - 1];
  VectorNf p[NHORIZON];
  VectorMf YU[NHORIZON - 1];
  VectorMf ZU[NHORIZON - 1];
  VectorMf ZU_new[NHORIZON - 1];
  VectorNf q[NHORIZON-1];
  VectorMf r[NHORIZON-1];
  VectorMf r_tilde[NHORIZON-1];
  MatrixNf A;
  MatrixNMf B;
  VectorNf f;

  for (int i = 0; i < NSIM; ++i) {
    X[i] = Map<VectorNf>(&X_data[i * NSTATES]);
  }

  for (int i = 0; i < NHORIZON; ++i) {
    if (i < NHORIZON - 1) {
      Uref[i] = Map<VectorMf>(ug_data);
    }
    Xref[i] = Map<VectorNf>(&Xref_data[i * NSTATES]);
  }

  for (int i = 0; i < NHORIZON; ++i) {
    for (int j = 0; j < 3; ++j) {
      Xref_data[i*NSTATES + j] = X_ref_data[(i)*3+j];
      printf("%f\n", Xref_data[i*NSTATES + j]);
    }
  }
  PrintMatrix(Xref[0]);
  printf("%d\n", Xref[0].data());
  printf("%d\n\n", &Xref_data[0]);
  printf("%d\n", Uref[0].data());
  printf("%d\n", ug_data);
  /* Create TinyMPC struct and problem data*/
  tiny_Model model;
  tiny_InitModel(&model, NSTATES, NINPUTS, NHORIZON, 0, 0, 0.01);
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
  data.Xref = Xref;
  data.Uref = Uref;

  /* Set up LQR cost */
  tiny_InitDataQuadCostFromArray(&work, Q_data, R_data);
  // slap_AddIdentity(data.R, work.rho); // \tilde{R}
  tiny_InitDataLinearCostFromArray(&work, q, r, r_tilde, q_data, r_data, r_tilde_data);

  /* Set up constraints */
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
    PrintMatrixT(work.data->Xref[NHORIZON-5]);
    PrintMatrixT(work.data->Uref[NHORIZON-5]);
    PrintMatrixT(work.data->q[NHORIZON-5]);
    PrintMatrixT(work.data->r[NHORIZON-5]);
  }

  /* Solver settings */
  stgs.en_cstr_goal = 0;
  stgs.en_cstr_inputs = 1;
  stgs.en_cstr_states = 0;
  stgs.max_iter = 50;           // limit this if needed
  stgs.verbose = 0;
  stgs.check_termination = 2;
  stgs.tol_abs_dual = 1e-2;
  stgs.tol_abs_prim = 1e-2;

  // Absolute formulation:
  // Warm-starting since horizon data is reused
  // Stop earlier as horizon exceeds the end
  X[0] = work.data->x0;  
  srand(1);  // random seed

  /* End of MPC initialization*/

  /* Start MPC loop */

  for (int k = 0; k < 1; ++k) {
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

    // Update reference: 3 positions, k counts each MPC step
    for (int i = 0; i < NHORIZON; ++i) {
      for (int j = 0; j < 3; ++j) {
        Xref_data[i*NSTATES + j] = X_ref_data[(k+i)*3+j];
        // printf("%f\n", Xref_data[i*NSTATES + j]);
      }
    }
    tiny_UpdateLinearCost(&work);

    for (int i = 0; i < NHORIZON; ++i) { 
      // PrintMatrixT(Xref[i]);
    }

    // Warm-start by previous solution
    // tiny_ShiftFill(Uhrz, T_ARRAY_SIZE(Uhrz));

    // Solve optimization problem using Augmented Lagrangian TVLQR
    tiny_SolveAdmm(&work);

    end = clock();
    cpu_time_used = ((double)(end - start)) * 1000 / CLOCKS_PER_SEC;  // ms
    printf("solve time:        %f\n", cpu_time_used);
    // printf("%f\n", cpu_time_used);

    if(work.info->status_val != TINY_SOLVED) {
      printf("!!! STOP AS SOLVER FAILED !!!\n");
      // return 0;
    }

    // PrintMatrixT(Uhrz[0]);

    // Matrix pos = slap_CreateSubMatrix(X[k], 0, 0, 3, 1);
    // PrintMatrixT(pose);
    // PrintMatrixT(pose_ref);

    // === 2. Simulate dynamics using the first control solution ===
    // tiny_Clamp(ZU_new[0].data, umin_data[0], umax_data[0], NINPUTS);
    tiny_EvalModel(&X[k + 1], &X[k], &ZU_new[0], &model, 0);
  }

  return 0;
}
