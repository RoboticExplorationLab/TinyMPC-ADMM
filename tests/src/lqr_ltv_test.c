// Test LQR
// Scenerio: Drive planar quadrotor to arbitrary goal state.

#include "tinympc/lqr_ltv.h"

#include "planar_quadrotor.h"
#include "simpletest.h"
#include "slap/slap.h"
#include "test_utils.h"
#include "tinympc/lqr_lti.h"
#include "tinympc/utils.h"

#define H 0.1
#define NSTATES 6
#define NINPUTS 2
#define NHORIZON 40
// GRADIENT EXPLOSION/VANISHING WHEN NHORIZON > 60 => LS FORMULATION

float x0_data[NSTATES] = {1, 1, 0.1, 0.1, 0.1, 0.1};  // initial state
float ug_data[NINPUTS] = {0};                         // initial state
float xg_data[NSTATES] = {3, 2, 0, 0, 0, 0};          // goal state
float xhover_data[NSTATES] = {0, 0, 0, 0, 0, 0};
float uhover_data[NINPUTS] = {4.905, 4.905};
float Q_data[NSTATES * NSTATES] = {0};
float R_data[NINPUTS * NINPUTS] = {0};
float Qf_data[NSTATES * NSTATES] = {0};

void DeltaLqrLtvTest() {
  float A_data[NSTATES * NSTATES * (NHORIZON - 1)] = {0};
  float B_data[NSTATES * NINPUTS * (NHORIZON - 1)] = {0};
  float f_data[NSTATES * (NHORIZON - 1)] = {0};
  float X_data[NSTATES * NHORIZON] = {0};
  float U_data[NINPUTS * (NHORIZON - 1)] = {0};
  float K_data[NINPUTS * NSTATES * (NHORIZON - 1)] = {0};
  float d_data[NINPUTS * (NHORIZON - 1)] = {0};
  float P_data[NSTATES * NSTATES * (NHORIZON)] = {0};
  float p_data[NSTATES * NHORIZON] = {0};
  tiny_LtvModel model;
  tiny_InitLtvModel(&model);
  tiny_ProblemData prob;
  tiny_InitProblemData(&prob);
  tiny_AdmmSettings solver;
  tiny_InitSettings(&solver);

  Matrix X[NHORIZON];
  Matrix U[NHORIZON - 1];
  Matrix Xref[NHORIZON];
  Matrix Uref[NHORIZON - 1];
  Matrix K[NHORIZON - 1];
  Matrix d[NHORIZON - 1];
  Matrix P[NHORIZON];
  Matrix p[NHORIZON];
  Matrix A[NHORIZON - 1];
  Matrix B[NHORIZON - 1];
  Matrix f[NHORIZON - 1];

  float* Xptr = X_data;
  // float* Xref_ptr = Xref_data;
  float* Uptr = U_data;
  // float* Uref_ptr = Uref_data;
  float* Kptr = K_data;
  float* dptr = d_data;
  float* Pptr = P_data;
  float* pptr = p_data;
  float* Aptr = A_data;
  float* Bptr = B_data;
  float* fptr = f_data;

  Matrix xhover = slap_Matrix(NSTATES, 1, xhover_data);
  Matrix uhover = slap_Matrix(NINPUTS, 1, uhover_data);

  for (int i = 0; i < NHORIZON; ++i) {
    if (i < NHORIZON - 1) {
      A[i] = slap_Matrix(NSTATES, NSTATES, Aptr);
      Aptr += NSTATES * NSTATES;
      B[i] = slap_Matrix(NSTATES, NINPUTS, Bptr);
      Bptr += NSTATES * NINPUTS;
      f[i] = slap_Matrix(NSTATES, 1, fptr);
      fptr += NSTATES;
      U[i] = slap_Matrix(NINPUTS, 1, Uptr);
      // slap_SetConst(U[i], 0.01);
      Uptr += NINPUTS;
      Uref[i] = slap_Matrix(NINPUTS, 1, ug_data);
      // Uref_ptr += NINPUTS;
      K[i] = slap_Matrix(NINPUTS, NSTATES, Kptr);
      Kptr += NINPUTS * NSTATES;
      d[i] = slap_Matrix(NINPUTS, 1, dptr);
      dptr += NINPUTS;
    }
    X[i] = slap_Matrix(NSTATES, 1, Xptr);
    Xptr += NSTATES;
    Xref[i] = slap_Matrix(NSTATES, 1, xg_data);
    // Xref_ptr += NSTATES;
    P[i] = slap_Matrix(NSTATES, NSTATES, Pptr);
    Pptr += NSTATES * NSTATES;
    p[i] = slap_Matrix(NSTATES, 1, pptr);
    pptr += NSTATES;
  }

  model.ninputs = NSTATES;
  model.nstates = NINPUTS;
  model.x0 = slap_Matrix(NSTATES, 1, x0_data);
  model.get_jacobians = tiny_PQuadGetJacobians;  // from planar quadrotor
  model.A = A;
  model.B = B;
  model.f = f;
  slap_Copy(X[0], model.x0);

  prob.ninputs = NINPUTS;
  prob.nstates = NSTATES;
  prob.nhorizon = NHORIZON;
  prob.Q = slap_Matrix(NSTATES, NSTATES, Q_data);
  slap_SetIdentity(prob.Q, 10e-1);
  prob.R = slap_Matrix(NINPUTS, NINPUTS, R_data);
  slap_SetIdentity(prob.R, 1e-1);
  prob.Qf = slap_Matrix(NSTATES, NSTATES, Qf_data);
  slap_SetIdentity(prob.Qf, 10e-1);
  prob.Xref = Xref;
  prob.Uref = Uref;
  prob.x0 = model.x0;
  prob.K = K;
  prob.d = d;
  prob.P = P;
  prob.p = p;

  float Q_temp_data[(NSTATES + NINPUTS) * (NSTATES + NINPUTS + 1)] = {0};
  Matrix Q_temp = slap_Matrix(NSTATES + NINPUTS, NSTATES + NINPUTS + 1,
                                       Q_temp_data);

  // Compute and store A, B offline
  for (int i = 0; i < NHORIZON - 1; ++i) {
    model.get_jacobians(&(model.A[i]), &(model.B[i]), xhover, uhover);
    tiny_PQuadNonlinearDynamics(&(model.f[i]), xhover, uhover);
    slap_MatrixAddition(model.f[i], model.f[i], xhover, -1);  // = 0
    // PrintMatrix(f[i]);
  }

  tiny_BackwardPassLtv(&prob, solver, model, &Q_temp);
  tiny_ForwardPassLtv(X, U, prob, model);
  // for (int k = 0; k < NHORIZON - 1; ++k) {
  //   // Control input: du = - d - K*dx
  //   slap_Copy(U[k], prob.d[k]);              // du[k] = -d[k]
  //   slap_MatMulAdd(U[k], prob.Kinf, X[k], -1, -1);   // du[k] -= Kinf * dx[k]
  //   // Next state: dx = A*dx + B*du + f
  //   tiny_DynamicsLtv(&X[k + 1], X[k], U[k], model, k);
  //   // tiny_PQuadNonlinearDynamics(&X[k+1], X[k], U[k]);
  // }
  for (int k = 0; k < NHORIZON - 1; ++k) {
    // printf("ex[%d] = %.4f\n", k, slap_NormedDifference(X[k], Xref[k]));
    // PrintMatrix(X[k]);
  }

  for (int k = NHORIZON - 5; k < NHORIZON; ++k) {
    TEST(SumOfSquaredError(X[k].data, Xref[k].data, NSTATES) < 0.5);
  }
}

void AbsLqrLtvTest() {
  float A_data[NSTATES * NSTATES * (NHORIZON - 1)] = {0};
  float B_data[NSTATES * NINPUTS * (NHORIZON - 1)] = {0};
  float f_data[NSTATES * (NHORIZON - 1)] = {0};
  float X_data[NSTATES * NHORIZON] = {0};
  float U_data[NINPUTS * (NHORIZON - 1)] = {0};
  float K_data[NINPUTS * NSTATES * (NHORIZON - 1)] = {0};
  float d_data[NINPUTS * (NHORIZON - 1)] = {0};
  float P_data[NSTATES * NSTATES * (NHORIZON)] = {0};
  float p_data[NSTATES * NHORIZON] = {0};
  tiny_LtvModel model;
  tiny_InitLtvModel(&model);
  tiny_ProblemData prob;
  tiny_InitProblemData(&prob);
  tiny_AdmmSettings solver;
  tiny_InitSettings(&solver);

  Matrix X[NHORIZON];
  Matrix U[NHORIZON - 1];
  Matrix Xref[NHORIZON];
  Matrix Uref[NHORIZON - 1];
  Matrix K[NHORIZON - 1];
  Matrix d[NHORIZON - 1];
  Matrix P[NHORIZON];
  Matrix p[NHORIZON];
  Matrix A[NHORIZON - 1];
  Matrix B[NHORIZON - 1];
  Matrix f[NHORIZON - 1];

  float* Xptr = X_data;
  // float* Xref_ptr = Xref_data;
  float* Uptr = U_data;
  // float* Uref_ptr = Uref_data;
  float* Kptr = K_data;
  float* dptr = d_data;
  float* Pptr = P_data;
  float* pptr = p_data;
  float* Aptr = A_data;
  float* Bptr = B_data;
  float* fptr = f_data;

  Matrix xhover = slap_Matrix(NSTATES, 1, xhover_data);
  Matrix uhover = slap_Matrix(NINPUTS, 1, uhover_data);

  for (int i = 0; i < NHORIZON; ++i) {
    if (i < NHORIZON - 1) {
      A[i] = slap_Matrix(NSTATES, NSTATES, Aptr);
      Aptr += NSTATES * NSTATES;
      B[i] = slap_Matrix(NSTATES, NINPUTS, Bptr);
      Bptr += NSTATES * NINPUTS;
      f[i] = slap_Matrix(NSTATES, 1, fptr);
      fptr += NSTATES;
      U[i] = slap_Matrix(NINPUTS, 1, Uptr);
      // slap_SetConst(U[i], 0.01);
      Uptr += NINPUTS;
      Uref[i] = slap_Matrix(NINPUTS, 1, uhover_data);
      // Uref_ptr += NINPUTS;
      K[i] = slap_Matrix(NINPUTS, NSTATES, Kptr);
      Kptr += NINPUTS * NSTATES;
      d[i] = slap_Matrix(NINPUTS, 1, dptr);
      dptr += NINPUTS;
    }
    X[i] = slap_Matrix(NSTATES, 1, Xptr);
    Xptr += NSTATES;
    Xref[i] = slap_Matrix(NSTATES, 1, xg_data);
    // Xref_ptr += NSTATES;
    P[i] = slap_Matrix(NSTATES, NSTATES, Pptr);
    Pptr += NSTATES * NSTATES;
    p[i] = slap_Matrix(NSTATES, 1, pptr);
    pptr += NSTATES;
  }

  model.ninputs = NSTATES;
  model.nstates = NINPUTS;
  model.x0 = slap_Matrix(NSTATES, 1, x0_data);
  model.get_jacobians = tiny_PQuadGetJacobians;  // from planar quadrotor
  model.A = A;
  model.B = B;
  model.f = f;
  slap_Copy(X[0], model.x0);

  prob.ninputs = NINPUTS;
  prob.nstates = NSTATES;
  prob.nhorizon = NHORIZON;
  prob.Q = slap_Matrix(NSTATES, NSTATES, Q_data);
  slap_SetIdentity(prob.Q, 10e-1);
  prob.R = slap_Matrix(NINPUTS, NINPUTS, R_data);
  slap_SetIdentity(prob.R, 1e-1);
  prob.Qf = slap_Matrix(NSTATES, NSTATES, Qf_data);
  slap_SetIdentity(prob.Qf, 10e-1);
  prob.Xref = Xref;
  prob.Uref = Uref;
  prob.x0 = model.x0;
  prob.K = K;
  prob.d = d;
  prob.P = P;
  prob.p = p;

  float Q_temp_data[(NSTATES + NINPUTS) * (NSTATES + NINPUTS + 1)] = {0};
  Matrix Q_temp = slap_Matrix(NSTATES + NINPUTS, NSTATES + NINPUTS + 1,
                                       Q_temp_data);

  // Compute and store A, B offline
  for (int i = 0; i < NHORIZON - 1; ++i) {
    model.get_jacobians(&(model.A[i]), &(model.B[i]), xhover, uhover);
    tiny_PQuadNonlinearDynamics(&(model.f[i]), xhover, uhover);
    slap_MatMulAdd(model.f[i], model.A[i], xhover, -1, 1);
    slap_MatMulAdd(model.f[i], model.B[i], uhover, -1, 1);
    // PrintMatrix(f[i]);
  }

  tiny_BackwardPassLtv(&prob, solver, model, &Q_temp);
  tiny_ForwardPassLtv(X, U, prob, model);
  // for (int k = 0; k < NHORIZON - 1; ++k) {
  //   // Control input: du = - d - K*(x - xf)
  //   slap_Copy(U[k], prob.d[k]);              // u[k] = -d[k]
  //   slap_MatMulAdd(U[k], prob.Kinf, X[k], -1, -1);   // u[k] -= Kinf * x[k]
  //   tiny_PQuadNonlinearDynamics(&X[k+1], X[k], U[k]);
  // }
  for (int k = 0; k < NHORIZON - 1; ++k) {
    // printf("ex[%d] = %.4f\n", k, slap_NormedDifference(X[k], Xref[k]));
    // PrintMatrix(X[k]);
  }

  for (int k = NHORIZON - 5; k < NHORIZON; ++k) {
    TEST(SumOfSquaredError(X[k].data, Xref[k].data, NSTATES) < 0.5);
  }
}

void DeltaLqrLtiTest() {
  float A_data[NSTATES * NSTATES] = {0};
  float B_data[NSTATES * NINPUTS] = {0};
  float f_data[NSTATES] = {0};
  float X_data[NSTATES * NHORIZON] = {0};
  float U_data[NINPUTS * (NHORIZON - 1)] = {0};
  float K_data[NINPUTS * NSTATES * (NHORIZON - 1)] = {0};
  float d_data[NINPUTS * (NHORIZON - 1)] = {0};
  float P_data[NSTATES * NSTATES * (NHORIZON)] = {0};
  float p_data[NSTATES * NHORIZON] = {0};
  tiny_LtiModel model;
  tiny_InitLtiModel(&model);
  tiny_ProblemData prob;
  tiny_InitProblemData(&prob);
  tiny_AdmmSettings solver;
  tiny_InitSettings(&solver);

  Matrix X[NHORIZON];
  Matrix U[NHORIZON - 1];
  Matrix Xref[NHORIZON];
  Matrix Uref[NHORIZON - 1];
  Matrix K[NHORIZON - 1];
  Matrix d[NHORIZON - 1];
  Matrix P[NHORIZON];
  Matrix p[NHORIZON];

  float* Xptr = X_data;
  // float* Xref_ptr = Xref_data;
  float* Uptr = U_data;
  // float* Uref_ptr = Uref_data;
  float* Kptr = K_data;
  float* dptr = d_data;
  float* Pptr = P_data;
  float* pptr = p_data;

  Matrix xhover = slap_Matrix(NSTATES, 1, xhover_data);
  Matrix uhover = slap_Matrix(NINPUTS, 1, uhover_data);

  for (int i = 0; i < NHORIZON; ++i) {
    if (i < NHORIZON - 1) {
      U[i] = slap_Matrix(NINPUTS, 1, Uptr);
      // slap_SetConst(U[i], 0.01);
      Uptr += NINPUTS;
      Uref[i] = slap_Matrix(NINPUTS, 1, ug_data);
      // Uref_ptr += NINPUTS;
      K[i] = slap_Matrix(NINPUTS, NSTATES, Kptr);
      Kptr += NINPUTS * NSTATES;
      d[i] = slap_Matrix(NINPUTS, 1, dptr);
      dptr += NINPUTS;
    }
    X[i] = slap_Matrix(NSTATES, 1, Xptr);
    Xptr += NSTATES;
    Xref[i] = slap_Matrix(NSTATES, 1, xg_data);
    // Xref_ptr += NSTATES;
    P[i] = slap_Matrix(NSTATES, NSTATES, Pptr);
    Pptr += NSTATES * NSTATES;
    p[i] = slap_Matrix(NSTATES, 1, pptr);
    pptr += NSTATES;
  }

  model.ninputs = NSTATES;
  model.nstates = NINPUTS;
  model.x0 = slap_Matrix(NSTATES, 1, x0_data);
  model.f = slap_Matrix(NSTATES, 1, f_data);
  model.A = slap_Matrix(NSTATES, NSTATES, A_data);
  model.B = slap_Matrix(NSTATES, NINPUTS, B_data);

  slap_Copy(X[0], model.x0);

  prob.ninputs = NINPUTS;
  prob.nstates = NSTATES;
  prob.nhorizon = NHORIZON;
  prob.Q = slap_Matrix(NSTATES, NSTATES, Q_data);
  slap_SetIdentity(prob.Q, 10e-1);
  prob.R = slap_Matrix(NINPUTS, NINPUTS, R_data);
  slap_SetIdentity(prob.R, 1e-1);
  prob.Qf = slap_Matrix(NSTATES, NSTATES, Qf_data);
  slap_SetIdentity(prob.Qf, 10e-1);
  prob.Xref = Xref;
  prob.Uref = Uref;
  prob.x0 = model.x0;
  prob.K = K;
  prob.d = d;
  prob.P = P;
  prob.p = p;

  float Q_temp_data[(NSTATES + NINPUTS) * (NSTATES + NINPUTS + 1)] = {0};
  Matrix Q_temp = slap_Matrix(NSTATES + NINPUTS, NSTATES + NINPUTS + 1,
                                       Q_temp_data);
  // Formulate problem as delta_x = x - x_hover = x (xhover = 0);
  // delta_u = u - u_hover
  // delta_x_ref = x_ref, delta_u_ref = 0

  // Compute and store A, B offline: dx = Adx + Bdu
  tiny_PQuadGetJacobians(&(model.A), &(model.B), xhover, uhover);

  tiny_BackwardPassLti(&prob, solver, model, &Q_temp);
  tiny_ForwardPassLti(X, U, prob, model);
  // for (int k = 0; k < NHORIZON - 1; ++k) {
  //   // Control input: du = - d - K*dx
  //   slap_Copy(U[k], prob.d[k]);              // u[k] = -d[k]
  //   slap_MatMulAdd(U[k], prob.Kinf, X[k], -1, -1);   // u[k] -= Kinf * dx[k]
  //   // Next state: dx = A*dx + B*du + f
  //   tiny_DynamicsLti(&X[k + 1], X[k], U[k], model);
  // }
  for (int k = 0; k < NHORIZON - 1; ++k) {
    // printf("ex[%d] = %.4f\n", k, slap_NormedDifference(X[k], Xref[k]));
    // PrintMatrix(X[k]);
  }

  for (int k = NHORIZON - 5; k < NHORIZON; ++k) {
    TEST(SumOfSquaredError(X[k].data, Xref[k].data, NSTATES) < 0.1);
  }
}

void AbsLqrLtiTest() {
  float X_data[NSTATES * NHORIZON] = {0};
  float U_data[NINPUTS * (NHORIZON - 1)] = {0};
  float K_data[NINPUTS * NSTATES * (NHORIZON - 1)] = {0};
  float d_data[NINPUTS * (NHORIZON - 1)] = {0};
  float P_data[NSTATES * NSTATES * (NHORIZON)] = {0};
  float p_data[NSTATES * NHORIZON] = {0};
  float A_data[NSTATES * NSTATES] = {0};
  float B_data[NSTATES * NINPUTS] = {0};
  float f_data[NSTATES] = {0};

  tiny_LtiModel model;
  tiny_InitLtiModel(&model);
  tiny_ProblemData prob;
  tiny_InitProblemData(&prob);
  tiny_AdmmSettings solver;
  tiny_InitSettings(&solver);

  Matrix X[NHORIZON];
  Matrix U[NHORIZON - 1];
  Matrix Xref[NHORIZON];
  Matrix Uref[NHORIZON - 1];
  Matrix K[NHORIZON - 1];
  Matrix d[NHORIZON - 1];
  Matrix P[NHORIZON];
  Matrix p[NHORIZON];

  float* Xptr = X_data;
  // float* Xref_ptr = Xref_data;
  float* Uptr = U_data;
  // float* Uref_ptr = Uref_data;
  float* Kptr = K_data;
  float* dptr = d_data;
  float* Pptr = P_data;
  float* pptr = p_data;

  Matrix xhover = slap_Matrix(NSTATES, 1, xhover_data);
  Matrix uhover = slap_Matrix(NINPUTS, 1, uhover_data);

  for (int i = 0; i < NHORIZON; ++i) {
    if (i < NHORIZON - 1) {
      U[i] = slap_Matrix(NINPUTS, 1, Uptr);
      // slap_SetConst(U[i], 0.01);
      Uptr += NINPUTS;
      Uref[i] = slap_Matrix(NINPUTS, 1, uhover_data);
      // Uref_ptr += NINPUTS;
      K[i] = slap_Matrix(NINPUTS, NSTATES, Kptr);
      Kptr += NINPUTS * NSTATES;
      d[i] = slap_Matrix(NINPUTS, 1, dptr);
      dptr += NINPUTS;
    }
    X[i] = slap_Matrix(NSTATES, 1, Xptr);
    Xptr += NSTATES;
    Xref[i] = slap_Matrix(NSTATES, 1, xg_data);
    // Xref_ptr += NSTATES;
    P[i] = slap_Matrix(NSTATES, NSTATES, Pptr);
    Pptr += NSTATES * NSTATES;
    p[i] = slap_Matrix(NSTATES, 1, pptr);
    pptr += NSTATES;
  }

  model.ninputs = NSTATES;
  model.nstates = NINPUTS;
  model.x0 = slap_Matrix(NSTATES, 1, x0_data);
  model.f = slap_Matrix(NSTATES, 1, f_data);
  model.A = slap_Matrix(NSTATES, NSTATES, A_data);
  model.B = slap_Matrix(NSTATES, NINPUTS, B_data);

  slap_Copy(X[0], model.x0);

  prob.ninputs = NINPUTS;
  prob.nstates = NSTATES;
  prob.nhorizon = NHORIZON;
  prob.Q = slap_Matrix(NSTATES, NSTATES, Q_data);
  slap_SetIdentity(prob.Q, 100e-1);
  prob.R = slap_Matrix(NINPUTS, NINPUTS, R_data);
  slap_SetIdentity(prob.R, 1e-1);
  prob.Qf = slap_Matrix(NSTATES, NSTATES, Qf_data);
  slap_SetIdentity(prob.Qf, 100e-1);
  prob.Xref = Xref;
  prob.Uref = Uref;
  prob.x0 = model.x0;
  prob.K = K;
  prob.d = d;
  prob.P = P;
  prob.p = p;

  float Q_temp_data[(NSTATES + NINPUTS) * (NSTATES + NINPUTS + 1)] = {0};
  Matrix Q_temp = slap_Matrix(NSTATES + NINPUTS, NSTATES + NINPUTS + 1,
                                       Q_temp_data);

  // Formulate problem as absolute: x, u
  // Compute and store A, B offline: x = Ax + Bu + f
  tiny_PQuadGetJacobians(&(model.A), &(model.B), xhover, uhover);
  tiny_PQuadNonlinearDynamics(&(model.f), xhover, uhover);
  slap_MatMulAdd(model.f, model.A, xhover, -1, 1);
  slap_MatMulAdd(model.f, model.B, uhover, -1, 1);

  tiny_BackwardPassLti(&prob, solver, model, &Q_temp);
  tiny_ForwardPassLti(X, U, prob, model);
  // for (int k = 0; k < NHORIZON - 1; ++k) {
  //   // Control input: u = - d - K*x
  //   slap_Copy(U[k], prob.d[k]);              // u[k] = -d[k]
  //   slap_MatMulAdd(U[k], prob.Kinf, X[k], -1, -1);   // u[k] -= Kinf * x[k]
  //   // Next state: x = A*x + B*u + f
  //   tiny_DynamicsLti(&X[k + 1], X[k], U[k], model);
  // }
  for (int k = 0; k < NHORIZON - 1; ++k) {
    printf("ex[%d] = %.4f\n", k, slap_NormedDifference(X[k], Xref[k]));
    // PrintMatrix(X[k]);
  }

  for (int k = NHORIZON - 5; k < NHORIZON; ++k) {
    TEST(SumOfSquaredError(X[k].data, Xref[k].data, NSTATES) < 0.1);
  }
}

int main() {
  printf("=== LQR LTV Test ===\n");
  // DeltaLqrLtvTest();
  // AbsLqrLtvTest();
  // DeltaLqrLtiTest();
  AbsLqrLtiTest();  // Most robust
  PrintTestResult();
  return TestResult();
}
