// Test all inequality-related functions and tiny_RiccatiConvergence

#include "simpletest.h"
#include "slap/slap.h"
#include "test_utils.h"
#include "tinympc/constraint_linear.h"
#include "tinympc/auxil.h"
#include "tinympc/model.h"
#include "tinympc/al_lqr.h"

#define NSTATES 4
#define NINPUTS 2
#define NHORIZON 3

// void IneqInputsTest() {
//   const float tol = 1e-6;
//   float u_max_data[NINPUTS] = {1, 1};
//   float u_min_data[NINPUTS] = {-1, -1};
//   float u_data[NINPUTS] = {1.1, 0.8};
//   float ans[NINPUTS * 2] = {
//       u_data[0] - u_max_data[0], u_data[1] - u_max_data[1],
//       -u_data[0] + u_min_data[0], -u_data[1] + u_min_data[1]};
//   float ineq_data[NINPUTS * 2];

//   tiny_ProblemData prob;
//   tiny_InitProblemData(&prob);

//   prob.nstates = NSTATES;
//   prob.ninputs = NINPUTS;
//   prob.nhorizon = NHORIZON;
//   prob.ncstr_states = 2 * NSTATES;
//   prob.ncstr_inputs = 2 * NINPUTS;
//   prob.ncstr_goal = NSTATES;
//   prob.u_max = slap_Matrix(NINPUTS, 1, u_max_data);
//   prob.u_min = slap_Matrix(NINPUTS, 1, u_min_data);
//   Matrix u = slap_Matrix(NINPUTS, 1, u_data);
//   Matrix mat = slap_Matrix(NINPUTS * 2, 1, ineq_data);
//   tiny_EvalInputConstraint(&mat, prob, u);
//   // slap_PrintMatrix(mat);
//   TEST(SumOfSquaredError(mat.data, ans, NINPUTS * 2) < tol);
// }

// void IneqInputsOffsetTest() {
//   const float tol = 1e-6;
//   float u_max_data[NINPUTS] = {1, 1};
//   float u_min_data[NINPUTS] = {-1, -1};
//   float ans[NINPUTS * 2] = {u_max_data[0], u_max_data[1], -u_min_data[0],
//                              -u_min_data[1]};
//   float ineq_data[NINPUTS * 2];

//   tiny_ProblemData prob;
//   tiny_InitProblemData(&prob);

//   prob.nstates = NSTATES;
//   prob.ninputs = NINPUTS;
//   prob.nhorizon = NHORIZON;
//   prob.ncstr_states = 2 * NSTATES;
//   prob.ncstr_inputs = 2 * NINPUTS;
//   prob.ncstr_goal = NSTATES;
//   prob.u_max = slap_Matrix(NINPUTS, 1, u_max_data);
//   prob.u_min = slap_Matrix(NINPUTS, 1, u_min_data);
//   Matrix mat = slap_Matrix(NINPUTS * 2, 1, ineq_data);
//   tiny_EvalInputConstraintOffset(&mat, prob);
//   // slap_PrintMatrix(mat);
//   TEST(SumOfSquaredError(mat.data, ans, NINPUTS * 2) < tol);
// }

// void IneqInputsJacobianTest() {
//   const float tol = 1e-6;
//   float u_max_data[NINPUTS] = {1, 1};
//   float u_min_data[NINPUTS] = {-1, -1};
//   float ans[NINPUTS * 2 * NINPUTS] = {1, 0, -1, 0, 0, 1, 0, -1};
//   float jac_data[NINPUTS * 2 * NINPUTS];

//   tiny_ProblemData prob;
//   tiny_InitProblemData(&prob);

//   prob.nstates = NSTATES;
//   prob.ninputs = NINPUTS;
//   prob.nhorizon = NHORIZON;
//   prob.ncstr_states = 2 * NSTATES;
//   prob.ncstr_inputs = 2 * NINPUTS;
//   prob.ncstr_goal = NSTATES;
//   prob.u_max = slap_Matrix(NINPUTS, 1, u_max_data);
//   prob.u_min = slap_Matrix(NINPUTS, 1, u_min_data);
//   Matrix mat = slap_Matrix(NINPUTS * 2, NINPUTS, jac_data);
//   tiny_EvalInputConstraintJacobian(&mat, prob);
//   // slap_PrintMatrix(mat);
//   TEST(SumOfSquaredError(mat.data, ans, NINPUTS * 2 * NINPUTS) < tol);
// }

void ActiveIneqMaskTest() {
  // const float tol = 1e-6;
  // float u_max_data[NINPUTS] = {1, 1};
  // float u_min_data[NINPUTS] = {-1, -1};
  float Acstr_input_data[2 * NINPUTS * NINPUTS] = {0};
  // [u_max, -u_min]
  float bcstr_input_data[2 * NINPUTS] = {1, 1, 1, 1};

  float u_data[NINPUTS] = {-2, 2};
  // float ans1[NINPUTS * 2] = {u_data[0] - u_max_data[0],    //-3
  //                             u_data[1] - u_max_data[1],    // 1
  //                             -u_data[0] + u_min_data[0],   // 1
  //                             -u_data[1] + u_min_data[1]};  // -3
  float dual_data[NINPUTS * 2] = {1, 0, 2, 0};
  // float ans2[NINPUTS * 2 * NINPUTS * 2] = {1, 0, 0, 0, 0, 1, 0, 0,
  //                                           0, 0, 1, 0, 0, 0, 0, 0};
  Matrix U[NHORIZON-1];

  tiny_Model model;
  tiny_InitModel(&model, NSTATES, NINPUTS, NHORIZON, 0, 0, 0.1);
  // tiny_InitModel(&model, NSTATES, NINPUTS, NHORIZON, 0, 1, 0.1);
  tiny_AdmmSettings stgs;
  tiny_InitSettings(&stgs);  //if switch on/off during run, initialize all
  tiny_AdmmData data;
  tiny_AdmmInfo info;
  tiny_AdmmSolution soln;
  tiny_AdmmWorkspace work;
  tiny_InitWorkspace(&work, &info, &model, &data, &soln, &stgs);
  
  float temp_data[work.data_size];
  T_INIT_ZEROS(temp_data);

  tiny_InitWorkspaceTempData(&work, temp_data);

  soln.U = U;
  U[0] = slap_Matrix(NINPUTS, 1, u_data);

  data.Acu = slap_Matrix(2 * NINPUTS, NINPUTS, Acstr_input_data);
  Matrix upper_half = slap_CreateSubMatrix(data.Acu, 0, 0, NINPUTS, NINPUTS);
  Matrix lower_half = slap_CreateSubMatrix(data.Acu, NINPUTS, 0, NINPUTS, NINPUTS);
  slap_SetIdentity(upper_half, 1);
  slap_SetIdentity(lower_half, -1);
  data.bcu = slap_Matrix(2 * NINPUTS, 1, bcstr_input_data);
  tiny_EvalInputConstraint(&work, 0);
  work.YU_hat = slap_Matrix(NINPUTS * 2, 1, dual_data);
  // Matrix mask = slap_Matrix(NINPUTS * 2, NINPUTS * 2, mask_data);
  tiny_ActiveIneqMask(&(work.cu_mask), work.YU_hat, work.cu);
  // PrintMatrix(work.cu_mask);
  // PrintMatrixT(work.cu);
  // TEST(SumOfSquaredError(work.cu.data, ans1, NINPUTS * 2) < tol);
  // TEST(SumOfSquaredError(work.cu_mask.data, ans2, NINPUTS * 2 * NINPUTS * 2) < tol);
}

void RiccatiConvergenceTest() {
  const float tol = 1e-6;
  tiny_Model model;
  tiny_InitModel(&model, NSTATES, 2, 3, 0, 0, 0.1);
  // tiny_InitModel(&model, NSTATES, NINPUTS, NHORIZON, 0, 1, 0.1);
  tiny_AdmmSettings stgs;
  tiny_InitSettings(&stgs);  //if switch on/off during run, initialize all
  tiny_AdmmData data;
  tiny_AdmmInfo info;
  tiny_AdmmSolution soln;
  tiny_AdmmWorkspace work;
  tiny_InitWorkspace(&work, &info, &model, &data, &soln, &stgs);
  
  float temp_data[work.data_size];
  T_INIT_ZEROS(temp_data);

  tiny_InitWorkspaceTempData(&work, temp_data);

  float d_data[4] = {1.2, -0.3, -2.1, 3.1};
  float ans = 3.744329045369811;
  Matrix d[2];
  float* dptr = d_data;
  for (int k = 0; k < 3 - 1; ++k) {
    d[k] = slap_Matrix(2, 1, dptr);
    dptr += 2;
  }
  soln.d = d;
  int res = tiny_CheckRiccati(&work);
  float norm_d_max = work.info->pri_res;
  TESTAPPROX(norm_d_max, ans, tol);
  TEST(res == 0);
}

void ProjectOrthantDualsTest() {
  float dual_data[2] = {-2, -1};
  float new_dual_data[2] = {10, -10};
  float ans[2] = {10, 0};
  Matrix dual = slap_Matrix(2, 1, dual_data);
  Matrix new_dual = slap_Matrix(2, 1, new_dual_data);
  tiny_ProjectOrthantDuals(&dual, new_dual);
  TEST(SumOfSquaredError(dual.data, ans, 2) < 1e-8);
}

int main() {
  printf("=== LQR Inequality Utility Test ===\n");
  // IneqInputsTest();
  // IneqInputsOffsetTest();
  // IneqInputsJacobianTest();
  ActiveIneqMaskTest();
  RiccatiConvergenceTest();
  ProjectOrthantDualsTest();
  PrintTestResult();
  return TestResult();
}
