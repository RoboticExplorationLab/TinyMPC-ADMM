#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "simpletest.h"
#include "slap/slap.h"
#include "test_utils.h"
#include "tinympc/cost_lqr.h"
#include "tinympc/model.h"
#include "tinympc/auxil.h"

#define NSTATES 4
#define NINPUTS 2
#define NHORIZON 2

float x_data[NSTATES] = {1.1, 1.2, 1.3, -4.3};
float u_data[NINPUTS] = {-2.1, 1.1};
float x_ref_data[NSTATES * NHORIZON] = {1.1, 1.2, 1.3, -4.2,
                                         1.2, 1.3, 1.3, -4.3};
float u_ref_data[NINPUTS * NHORIZON] = {-2.1, 1.4, 1};
float Q_data[NSTATES * NSTATES] = {0};  // NOLINT
float R_data[NINPUTS * NINPUTS] = {0};  // NOLINT
float q_data[NSTATES*NHORIZON] = {0};            // NOLINT
float r_data[NINPUTS*(NHORIZON-1)] = {0};            // NOLINT
float Pinf_data[NSTATES * NSTATES] = {0};
float ans_stage[2] = {0.04549999999999994, 0.1314999999999999};
float ans_term = 0.0049999999999999975;
float ans_gradx[NSTATES] = {-0.11, -0.12, -0.13, 0.42};
float ans_gradu[NINPUTS] = {0.21, -0.14};
float ans_gradxf[NSTATES] = {-0.6, -0.65, -0.65, 2.15};

void AddCostTest() {
  const float tol = 1e-6;
  Matrix Uref[NHORIZON];
  Matrix Xref[NHORIZON];
  Matrix U[NHORIZON-1];
  Matrix X[NHORIZON];
  float* uptr = u_ref_data;
  float* xptr = x_ref_data;
  for (int i = 0; i < NHORIZON; ++i) {
    Uref[i] = slap_Matrix(NINPUTS, 1, uptr);
    uptr += NINPUTS;
    Xref[i] = slap_Matrix(NSTATES, 1, xptr);
    xptr += NSTATES;
  }
  Matrix x = slap_Matrix(NSTATES, 1, x_data);
  Matrix u = slap_Matrix(NINPUTS, 1, u_data);

  tiny_Model model;
  tiny_InitModel(&model, NSTATES, NINPUTS, NHORIZON, 0, 0, 0.1);
  tiny_AdmmSettings stgs;
  tiny_InitSettings(&stgs);  //if switch on/off during run, initialize all
  tiny_AdmmData data;
  tiny_AdmmInfo info;
  tiny_AdmmSolution soln;
  tiny_AdmmWorkspace work;
  tiny_InitWorkspace(&work, &info, &model, &data, &soln, &stgs);
  
  float temp_data[work.data_size];
  T_INIT_ZEROS(temp_data);

  tiny_InitWorkspaceTempData(&work, 0, 0, 0, 0, temp_data);

  data.Q = slap_Matrix(NSTATES, NSTATES, Q_data);
  slap_SetIdentity(data.Q, 0.1);
  data.R = slap_Matrix(NINPUTS, NINPUTS, R_data);
  slap_SetIdentity(data.R, 1);
  soln.Pinf = slap_Matrix(NSTATES, NSTATES, Pinf_data);
  slap_SetIdentity(soln.Pinf, 0.5);
  data.Xref = Xref;
  data.Uref = Uref;

  soln.X = X;
  soln.U = U;

  work.info->obj_val = 0.0;
  for (int k = 0; k < NHORIZON - 1; ++k) {
    soln.X[k] = x;
    soln.U[k] = u;
    tiny_AddStageCost(&work, k);
    TESTAPPROX(work.info->obj_val, ans_stage[k], tol);
  }
  soln.X[NHORIZON - 1] = x;
  work.info->obj_val = 0.0;
  tiny_AddTerminalCost(&work);
  TESTAPPROX(work.info->obj_val, ans_term, tol);
}

int main() {
  printf("=== Cost Test ===\n");
  AddCostTest();
  PrintTestResult();
  return TestResult();
}
