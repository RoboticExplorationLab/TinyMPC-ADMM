#include "simpletest.h"
// #include "test_utils.h"
#include "tinympc/model.h"

# ifdef __cplusplus
extern "C" {
# endif // ifdef __cplusplus

void tiny_InitModel_Test() {
  // Test good 
  int NSTATES = 2;
  int NINPUTS = 1;
  int NHORIZON = 3;
  tiny_Model ltv_model;
  int true_size = NSTATES*(NSTATES + NINPUTS + 1)*(NHORIZON-1);
  TEST(tiny_InitModel(&ltv_model, NSTATES, NINPUTS, NHORIZON, 1, 1, 0.1) == TINY_NO_ERROR);

  tiny_Model lti_model;
  true_size = NSTATES*(NSTATES + NINPUTS + 1);
  TEST(tiny_InitModel(&lti_model, NSTATES, NINPUTS, NHORIZON, 0, 1, 0.1) == TINY_NO_ERROR);

  tiny_Model ltv_model2;
  true_size = NSTATES*(NSTATES + NINPUTS)*(NHORIZON-1);
  TEST(tiny_InitModel(&ltv_model2, NSTATES, NINPUTS, NHORIZON, 1, 0, 0.1) == TINY_NO_ERROR);
  TEST(ltv_model2.dt == (float)(0.1) && ltv_model2.affine == 0);

  tiny_Model lti_model2;
  true_size = NSTATES*(NSTATES + NINPUTS);
  TEST(tiny_InitModel(&lti_model2, NSTATES, NINPUTS, NHORIZON, 0, 0, 0.1) == TINY_NO_ERROR);     

  // Test bad (there are more bad cases)
  // NSTATES = 1;
  // NINPUTS = 0;
  // int NHORIZON = 0;
  // TEST(tiny_InitModel(&model, NSTATES, NINPUTS, NHORIZON, 1, 0.1) == TINY_SLAP_ERROR);
}

void tiny_InitModel_Test() {
  const float tol = 1e-6;
  const int NSTATES = 2;
  const int NINPUTS = 2;
  const int NHORIZON = 3;
  float A_data[] = {1, 2, 3, 4, 5, 6, 7, 8};  
  float B_data[] = {1.1, 2.2, 1.3, 2.3};        
  float f_data[] = {2.1, 1.2, 3.3, 1.3};     
  Eigen::MatrixNf A[NHORIZON-1];
  Eigen::MatrixNMf B[NHORIZON-1];
  Eigen::VectorNf f[NHORIZON-1];  
  tiny_Model ltv_model;
  tiny_InitModel(&ltv_model, NSTATES, NINPUTS, NHORIZON, 1, 1, 0.1);
  TEST(ltv_model.nhorizon == NHORIZON);
  tiny_InitModel(&ltv_model, A, B, f, A_data, B_data, f_data);
  for (int k = 0; k < NHORIZON-1; ++k) {
    TEST(ltv_model.A[k].rows() == ltv_model.nstates);
    TEST(ltv_model.A[k].cols() == ltv_model.nstates);

    TEST(ltv_model.B[k].rows() == ltv_model.nstates);
    TEST(ltv_model.B[k].cols() == ltv_model.ninputs);

    TEST(ltv_model.f[k].rows() == ltv_model.nstates);
    TEST(ltv_model.f[k].cols() == 1);
  }
  printf("%d\n", ltv_model.A[0].rows());
  printf("%d\n", ltv_model.nstates);
  // TEST(SumOfSquaredErrorMatrices(A_data, ltv_model.A, NHORIZON-1) < tol);
  // TEST(SumOfSquaredErrorMatrices(B_data, ltv_model.B, NHORIZON-1) < tol);
  // TEST(SumOfSquaredErrorMatrices(f_data, ltv_model.f, NHORIZON-1) < tol);

  tiny_Model ltv_model2;
  tiny_InitModel(&ltv_model2, NSTATES, NINPUTS, NHORIZON, 1, 0, 0.1);
  tiny_InitModel(&ltv_model2, A, B, TINY_NULL, A_data, B_data, TINY_NULL);
  for (int k = 0; k < NHORIZON-1; ++k) {
    TEST(ltv_model2.A[k].rows() == ltv_model2.nstates);
    TEST(ltv_model2.A[k].cols() == ltv_model2.nstates);

    TEST(ltv_model2.B[k].rows() == ltv_model2.nstates);
    TEST(ltv_model2.B[k].cols() == ltv_model2.ninputs);
  }
  // TEST(SumOfSquaredErrorMatrices(A_data, ltv_model2.A, NHORIZON-1) < tol);
  // TEST(SumOfSquaredErrorMatrices(B_data, ltv_model2.B, NHORIZON-1) < tol);
  // TEST(ltv_model2.f == TINY_NULL);

  float A2_data[] = {1, 2, 3, 4};  
  float B2_data[] = {1.1, 2.2};        
  float f2_data[] = {2.1, 1.2};     
  Eigen::MatrixNf A2;
  Eigen::MatrixNMf B2;     
  Eigen::VectorNf f2;      
  tiny_Model lti_model;
  tiny_InitModel(&lti_model, NSTATES, NINPUTS, NHORIZON, 0, 1, 0.1);
  TEST(lti_model.nhorizon == NHORIZON);
  tiny_InitModel(&lti_model, &A2, &B2, &f2, A2_data, B2_data, f2_data);

  TEST(lti_model.A[0].rows() == lti_model.nstates);
  TEST(lti_model.A[0].cols() == lti_model.nstates);

  TEST(lti_model.B[0].rows() == lti_model.nstates);
  TEST(lti_model.B[0].cols() == lti_model.ninputs);

  TEST(lti_model.f[0].rows() == lti_model.nstates);
  TEST(lti_model.f[0].cols() == 1);

  // TEST(SumOfSquaredErrorMatrices(A2_data, lti_model.A, 1) < tol);
  // TEST(SumOfSquaredErrorMatrices(B2_data, lti_model.B, 1) < tol);
  // TEST(SumOfSquaredErrorMatrices(f2_data, lti_model.f, 1) < tol);

  tiny_Model lti_model2;
  tiny_InitModel(&lti_model2, NSTATES, NINPUTS, NHORIZON, 0, 0, 0.1);
  TEST(lti_model2.nhorizon == NHORIZON);
  tiny_InitModel(&lti_model2, &A2, &B2, &f2, A2_data, B2_data, f2_data);

  TEST(lti_model2.A[0].rows() == lti_model.nstates);
  TEST(lti_model2.A[0].cols() == lti_model.nstates);

  TEST(lti_model2.B[0].rows() == lti_model.nstates);
  TEST(lti_model2.B[0].cols() == lti_model.ninputs);

  // TEST(SumOfSquaredErrorMatrices(A2_data, lti_model2.A, 1) < tol);
  // TEST(SumOfSquaredErrorMatrices(B2_data, lti_model2.B, 1) < tol);
  // TEST(ltv_model2.f == TINY_NULL);
}

// void tiny_InitModelDataMatrix_Test() {
//   const float tol = 1e-8;
//   const int NSTATES = 2;
//   const int NINPUTS = 2;
//   const int NHORIZON = 3;
//   float A_data[] = {1, 2, 3, 4, 5, 6, 7, 8};  
//   float B_data[] = {1.1, 2.2, 1.3, 2.3};        
//   float f_data[] = {2.1, 1.2, 3.3, 1.3};     
//   Matrix A[NHORIZON-1];
//   Matrix B[NHORIZON-1];     
//   Matrix f[NHORIZON-1];        
//   tiny_Model model;
//   tiny_InitModel(&model, NSTATES, NINPUTS, NHORIZON, 1, 0.1);

//   float* A_ptr = A_data;
//   float* B_ptr = B_data;
//   float* f_ptr = f_data;

//   for (int k = 0; k < model.nhorizon-1; ++k) {
//     A[k] = slap_Matrix(model.nstates, model.nstates, A_ptr);
//     A_ptr += model.nstates * model.nstates;
//     B[k] = slap_Matrix(model.nstates, model.ninputs, B_ptr);
//     B_ptr += model.nstates * model.ninputs;
//     f[k] = slap_Matrix(model.nstates, 1, f_ptr);
//     f_ptr += model.nstates;      
//   }  

//   tiny_InitModelDataMatrix(&model, A, B, f);

//   for (int k = 0; k < NHORIZON-1; ++k) {
//     TEST(model.A[k].rows() == model.nstates);
//     TEST(model.A[k].cols() == model.nstates);

//     TEST(model.B[k].rows() == model.nstates);
//     TEST(model.B[k].cols() == model.ninputs);

//     TEST(model.f[k].rows() == model.nstates);
//     TEST(model.f[k].cols() == 1);

//   }
//   TEST(SumOfSquaredErrorMatrices(A_data, model.A, NHORIZON-1) < tol);
//   TEST(SumOfSquaredErrorMatrices(B_data, model.B, NHORIZON-1) < tol);
//   TEST(SumOfSquaredErrorMatrices(f_data, model.f, NHORIZON-1) < tol);

//   tiny_Model model_not_affine;
//   tiny_InitModel(&model_not_affine, NSTATES, NINPUTS, NHORIZON, 0, 0.1);
//   tiny_InitModelDataMatrix(&model_not_affine, A, B, f);
//   TEST(SumOfSquaredErrorMatrices(A_data, model_not_affine.A, NHORIZON-1) < tol);
//   TEST(SumOfSquaredErrorMatrices(B_data, model_not_affine.B, NHORIZON-1) < tol);
//   TEST(model_not_affine.f == TINY_NULL);
// }

// void tiny_InitModelMemory_Test() {
//   const float tol = 1e-8;
//   const int NSTATES = 2;
//   const int NINPUTS = 2;
//   const int NHORIZON = 3;
//   float A_data[] = {1, 2, 3, 4, 5, 6, 7, 8};  
//   float B_data[] = {1.1, 2.2, 1.3, 2.3};        
//   float f_data[] = {2.1, 1.2, 3.3, 1.3};  
                    
//   tiny_Model model;
//   tiny_InitModel(&model, NSTATES, NINPUTS, NHORIZON, 1, 0.1);
//   float data[model.data_size];  // data
//   Matrix mats[(NHORIZON-1)*3];  // array of matrices A, B, f
//   tiny_InitModelMemory(&model, mats, data);

//   // just for testing, use UpdateJacobian isFillModeltead
//   tiny_FillModelMemory(&model, A_data, B_data, f_data);  

//   for (int k = 0; k < NHORIZON-1; ++k) {
//     TEST(model.A[k].rows() == model.nstates);
//     TEST(model.A[k].cols() == model.nstates);

//     TEST(model.B[k].rows() == model.nstates);
//     TEST(model.B[k].cols() == model.ninputs);

//     TEST(model.f[k].rows() == model.nstates);
//     TEST(model.f[k].cols() == 1);

//   }
//   TEST(SumOfSquaredErrorMatrices(A_data, model.A, NHORIZON-1) < tol);
//   TEST(SumOfSquaredErrorMatrices(B_data, model.B, NHORIZON-1) < tol);
//   TEST(SumOfSquaredErrorMatrices(f_data, model.f, NHORIZON-1) < tol);

//   tiny_Model model_not_affine;
//   tiny_InitModel(&model_not_affine, NSTATES, NINPUTS, NHORIZON, 0, 0.1);
//   float data2[model_not_affine.data_size];  // data
//   Matrix mats2[(NHORIZON-1)*2];  // array of matrices A and B
//   tiny_InitModelMemory(&model_not_affine, mats2, data2);
//   tiny_FillModelMemory(&model_not_affine, A_data, B_data, f_data);  
//   for (int k = 0; k < NHORIZON-1; ++k) {
//     TEST(model_not_affine.A[k].rows() == model.nstates);
//     TEST(model_not_affine.A[k].cols() == model.nstates);

//     TEST(model_not_affine.B[k].rows() == model.nstates);
//     TEST(model_not_affine.B[k].cols() == model.ninputs);
//   }
//   TEST(SumOfSquaredErrorMatrices(A_data, model_not_affine.A, NHORIZON-1) < tol);
//   TEST(SumOfSquaredErrorMatrices(B_data, model_not_affine.B, NHORIZON-1) < tol);
//   TEST(model_not_affine.f == TINY_NULL);
// }

// void tiny_SetModelJacFunc_Test() {
//   int NSTATES = 2;
//   int NINPUTS = 1;
//   int NHORIZON = 3;
//   tiny_Model model;
//   tiny_InitModel(&model, NSTATES, NINPUTS, NHORIZON, 1, 1, 0.1);
//   tiny_SetModelJacFunc(&model, tiny_Bicycle5dGetJacobians);
//   TEST(model.get_jacobians == tiny_Bicycle5dGetJacobians);
// }

// void tiny_SetModelNonlFunc_Test() {
//   int NSTATES = 2;
//   int NINPUTS = 1;
//   int NHORIZON = 3;
//   tiny_Model model;
//   tiny_InitModel(&model, NSTATES, NINPUTS, NHORIZON, 0, 1, 0.1);
//   tiny_SetModelNonlFunc(&model, tiny_Bicycle5dNonlinearDynamics);
//   TEST(model.get_nonl_model == tiny_Bicycle5dNonlinearDynamics);
// }

int main() {
  tiny_InitModel_Test();
  tiny_InitModel_Test();
  // tiny_InitModelDataMatrix_Test();
  // tiny_InitModelMemory_Test();
  // tiny_SetModelJacFunc_Test();
  // tiny_SetModelNonlFunc_Test();
  PrintTestResult();
  return TestResult();
}

# ifdef __cplusplus
}
# endif // ifdef __cplusplus