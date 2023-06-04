#include "constraint_linear.h"

#ifdef __cplusplus
extern "C" {
#endif

enum tiny_ErrorCode tiny_SetInputBound(tiny_AdmmWorkspace* work, float* Ac_data, float* lc_data, float* uc_data) {
  int n = work->data->model->ninputs;
  work->stgs->en_cstr_inputs = EN_CSTR_INPUTS;
  work->data->Acu = Eigen::Map<Eigen::Matrix4f>(Ac_data);
  (work->data->Acu).setIdentity();
  work->data->lcu = Eigen::Map<Eigen::Vector4f>(lc_data);
  work->data->ucu = Eigen::Map<Eigen::Vector4f>(uc_data);
  return TINY_NO_ERROR;
}

enum tiny_ErrorCode tiny_SetStateBound(tiny_AdmmWorkspace* work, float* Ac_data, float* lc_data, float* uc_data) {
  int n = work->data->model[0].nstates;
  work->stgs->en_cstr_states = EN_CSTR_STATES;
  work->data->Acx = Eigen::Map<Eigen::Matrix12f>(Ac_data);
  (work->data->Acx).setIdentity();
  work->data->lcx = Eigen::Map<Eigen::Vector12f>(lc_data);
  work->data->ucx = Eigen::Map<Eigen::Vector12f>(uc_data);
  return TINY_NO_ERROR;
}

// enum tiny_ErrorCode tiny_ProjectInput(tiny_AdmmWorkspace* work) {
//   int n = work->data->model[0].ninputs;
//   int N = work->data->model[0].ninputs;

//   for (int k = 0; k < N - 1; ++k) {
//     for (int i = 0; i < n; ++i) {

//       work->ZU_new[k].data[i] = T_MIN(T_MAX(z[i],
//                                 work->data->lcu[i]),  // Between lower
//                                 work->data->ucu[i]);  // and upper bounds
//     } 
//   }
//   return TINY_NO_ERROR;
// }

int IsConstrained(tiny_AdmmWorkspace* work) {
  if (!work->stgs->en_cstr_goal && 
      !work->stgs->en_cstr_inputs && 
      !work->stgs->en_cstr_states) {
    return 0; // unconstrained
  }
  return 1;    
}

#ifdef __cplusplus
}
#endif