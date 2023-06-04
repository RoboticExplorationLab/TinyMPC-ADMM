#ifndef TYPES_H
# define TYPES_H

# include <Eigen.h>
#include "constants.h"
#include "errors.h"

# ifdef __cplusplus
extern "C" {
# endif // ifdef __cplusplus

/**
 * Linear algebra type
 */
namespace Eigen
{ 
    typedef Matrix<float, 12, 12> Matrix12f;
    typedef Matrix<float, 12, 4>  Matrix12x4f;
    typedef Matrix<float, 4, 12>  Matrix4x12f;
    typedef Vector<float, 12>     Vector12f;  
}

// for a horizon of N x(0)->x(N-1), need N-1 matrices
typedef struct {
  int nstates;
  int ninputs;
  int nhorizon;

  int ltv;            ///< Boolean, true if model is LTV  
  int affine;         ///< Boolean, true if model is affine
  float dt;          ///< Sample time Ts of the discrete model

  Eigen::Matrix12f*     A;
  Eigen::Matrix12x4f*   B;
  Eigen::Vector12f*     f;

  void (*get_jacobians)(Eigen::Matrix12f*, Eigen::Matrix12x4f*, Eigen::Vector12f*, Eigen::Vector4f*);
  void (*get_nonl_model)(Eigen::Vector12f*, Eigen::Vector12f*, Eigen::Vector4f*);
  int data_size;
} tiny_Model;


/**
 * Solution structure
 */
typedef struct {
  Eigen::Vector12f* X;      ///< State trajectory solution 
  Eigen::Vector4f* U;      ///< Input trajectory solution

  Eigen::Matrix4x12f Kinf;    ///< Feedback gain of IHLQR
  Eigen::Vector4f* d;      ///< Feedforward gain
  Eigen::Matrix12f Pinf;    ///< Terminal cost Hessian of IHLQR
  Eigen::Vector12f* p;      ///< Terminal cost gradient
  
  Eigen::Vector4f* YU;     ///< Dual variables for input constraints
  Eigen::Vector12f* YX;     ///< Dual variables for state constraints
  Eigen::Vector12f YG;      ///< Dual variables for goal constraint

  int data_size;
} tiny_AdmmSolution;


/**
 * Solver return information
 */
typedef struct {
  int iter;           ///< Number of AL iterations taken
  int iter_riccati;   ///< Number of Riccati iterations taken
  int status_val;     ///< Integer, status defined in constants.h

  float obj_val;     ///< primal objective
  float pri_res;     ///< norm of primal residual
  float dua_res;     ///< norm of dual residual
} tiny_AdmmInfo;


/**********************************
* Main structures and Data Types *
**********************************/

/**
 * Settings struct
 */
typedef struct {
  float reg_min;             ///< Minimum regularization
  float reg_max;             ///< Maximum regularization
  float reg_mul;             ///< Regularization update multiplier
  int    en_reg_update;       ///< Boolean, enable regularization update (tighter solve)
  
  float rho_init;            ///< Initial rho
  float rho_max;             ///< Maximum rho
  float rho_mul;             ///< Penalty multiplier

  float alpha_mul;           ///< Line-search step multiplier

  int    max_iter;            ///< Maximum number of AL iterations
  int    max_iter_riccati;    ///< Maximum number of Riccati solve iterations
  int    max_iter_ls;         ///< Maximum number of line-search iterations

  float tol_abs_prim;        ///< Riccati solve tolerance
  float tol_abs_dual;        ///< Constraint tolerance

  int    en_cstr_states;      ///< Boolean, enable inequality constraints on states
  int    en_cstr_inputs;      ///< Boolean, enable inequality constraints on inputs
  int    en_cstr_goal;        ///< Boolean, enable equality constraint on goal

  int    verbose;             ///< Integer, level to write out progress
  int    adaptive_horizon;    ///< Integer, after `adaptive_horizon` steps, use the second model with longer interval; if 0, disabled 
  int    check_riccati;       ///< Boolean, if 0, then termination checking is disabled
  int    check_termination;   ///< Integer, check termination interval; if 0, then termination checking is disabled
  int    warm_start;          ///< boolean, enable warm start
  float time_limit;          ///< Time limit of each MPC step; if 0, disabled
} tiny_AdmmSettings;

// void tiny_InitSettings(tiny_AdmmSettings* solver);


/**
 * Data structure
 */
typedef struct {
  tiny_Model* model;    ///< System model
  Eigen::Vector12f  x0;

  Eigen::Matrix12f  Q;
  Eigen::Matrix4f   R;
  Eigen::Vector12f* q;
  Eigen::Vector4f*  r;
  Eigen::Vector4f*  r_tilde;
  
  Eigen::Vector12f* Xref;
  Eigen::Vector4f*  Uref;

  Eigen::Matrix12f  Acx;
  Eigen::Vector12f  ucx;
  Eigen::Vector12f  lcx;
  Eigen::Matrix4f   Acu;
  Eigen::Vector4f   ucu;
  Eigen::Vector4f   lcu;
  
  int data_size;
} tiny_AdmmData;

// void tiny_InitProblemData(tiny_ProblemData* prob);

typedef struct {
  tiny_AdmmData*        data;      ///< problem data
  tiny_AdmmSettings*    stgs;      ///< problem settings
  tiny_AdmmSolution*    soln;      ///< problem solution
  tiny_AdmmInfo*        info;      ///< solver information

  float reg;
  float alpha;
  float rho;

  // Temporary data
  Eigen::Vector4f     Qu;          ///< temporary 
  Eigen::Matrix4f     Quu_inv;     ///< mxm cache for (R + B'*Pinf*B)\I 
  Eigen::Matrix12f    AmBKt;       ///< nxn cache for (A - BKinf)'
  Eigen::Matrix12x4f  coeff_d2p;   ///< nxm cache for Kinf'*R - AmBKt*Pinf*B
  
  Eigen::Vector4f*  ZU;         ///< Slack variable for input
  Eigen::Vector4f*  ZU_new;     ///< Updated slack variable for input
  Eigen::Vector12f* ZX;         ///< Slack variable for input
  Eigen::Vector12f* ZX_new;     ///< Updated slack variable for input

  int data_size;      ///< sum data size of all temporary data //TODO: + model + solution 
  int first_run;      ///< flag indicating whether the solve function has been run before
} tiny_AdmmWorkspace;


# ifdef __cplusplus
}
# endif // ifdef __cplusplus

#endif // ifndef TYPES_H
