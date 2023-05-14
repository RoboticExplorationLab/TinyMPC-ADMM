#ifndef TYPES_H
# define TYPES_H

# ifdef __cplusplus
extern "C" {
# endif // ifdef __cplusplus


#include "slap/slap.h"
#include "constants.h"
#include "errors.h"

// for a horizon of N x(0)->x(N-1), need N-1 matrices
typedef struct {
  int nstates;
  int ninputs;
  int nhorizon;

  int ltv;            ///< Boolean, true if model is LTV  
  int affine;         ///< Boolean, true if model is affine
  sfloat dt;          ///< Sample time Ts of the discrete model

  Matrix* A;
  Matrix* B;
  Matrix* f;

  void (*get_jacobians)(Matrix*, Matrix*, const Matrix, const Matrix);
  void (*get_nonl_model)(Matrix*, const Matrix, const Matrix);
  int data_size;
} tiny_Model;

/**
 * Solution structure
 */
typedef struct {
  Matrix* X;      ///< State trajectory solution 
  Matrix* U;      ///< Input trajectory solution

  Matrix Kinf;    ///< Feedback gain of IHLQR
  Matrix* d;      ///< Feedforward gain
  Matrix Pinf;    ///< Terminal cost Hessian of IHLQR
  Matrix* p;      ///< Terminal cost gradient
  
  Matrix* YU;     ///< Dual variables for input constraints
  Matrix* YX;     ///< Dual variables for state constraints
  Matrix YG;      ///< Dual variables for goal constraint

  int data_size;
} tiny_ADMMSolution;


/**
 * Solver return information
 */
typedef struct {
  int iter;           ///< Number of AL iterations taken
  int iter_riccati;   ///< Number of Riccati iterations taken
  int status_val;     ///< Integer, status defined in constants.h

  sfloat obj_val;     ///< primal objective
  sfloat pri_res;     ///< norm of primal residual
  sfloat dua_res;     ///< norm of dual residual
} tiny_ADMMInfo;


/**********************************
* Main structures and Data Types *
**********************************/

/**
 * Settings struct
 */
typedef struct {
  sfloat reg_min;             ///< Minimum regularization
  sfloat reg_max;             ///< Maximum regularization
  sfloat reg_mul;             ///< Regularization update multiplier
  int    en_reg_update;       ///< Boolean, enable regularization update (tighter solve)
  
  sfloat rho_init;            ///< Initial rho
  sfloat rho_max;             ///< Maximum rho
  sfloat rho_mul;             ///< Penalty multiplier

  sfloat alpha_mul;           ///< Line-search step multiplier

  int    max_iter;            ///< Maximum number of AL iterations
  int    max_iter_riccati;    ///< Maximum number of Riccati solve iterations
  int    max_iter_ls;         ///< Maximum number of line-search iterations

  sfloat tol_abs_prim;        ///< Riccati solve tolerance
  sfloat tol_abs_dual;        ///< Constraint tolerance

  int    en_cstr_states;      ///< Boolean, enable inequality constraints on states
  int    en_cstr_inputs;      ///< Boolean, enable inequality constraints on inputs
  int    en_cstr_goal;        ///< Boolean, enable equality constraint on goal

  int    verbose;             ///< Integer, level to write out progress
  int    adaptive_horizon;    ///< Integer, after `adaptive_horizon` steps, use the second model with longer interval; if 0, disabled 
  int    check_riccati;       ///< Boolean, if 0, then termination checking is disabled
  int    check_termination;   ///< Integer, check termination interval; if 0, then termination checking is disabled
  int    warm_start;          ///< boolean, enable warm start
  sfloat time_limit;          ///< Time limit of each MPC step; if 0, disabled
} tiny_ADMMSettings;

// void tiny_InitSettings(tiny_ADMMSettings* solver);


/**
 * Data structure
 */
typedef struct {
  tiny_Model* model;    ///< System model
  Matrix x0;

  Matrix  Q;
  Matrix  R;
  Matrix* q;
  Matrix* r;
  
  Matrix* X_ref;
  Matrix* U_ref;

  Matrix Acx;
  Matrix ucx;
  Matrix lcx;
  Matrix Acu;
  Matrix ucu;
  Matrix lcu;
  
  int data_size;
} tiny_ADMMData;

// void tiny_InitProblemData(tiny_ProblemData* prob);

typedef struct {
  tiny_ADMMData*        data;      ///< problem data
  tiny_ADMMSettings*    stgs;      ///< problem settings
  tiny_ADMMSolution*    soln;      ///< problem solution
  tiny_ADMMInfo*        info;      ///< solver information

  sfloat reg;
  sfloat alpha;
  sfloat rho;

  // Temporary data
  Matrix Qu;          ///< temporary 
  Matrix Quu_inv;     ///< mxm cache for (R + B'*Pinf*B)\I 
  Matrix AmBKt;       ///< nxn cache for (A - BKinf)'
  Matrix coeff_d2p;   ///< nxm cache for Kinf'*R - AmBKt*Pinf*B
  
  Matrix* ZU;         ///< Slack variable for input
  Matrix* ZU_new;     ///< Updated slack variable for input
  Matrix* ZX;         ///< Slack variable for input
  Matrix* ZX_new;     ///< Updated slack variable for input

  int data_size;      ///< sum data size of all temporary data //TODO: + model + solution 
  int first_run;      ///< flag indicating whether the solve function has been run before
} tiny_ADMMWorkspace;


# ifdef __cplusplus
}
# endif // ifdef __cplusplus

#endif // ifndef TYPES_H
