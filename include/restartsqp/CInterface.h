/* Copyright (C) 2020
 * All Rights Reserved.
 * This code is published under the Eclipse Public License.
 *
 * This C interface is derived from the file IpCInterface.h from the Ipopt project
 */

#ifndef __RESTARTSQP_CINTERFACE_H__
#define __RESTARTSQP_CINTERFACE_H__

#ifdef __cplusplus
extern "C"
{
#endif

/* This includes the SqpSolverExitStatus enum type */
#include "restartsqp/ReturnCodes.h"

/** Structure collecting all information about the problem definition and solve statistics etc. */
struct SqpProblemInfo;

/** Pointer to an RestartSqp Problem. */
typedef struct SqpProblemInfo* SqpProblem;

/** define a boolean type for C */
typedef char Bool;
#ifndef TRUE
# define TRUE (1)
#endif
#ifndef FALSE
# define FALSE (0)
#endif

/** A pointer for anything that is to be passed between the called and individual callback function */
typedef void* UserDataPtr;

/** Type defining the callback function for evaluating the value of the objective function.
 *
 *  Return value should be set to false if there was a problem doing the evaluation.
 *
 *  See also CrossoverSqpTNlp::eval_objective_value.
 */
typedef Bool (*Eval_Obj_CB)(
   int         num_variables,
   double*     primal_variables,
   Bool        new_primal_variables,
   double*     objective_value,
   UserDataPtr user_data
);

/** Type defining the callback function for evaluating the gradient of the objective function.
 *
 *  Return value should be set to false if there was a problem doing the evaluation.
 *
 *  See also CrossoverSqpTNlp::eval_objective_gradient.
 */
typedef Bool (*Eval_Obj_Grad_CB)(
   int         num_variables,
   double*     primal_variables,
   Bool        new_priml_variables,
   double*     objective_gradient,
   UserDataPtr user_data
);

/** Type defining the callback function for evaluating the value of the constraint functions.
 *
 *  Return value should be set to false if there was a problem doing the evaluation.
 *
 *  See also CrossoverSqpTNlp::eval_constraint_values.
 */
typedef Bool (*Eval_Constr_CB)(
   int         num_variables,
   double*     primal_variables,
   Bool        new_primal_variables,
   int         num_constraints,
   double*     constraint_values,
   UserDataPtr user_data
);

/** Type defining the callback function for evaluating the Jacobian of the constrant functions.
 *  Counting of row and column indices starts at 0.
 *
 *  Return value should be set to false if there was a problem doing the evaluation.
 *
 *  See also CrossoverSqpTNlp::eval_jac_g.
 */
typedef Bool (*Eval_Constr_Jac_CB)(
   int         num_variables,
   double*     primal_variables,
   Bool        new_primal_variables,
   int         num_constraints,
   int         num_nonzeros_jacobian,
   int*        row_indices,
   int*        column_indices,
   double*     nonzero_values,
   UserDataPtr user_data
);

/** Type defining the callback function for evaluating the Hessian of the Lagrangian function.
 *
 *  Return value should be set to false if there was a problem doing the evaluation.
 *
 *  IMPORTANT:The sign of the multipliers is reversed compared to Ipopt!!!
 *
 *  See also CrossoverSqpTNlp::eval_h.
 */
typedef Bool (*Eval_Lag_Hess_CB)(
   int         num_variables,
   double*     primal_variables,
   Bool        new_primal_variables,
   double      objective_scaling_factor,
   int         num_constraints,
   double*     constraint_multipliers,
   Bool        new_constraint_multipliers,
   int         num_nonzeros_hessian,
   int*        row_indices,
   int*        column_indices,
   double*     nonzero_values,
   UserDataPtr user_data
);

/** Function for creating a new SQP Problem object.
 *
 *  This function returns an object that can be passed to the CrossoverSqpSolve call.
 *  It contains the basic definition of the optimization problem, such
 *  as number of variables and constraints, bounds on variables and
 *  constraints, information about the derivatives, and the callback
 *  function for the computation of the optimization problem
 *  functions and derivatives.  During this call, the options file
 *  sqp.opt is read as well.
 *
 *  If NULL is returned, there was a problem with one of the inputs
 *  or reading the options file.
 *
 *  See also CrossoverSqpTNlp::get_nlp_info and CrossoverSqpTNlp::get_bounds_info.
 */
SqpProblem CreateSqpProblem(
   int         num_variables,           /**< Number of optimization variables */
   double*       variable_lower_bounds,         /**< Lower bounds on variables
                               *
                               * This array of size n is copied internally, so that the
                               * caller can change the incoming data after
                               * return without that SqpProblem is modified.
                               */
   double*       variable_upper_bounds,         /**< Upper bounds on variables
                               *
                               * This array of size n is copied internally,
                               * so that the caller can change the incoming data after
                               * return without that SqpProblem is modified.
                               */
   int          num_constraints,           /**< Number of constraints */
   double*        constraint_lower_bounds,         /**< Lower bounds on constraints
                               *
                               * This array of size m is copied internally,
                               * so that the caller can change the incoming data after
                               * return without that SqpProblem is modified.
                               */
   double*        constraint_upper_bounds,         /**< Upper bounds on constraints
                               *
                               * This array of size m is copied internally,
                               * so that the caller can change the incoming data after
                               * return without that IpoptProblem is modified.
                               */
   int          num_nonzeros_jacobian,    /**< Number of non-zero elements in constraint Jacobian */
   int          num_nonzeros_hessian,   /**< Number of non-zero elements in Hessian of Lagrangian */
   Eval_Obj_CB      eval_obj,      /**< Callback function for evaluating objective function */
   Eval_Obj_Grad_CB eval_obj_grad, /**< Callback function for evaluating gradient of objective function */
    Eval_Constr_CB      eval_constr,      /**< Callback function for evaluating constraint functions */
   Eval_Constr_Jac_CB  eval_constr_jac,  /**< Callback function for evaluating Jacobian of constraint functions */
   Eval_Lag_Hess_CB      eval_lag_hess       /**< Callback function for evaluating Hessian of Lagrangian function */
);

/** Method for freeing a previously created IpoptProblem.
 *
 * After freeing an IpoptProblem, it cannot be used anymore.
 */
void FreeSqpProblem(
   SqpProblem sqp_problem
);

/** Function for adding a string option.
 *
 * @return FALSE, if the option could not be set (e.g., if keyword is unknown)
 */
Bool AddSqpStrOption(
   SqpProblem sqp_problem,
   char*        keyword,
   char*        val
);

/** Function for adding a Number option.
 *
 * @return FALSE, if the option could not be set (e.g., if keyword is unknown)
 */
Bool AddSqpNumOption(
   SqpProblem   sqp_problem,
   char*        keyword,
   double       val
);

/** Function for adding an Int option.
 *
 * @return FALSE, if the option  could not be set (e.g., if keyword is unknown)
 @*/
Bool AddSqpIntOption(
   SqpProblem sqp_problem,
   char*        keyword,
   int          val
);

/** Function for opening an output file for a given name with given printlevel.
 *
 * @return FALSE, if there was a problem opening the file.
 */
Bool OpenSqpOutputFile(
   SqpProblem sqp_problem,
   char*        file_name,
   int          print_level
);

/** Function that returns the double array with the lower bounds for the
 *  variables.  The caller must not free this array.  This function can be
 *  used to change the values of the bounds.
 */
double* GetVariableLowerBounds(SqpProblem sqp_problem);

/** Function that returns the double array with the upper bounds for the
 *  variables.  The caller must not free this array.  This function can be
 *  used to change the values of the bounds.
 */
double* GetVariableUpperBounds(SqpProblem sqp_problem);

/** Function that returns the double array with the lower bounds for the
 *  constraints.  The caller must not free this array.  This function can be
 *  used to change the values of the bounds.
 */
double* GetConstraintLowerBounds(SqpProblem sqp_problem);

/** Function that returns the double array with the upper bounds for the
 *  constriants.  The caller must not free this array.  This function can be
 *  used to change the values of the bounds.
 */
double* GetConstraintUpperBounds(SqpProblem sqp_problem);

/**  Function calling the Sqp optimization algorithm for a problem
 *   previously defined with SqpProblem.  This does not assume that
 *   the problem has been solved before.  The user has to provide
 *   initial values for the primal variables.  Values for the dual variables
 *   of the initial active set are optional.  The option "starting_mode"
 *   should be set accordingly.
 *
 * @return outcome of the optimization procedure (e.g., success, failure etc).
 */
enum SqpSolverExitStatus SqpSolve(
   SqpProblem   sqp_problem, /**< Problem that is to be optimized.
                                *
                                * Sqp will use the options previously specified with
                                * AddSqpOption (etc) for this problem.
                                */
   double*      primal_variables,             /**< Input: Starting point; Output: Optimal solution */
   double*      constraint_values,             /**< Values of constraint at final point (output only; ignored if set to NULL) */
   double*      objective_value,       /**< Final value of objective function (output only; ignored if set to NULL) */
   double*      bound_multipliers,      /**< Input: Initial values for the multipliers for lower variable bounds (only if warm start option is chosen);
                                 *  Output: Final multipliers for lower variable bounds (ignored if set to NULL)
                                 */
   double*      constraint_multipliers,        /**< Input: Initial values for the constraint multipliers (only if warm start option is chosen);
                                *  Output: Final multipliers for constraints (ignored if set to NULL)
                                */
   char*        bound_activities,      /**< Input: Initial values for the multipliers for upper variable bounds (only if warm start option is chosen);
                                *  Output: Final multipliers for upper variable bounds (ignored if set to NULL)
                                */
   char*        constraint_activities,
   UserDataPtr  user_data      /**< Pointer to user data.
                                *
                                * This will be passed unmodified to the callback functions.
                                */
);

/**  Function calling the Sqp optimization algorithm for a problem
 *   previously defined with SqpIpoptProblem.  This is the crossover version
 *   that solves the problem with Ipopt first to get very close to the optimal
 *   solution and then reverts to the active-set algorithm to get an optimal
 *   active set that can be used for warm starts later.
 *
 * @return outcome of the optimization procedure (e.g., success, failure etc).
 */
enum SqpSolverExitStatus CrossoverSolve(
   SqpProblem   sqp_problem, /**< Problem that is to be optimized.
                                *
                                * Sqp will use the options previously specified with
                                * AddSqpOption (etc) for this problem.
                                */
   double*      primal_variables,             /**< Input: Starting point; Output: Optimal solution */
   double*      constraint_values,             /**< Values of constraint at final point (output only; ignored if set to NULL) */
   double*      objective_value,       /**< Final value of objective function (output only; ignored if set to NULL) */
   double*      bound_multipliers,      /**< Input: Initial values for the multipliers for lower variable bounds (only if warm start option is chosen);
                                 *  Output: Final multipliers for lower variable bounds (ignored if set to NULL)
                                 */
   double*      constraint_multipliers,        /**< Input: Initial values for the constraint multipliers (only if warm start option is chosen);
                                *  Output: Final multipliers for constraints (ignored if set to NULL)
                                */
   char*        bound_activities,      /**< Input: Initial values for the multipliers for upper variable bounds (only if warm start option is chosen);
                                *  Output: Final multipliers for upper variable bounds (ignored if set to NULL)
                                */
   char*        constraint_activities,
   UserDataPtr  user_data      /**< Pointer to user data.
                                *
                                * This will be passed unmodified to the callback functions.
                                */
);

/**  Function calling the Sqp optimization algorithm for a problem
 *   previously solve with SqpSolve or CrossoverSolve.  The problem can
 *   have changed, as long as the number and order of variables and
 *   constraints have not changed.  This will use the previous primal-dual
 *   solution and working set from the most recent solve as starting point.
 *
 * @return outcome of the optimization procedure (e.g., success, failure etc).
 */
enum SqpSolverExitStatus SqpResolve(
   SqpProblem   sqp_problem, /**< Problem that is to be optimized.
                                *
                                * Sqp will use the options previously specified with
                                * AddSqpOption (etc) for this problem.
                                */
   double*      primal_variables,             /**< Input: Starting point; Output: Optimal solution */
   double*      constraint_values,             /**< Values of constraint at final point (output only; ignored if set to NULL) */
   double*      objective_value,       /**< Final value of objective function (output only; ignored if set to NULL) */
   double*      bound_multipliers,      /**< Input: Initial values for the multipliers for lower variable bounds (only if warm start option is chosen);
                                 *  Output: Final multipliers for lower variable bounds (ignored if set to NULL)
                                 */
   double*      constraint_multipliers,        /**< Input: Initial values for the constraint multipliers (only if warm start option is chosen);
                                *  Output: Final multipliers for constraints (ignored if set to NULL)
                                */
   char*        bound_activities,      /**< Input: Initial values for the multipliers for upper variable bounds (only if warm start option is chosen);
                                *  Output: Final multipliers for upper variable bounds (ignored if set to NULL)
                                */
   char*        constraint_activities,
   UserDataPtr  user_data      /**< Pointer to user data.
                                *
                                * This will be passed unmodified to the callback functions.
                                */
);

#ifdef __cplusplus
} /* extern "C" { */
#endif

#endif
