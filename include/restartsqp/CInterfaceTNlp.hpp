/* Copyright (C) 2020
 * All Rights Reserved.
 * This code is published under the Eclipse Public License.
 *
 * This file is derived from the file IpStdInterfaceNLP.hpp from the Ipopt project
 */

#ifndef __SQP_CINTERFACETNLP_HPP__
#define __SQP_CINTERFACETNLP_HPP__

#include "restartsqp/CInterface.h"
#include "restartsqp/SqpTNlp.hpp"

namespace RestartSqp
{
/** Implementation of a SqpTNlp for the C interface.
 *
 *  The C interface is exposed to the user as a single C
 *  function that is given problem dimension, starting points, and
 *  pointers for functions that evaluate objective function etc.
 */
class CInterfaceTNlp : public SqpTNlp
{
public:
   /**@name Constructors/Destructors */
   //@{
   /** Constructor, given dimensions of problem, function pointers
    *  for evaluation callback functions, and starting points.
    *
    *  Note that the constructor does not make a copy of any of the double
    *  arrays, i.e. it is up to the called to keep them around.
    */
   CInterfaceTNlp(
      int           num_variables,
      const double*   variable_lower_bounds,
      const double*   variable_upper_bounds,
      int           num_constraints,
      const double*   constraint_lower_bounds,
      const double*   constraint_upper_bounds,
      int           num_nonzeros_jacobian,
      int           num_nonzeros_hessian,
      const double*   initial_primal_variables,
      const double*   initial_bound_multipliers,
      const double*   initial_constraint_multipliers,
       const char* initial_bound_activities,
       const char* initial_constraint_activities,
       Eval_Obj_CB      eval_obj,
       Eval_Obj_Grad_CB eval_obj_grad,
       Eval_Constr_CB      eval_constr,
       Eval_Constr_Jac_CB  eval_constr_jac,
       Eval_Lag_Hess_CB      eval_lag_hess,
      double*         primal_variables,
      double*         bound_multipliers,
      double*         constraint_multipliers,
      double*         constraint_values,
      double*         objective_value,
       char* bound_activities,
        char* constraint_activities,
      UserDataPtr     user_data
   );

   /** Default destructor */
   virtual ~CInterfaceTNlp();
   //@}

   /**@name Methods to gather information about the NLP.
    *
    * These methods are overloaded from TNLP. See TNLP for their more detailed documentation.
    */
   //@{
   bool get_nlp_info(
      int&          num_variables,
      int&          num_constraints,
      int&          num_nonzeros_jacobian,
      int&          num_nonzeros_hessian,
      std::string& nlp_name
   ) override;

   bool get_bounds_info(
      int   num_variables,
      double* variable_lower_bounds,
      double* variable_upper_bounds,
      int   num_constraints,
      double* constraint_lower_bounds,
      double* constraint_upper_bounds
   ) override;

   bool get_starting_point(
       int num_variables,
       bool init_primal_variables,
       double* primal_variables,
       bool init_bound_multipliers,
       double* bound_multipliers,
       int num_constraints,
       bool init_constraint_multipliers,
       double* constraint_multipliers
   ) override;

   bool eval_objective_value(int num_variables,
                                           const double* primal_variables,
                                           bool new_primal_variables,
                                           double& objective_value
   ) override;

   bool eval_objective_gradient(int num_variables,
                                          const double* primal_variables,
                                          bool new_primal_variables,
                                          double* objective_gradient
                                ) override;

   bool eval_constraint_values(int num_variables,
                              const double* primal_variables,
                              bool new_primal_variables,
                              int num_constraints,
                              double* constraint_values
   ) override;

   bool eval_constraint_jacobian(int num_variables, const double* primal_variables,
                                 bool new_primal_variables, int num_constraints,
                                 int num_nonzeros_jacobian, int* row_indices,
                                 int* column_indices, double* nonzero_values
   ) override;

   bool eval_lagrangian_hessian(
             int num_variables, const double* primal_variables,
             bool new_primal_variables, double objective_scaling_factor,
             int num_constraints, const double* constraint_multipliers,
             bool new_constraint_multipliers, int num_nonzeros_hessian,
             int* row_indices, int* column_indices, double* nonzero_values
   ) override;

   //@}

   /** @name Solution Methods */
   //@{
   void finalize_solution(
       SqpSolverExitStatus status, int num_variables,
       const double* primal_solution, const double* bound_multipliers,
       const ActivityStatus* bound_activity_status, int num_constraints,
       const double* constraint_values, const double* constraint_multipliers,
       const ActivityStatus* constraint_activity_status, double objective_value,
       std::shared_ptr<const Statistics> stats
   ) override;
   //@}

   bool use_initial_working_set() override;

     /** Set the initial working set. */
     bool get_initial_working_sets(int num_variables,
                                   ActivityStatus* bounds_working_set,
                                   int num_constraints,
                                   ActivityStatus* constraints_working_set);

private:
   /** @name Information about the problem */
   //@{
   /** number of variables */
   const int num_variables_;
   /** number of constraints */
   const int num_constraints_;
   /** Pointer to double array containing lower bounds for variables */
   const double* variable_lower_bounds_;
   /** Pointer to double array containing upper bounds for variables */
   const double* variable_upper_bounds_;
   /** Pointer to double array containing lower bounds for constraints */
   const double* constraint_lower_bounds_;
   /** Pointer to double array containing upper bounds for constraints */
   const double* constraint_upper_bounds_;
   /** double of non-zero elements in the constraint Jacobian */
   const int num_nonzeros_jacobian_;
   /** double of non-zero elements in the Hessian */
   const int num_nonzeros_hessian_;
   /** Pointer to double array containing starting point for variables */
   const double* initial_primal_variables_;
   /** Pointer to double array containing starting values for variable multipliers */
   const double* initial_bound_multipliers_;
   /** Pointer to double array containing starting values for constraint multipliers */
   const double* initial_constraint_multipliers_;
   /** Pointer to char array contining the initial activities for variable bounds. */
   const char* initial_bound_activities_;
   /** Pointer to char array contining the initial activities for constraints. */
   const char* initial_constraint_activities_;
   /** Pointer to callback function evaluating value of objective function */
   Eval_Obj_CB      eval_obj_;
   /** Pointer to callback function evaluating gradient of objective function */
   Eval_Obj_Grad_CB eval_obj_grad_;
   /** Pointer to callback function evaluating value of constraints */
   Eval_Constr_CB      eval_constr_;
   /** Pointer to callback function evaluating Jacobian of constraints */
   Eval_Constr_Jac_CB  eval_constr_jac_;
   /** Pointer to callback function evaluating Hessian of Lagrangian */
   Eval_Lag_Hess_CB      eval_lag_hess_;
   /** Pointer to user data */
   UserDataPtr user_data_;
   //@}

   /** A non-const copy of primal variables - this is kept up-to-date in apply_new_primal_variables */
   double* non_const_primal_variables_;

   /** @name Pointers to the user provided vectors for solution */
   //@{
   double* primal_variables_;
   double* bound_multipliers_;
   double* constraint_multipliers_;
   double* constraint_values_;
   double* objective_value_;
   char* bound_activities_;
   char* constraint_activities_;
   //@}

   /** Update the internal state if the value of primal variables changes */
   void apply_new_primal_variable_(
      bool          new_primal_variables,
      int         num_variables,
      const double* primal_variables
   );

   /**@name Default Compiler Generated Methods
    * (Hidden to avoid implicit creation/calling).
    * These methods are not implemented and
    * we do not want the compiler to implement
    * them for us, so we declare them private
    * and do not define them. This ensures that
    * they will not be implicitly created/called. */
   //@{
   /** Default Constructor */
   CInterfaceTNlp();

   /** Copy Constructor */
   CInterfaceTNlp(
      const CInterfaceTNlp&
   );

   /** Default Assignment Operator */
   void operator=(
      const CInterfaceTNlp&
   );
   //@}

};

} // namespace Ipopt

#endif
