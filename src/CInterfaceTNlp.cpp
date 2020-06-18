// Copyright (C) 2004, 2010 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-09-02

#include "restartsqp/CInterfaceTNlp.hpp"

/// @todo use more memcpy

namespace RestartSqp
{

CInterfaceTNlp::CInterfaceTNlp(
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
)
   : num_variables_(num_variables),
     num_constraints_(num_constraints),
     variable_lower_bounds_(variable_lower_bounds),
     variable_upper_bounds_(variable_upper_bounds),
     constraint_lower_bounds_(constraint_lower_bounds),
     constraint_upper_bounds_(constraint_upper_bounds),
     num_nonzeros_jacobian_(num_nonzeros_jacobian),
     num_nonzeros_hessian_(num_nonzeros_hessian),
     initial_primal_variables_(initial_primal_variables),
     initial_bound_multipliers_(initial_bound_multipliers),
     initial_constraint_multipliers_(initial_constraint_multipliers),
     initial_bound_activities_(initial_bound_activities),
     initial_constraint_activities_(initial_constraint_activities),
     eval_obj_(eval_obj),
     eval_obj_grad_(eval_obj_grad),
     eval_constr_(eval_constr),
     eval_constr_jac_(eval_constr_jac),
     eval_lag_hess_(eval_lag_hess),
     user_data_(user_data),
     non_const_primal_variables_(NULL),
     primal_variables_(primal_variables),
     bound_multipliers_(bound_multipliers),
     constraint_multipliers_(constraint_multipliers),
     constraint_values_(constraint_values),
     objective_value_(objective_value),
     bound_activities_(bound_activities),
     constraint_activities_(constraint_activities)
{
#if 0 // add sanity checks
   ASSERT_EXCEPTION(n_var_ > 0, INVALID_STDINTERFACE_NLP, "The number of variables must be at least 1.");
   ASSERT_EXCEPTION(n_con_ >= 0, INVALID_STDINTERFACE_NLP, "The number of constrains must be non-negative.");
   ASSERT_EXCEPTION(x_L_, INVALID_STDINTERFACE_NLP, "No lower bounds for variables provided.");
   ASSERT_EXCEPTION(x_U_, INVALID_STDINTERFACE_NLP, "No upper bounds for variables provided.");
   ASSERT_EXCEPTION(g_L_ || n_con_ == 0, INVALID_STDINTERFACE_NLP, "No lower bounds for constraints provided.");
   ASSERT_EXCEPTION(g_U_ || n_con_ == 0, INVALID_STDINTERFACE_NLP, "No upper bounds for constraints provided.");
   ASSERT_EXCEPTION(nele_jac_ >= 0, INVALID_STDINTERFACE_NLP, "Number of non-zero elements in constraint Jacobian must be non-negative.");
   ASSERT_EXCEPTION(nele_hess_ >= 0, INVALID_STDINTERFACE_NLP, "Number of non-zero elements in Hessian of Lagrangian must be non-negative.");
   ASSERT_EXCEPTION(index_style_ == 0 || index_style_ == 1, INVALID_STDINTERFACE_NLP, "Valid index styles are 0 (C style) or 1 (Fortran style)");
   ASSERT_EXCEPTION(start_x_, INVALID_STDINTERFACE_NLP, "No initial point for the variables provided.");
   ASSERT_EXCEPTION(eval_f_, INVALID_STDINTERFACE_NLP, "No callback function for evaluating the value of objective function provided.");
   ASSERT_EXCEPTION(eval_g_, INVALID_STDINTERFACE_NLP, "No callback function for evaluating the values of constraints provided.");
   ASSERT_EXCEPTION(eval_grad_f_, INVALID_STDINTERFACE_NLP, "No callback function for evaluating the gradient of objective function provided.");
   ASSERT_EXCEPTION(eval_jac_g_, INVALID_STDINTERFACE_NLP, "No callback function for evaluating the Jacobian of the constraints provided.");
   ASSERT_EXCEPTION(eval_h_, INVALID_STDINTERFACE_NLP, "No callback function for evaluating the Hessian of the constraints provided.");
#endif
}

CInterfaceTNlp::~CInterfaceTNlp()
{
   delete[] non_const_primal_variables_;
}

bool CInterfaceTNlp::get_nlp_info(
    int&          num_variables,
    int&          num_constraints,
    int&          num_nonzeros_jacobian,
    int&          num_nonzeros_hessian,
    std::string& nlp_name
)
{
   num_variables = num_variables_;
   num_constraints = num_constraints_;
   num_nonzeros_jacobian = num_nonzeros_jacobian_;
   num_nonzeros_hessian = num_nonzeros_hessian_;

   nlp_name = "CInterfaceNLP";

   return true;
}

bool CInterfaceTNlp::get_bounds_info(
    int   num_variables,
    double* variable_lower_bounds,
    double* variable_upper_bounds,
    int   num_constraints,
    double* constraint_lower_bounds,
    double* constraint_upper_bounds
)
{
   assert(num_variables == num_variables_);
   assert(num_constraints == num_constraints_);

   std::copy(variable_lower_bounds_,
             variable_lower_bounds_+num_variables,
             variable_lower_bounds);
   std::copy(variable_upper_bounds_,
             variable_upper_bounds_+num_variables,
             variable_upper_bounds);

   std::copy(constraint_lower_bounds_,
             constraint_lower_bounds_+num_constraints,
             constraint_lower_bounds);
   std::copy(constraint_upper_bounds_,
             constraint_upper_bounds_+num_constraints,
             constraint_upper_bounds);

   return true;
}


bool CInterfaceTNlp::get_starting_point(
    int num_variables,
    bool init_primal_variables,
    double* primal_variables,
    bool init_bound_multipliers,
    double* bound_multipliers,
    int num_constraints,
    bool init_constraint_multipliers,
    double* constraint_multipliers
)
{
   assert(num_variables == num_variables_);
   assert(num_constraints == num_constraints_);

   if( init_primal_variables ) {
     std::copy(initial_primal_variables_,
               initial_primal_variables_+num_variables,
               primal_variables);
   }

   if( init_bound_multipliers ) {
     if (initial_bound_multipliers_ == nullptr) {
       return false;
     }
     std::copy(initial_bound_multipliers_,
               initial_bound_multipliers_+num_variables,
               bound_multipliers);
   }

   if( init_constraint_multipliers ) {
     if (initial_constraint_multipliers_ == nullptr) {
       return false;
     }
     std::copy(initial_constraint_multipliers_,
               initial_constraint_multipliers_+num_constraints,
               constraint_multipliers);
   }

   return true;
}

bool CInterfaceTNlp::eval_objective_value(int num_variables,
                                          const double* primal_variables,
                                          bool new_primal_variables,
                                          double& objective_value
)
{
   assert(num_variables == num_variables_);

   apply_new_primal_variable_(new_primal_variables, num_variables, primal_variables);

   Bool retval = (*eval_obj_)(num_variables, non_const_primal_variables_, (Bool) new_primal_variables, &objective_value, user_data_);

   return (retval != 0);
}

bool CInterfaceTNlp::eval_objective_gradient(int num_variables,
                                             const double* primal_variables,
                                             bool new_primal_variables,
                                             double* objective_gradient
)
{
  assert(num_variables == num_variables_);

  apply_new_primal_variable_(new_primal_variables, num_variables, primal_variables);

   Bool retval = (*eval_obj_grad_)(num_variables, non_const_primal_variables_, (Bool) new_primal_variables, objective_gradient, user_data_);

   return (retval != 0);
}

bool CInterfaceTNlp::eval_constraint_values(int num_variables,
                                            const double* primal_variables,
                                            bool new_primal_variables,
                                            int num_constraints,
                                            double* constraint_values
)
{
  assert(num_variables == num_variables_);
  assert(num_constraints == num_constraints_);

  apply_new_primal_variable_(new_primal_variables, num_variables, primal_variables);

   Bool retval = (*eval_constr_)(num_variables, non_const_primal_variables_, (Bool) new_primal_variables, num_constraints, constraint_values, user_data_);

   return (retval != 0);
}

bool CInterfaceTNlp::eval_constraint_jacobian(
    int num_variables, const double* primal_variables,
                                     bool new_primal_variables, int num_constraints,
                                     int num_nonzeros_jacobian, int* row_indices,
                                     int* column_indices, double* nonzero_values)
{
  assert(num_variables == num_variables_);
  assert(num_constraints == num_constraints_);
  assert(num_nonzeros_jacobian == num_nonzeros_jacobian_);

   // need valid combination of iRow, jCol, and values pointers
   assert((row_indices != NULL && column_indices != NULL && nonzero_values == NULL) || (row_indices == NULL && column_indices == NULL && nonzero_values != NULL));

   apply_new_primal_variable_(new_primal_variables, num_variables, primal_variables);
   Bool retval = (*eval_constr_jac_)(num_variables, non_const_primal_variables_, (Bool) new_primal_variables, num_constraints, num_nonzeros_jacobian, row_indices, column_indices, nonzero_values, user_data_);

   return (retval != 0);
}

bool CInterfaceTNlp::eval_lagrangian_hessian(
    int num_variables, const double* primal_variables,
    bool new_primal_variables, double objective_scaling_factor,
    int num_constraints, const double* constraint_multipliers,
    bool new_constraint_multipliers, int num_nonzeros_hessian,
    int* row_indices, int* column_indices, double* nonzero_values
)
{
  assert(num_variables == num_variables_);
  assert(num_constraints == num_constraints_);
  assert(num_nonzeros_hessian == num_nonzeros_hessian_);

  // need valid combination of iRow, jCol, and values pointers
  assert((row_indices != NULL && column_indices != NULL && nonzero_values == NULL) || (row_indices == NULL && column_indices == NULL && nonzero_values != NULL));

  apply_new_primal_variable_(new_primal_variables, num_variables, primal_variables);

   double* non_const_constr_mult =  NULL;
   if( constraint_multipliers ) {
     non_const_constr_mult = new double[num_constraints];
     std::copy(constraint_multipliers,
               constraint_multipliers + num_constraints,
               non_const_constr_mult);
   }
   Bool retval = (*eval_lag_hess_)(num_variables, non_const_primal_variables_, (Bool) new_primal_variables, objective_scaling_factor, num_constraints, non_const_constr_mult, (Bool) new_constraint_multipliers, num_nonzeros_hessian, row_indices, column_indices, nonzero_values, user_data_);

   delete[] non_const_constr_mult;

   return (retval != 0);
}

bool CInterfaceTNlp::use_initial_working_set()
{
  if (initial_bound_activities_) {
    return true;
  }
  return false;
}

bool CInterfaceTNlp::get_initial_working_sets(int num_variables,
                                              ActivityStatus* bounds_working_set,
                                              int num_constraints,
                                              ActivityStatus* constraints_working_set)
{
  assert(num_variables == num_variables_);
  assert(num_constraints == num_constraints_);

  if (initial_bound_activities_ == NULL || (num_constraints > 0 && initial_constraint_activities_ == NULL)) {
    return false;
  }

  for (int i=0; i<num_variables; ++i) {
    bounds_working_set[i] = (ActivityStatus) initial_bound_activities_[i];
  }

  for (int i=0; i<num_constraints; ++i) {
    constraints_working_set[i] = (ActivityStatus) initial_constraint_activities_[i];
  }

  return true;
}

void CInterfaceTNlp::finalize_solution(
    SqpSolverExitStatus status, int num_variables,
    const double* primal_solution, const double* bound_multipliers,
    const ActivityStatus* bound_activity_status, int num_constraints,
    const double* constraint_values, const double* constraint_multipliers,
    const ActivityStatus* constraint_activity_status, double objective_value,
    std::shared_ptr<const Statistics> stats
)
{
  assert(num_variables == num_variables_);
  assert(num_constraints == num_constraints_);

  if ( primal_variables_ ) {
    std::copy(primal_solution, primal_solution+num_variables, primal_variables_);
  }

  if ( bound_multipliers_ ) {
    std::copy(bound_multipliers, bound_multipliers+num_variables,
              bound_multipliers_);
  }

  if ( constraint_multipliers_ ) {
    std::copy(constraint_multipliers, constraint_multipliers+num_constraints,
              constraint_multipliers_);
  }

  if (bound_activities_) {
    for (int i=0; i<num_variables; ++i) {
      bound_activities_[i] = (char)bound_activity_status[i];
    }
  }

  if (constraint_activities_) {
    for (int i=0; i<num_constraints; ++i) {
      constraint_activities_[i] = (char)constraint_activity_status[i];
    }
  }

  if( objective_value_ != NULL )
   {
      *objective_value_ = objective_value;
   }

   // don't need to store the status, we get the status from the OptimizeTNLP method
}

void CInterfaceTNlp::apply_new_primal_variable_(
    bool          new_primal_variables,
    int         num_variables,
    const double* primal_variables
)
{
   if( new_primal_variables )
   {
      assert(primal_variables != NULL);

      //copy the data to the non_const_x_
      if( !non_const_primal_variables_ )
      {

         non_const_primal_variables_ = new double[num_variables];
      }

      std::copy(primal_variables, primal_variables+num_variables,
                non_const_primal_variables_);

   }
}

} // namespace RestartSqp
