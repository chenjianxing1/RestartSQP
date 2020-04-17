/* Copyright (C) 2019
* All Rights Reserved.
*
* Authors: Xinyi Luo
* Date:2019-07
*/
#include "restartsqp/SqpIpoptTNlp.hpp"
#include <cassert>

using namespace std;
using namespace Ipopt;

namespace RestartSqp {

/** Default constructor*/
SqpIpoptTNlp::SqpIpoptTNlp(SmartPtr<TNLP> ipopt_tnlp, const string& nlp_name)
 : ipopt_tnlp_(ipopt_tnlp)
 , nlp_name_(nlp_name)
{
}

/** Default constructor*/
SqpIpoptTNlp::~SqpIpoptTNlp()
{
}

bool SqpIpoptTNlp::get_nlp_info(int& num_variables, int& num_constraints,
                               int& num_nonzeros_jacobian,
                               int& num_nonzeros_hessian, std::string& nlp_name)
{
  bool retval;
  retval = ipopt_tnlp_->get_nlp_info(num_variables, num_constraints,
                                     num_nonzeros_jacobian,
                                     num_nonzeros_hessian, index_style_);
  nlp_name = nlp_name_;

  return retval;
}

/**
 *@brief get the bounds information from the NLP object
 */
bool SqpIpoptTNlp::get_bounds_info(int num_variabes,
                                  double* variable_lower_bounds,
                                  double* variable_upper_bounds,
                                  int num_constraints,
                                  double* constraint_lower_bounds,
                                  double* constraint_upper_bounds)
{
  bool retval;
  retval = ipopt_tnlp_->get_bounds_info(
      num_variabes, variable_lower_bounds, variable_upper_bounds,
      num_constraints, constraint_lower_bounds, constraint_upper_bounds);
  return retval;
}

/*
 * @brief Get the starting point from the NLP object.
 */
bool SqpIpoptTNlp::get_starting_point(
    int num_variables, bool init_primal_variables, double* primal_variables,
    bool init_bound_multipliers, double* bound_multipliers, int num_constraints,
    bool init_constraint_multipliers, double* constraint_multipliers)
{
  //
  // TODO: There is a problem with getting bound multipliers in the Ipopt AMPL
  // interface (somewhat related to the suffix handler).  This might be an Ipopt
  // bug; need to look into this?
  const bool bug_fixed = false;

  bool retval;
  double* z_U = NULL;
  bool init_z_call = init_bound_multipliers;

  if (!bug_fixed) {
    init_z_call = false;
  }

  if (init_z_call) {
    z_U = new double[num_variables];
  }

  // Ipopt separates the bound multipliers into two parts
  retval = ipopt_tnlp_->get_starting_point(
      num_variables, init_primal_variables, primal_variables, init_z_call,
      bound_multipliers, z_U, num_constraints, init_constraint_multipliers,
      constraint_multipliers);

  // We need to collapse the bound multipliers into a single vector
  if (init_z_call) {
    for (int i = 0; i < num_variables; ++i) {
      bound_multipliers[i] -= z_U[i];
    }
  } else if (init_bound_multipliers) {
    // In this case we didn't get bound multipliers from the call
    for (int i = 0; i < num_variables; ++i) {
      bound_multipliers[i] = 0.;
    }
  }

  delete[] z_U;

  if (!retval) {
    return retval;
  }

  // We also need to reverse the sign of the constraint multipliers
  if (init_constraint_multipliers) {
    for (int i = 0; i < num_constraints; ++i) {
      constraint_multipliers[i] = -constraint_multipliers[i];
    }
  }

  return retval;
}

/**
 *@brief Evaluate the objective value
 */
bool SqpIpoptTNlp::eval_objective_value(int num_variables,
                                       const double* primal_variables,
                                       bool new_primal_variables,
                                       double& objective_value)
{
  return ipopt_tnlp_->eval_f(num_variables, primal_variables,
                             new_primal_variables, objective_value);
}

/**
 *@brief Evaluate gradient at point x
 */
bool SqpIpoptTNlp::eval_objective_gradient(int num_variables,
                                          const double* primal_variables,
                                          bool new_primal_variables,
                                          double* objective_gradient)
{
  return ipopt_tnlp_->eval_grad_f(num_variables, primal_variables,
                                  new_primal_variables, objective_gradient);
}

/**
 * @brief Evaluate the constraints at point x
 *
 */
bool SqpIpoptTNlp::eval_constraint_values(int num_variables,
                                         const double* primal_variables,
                                         bool new_primal_variables,
                                         int num_constraints,
                                         double* constraint_values)
{
  return ipopt_tnlp_->eval_g(num_variables, primal_variables,
                             new_primal_variables, num_constraints,
                             constraint_values);
}

/**
 * @brief Get the matrix structure of the Jacobian
 * Always call this before the first time using @Eval_Jacobian
 */

bool SqpIpoptTNlp::eval_constraint_jacobian(
    int num_variables, const double* primal_variables,
    bool new_primal_variables, int num_constraints, int num_nonzeros_jacobian,
    int* row_indices, int* column_indices, double* nonzero_values)
{
  bool retval = ipopt_tnlp_->eval_jac_g(
      num_variables, primal_variables, new_primal_variables, num_constraints,
      num_nonzeros_jacobian, row_indices, column_indices, nonzero_values);

  // Adjust numbering if necessary;
  if (retval && row_indices && index_style_ == TNLP::FORTRAN_STYLE) {
    for (int i=0; i<num_nonzeros_jacobian; ++i) {
      row_indices[i]-=1;
      column_indices[i]-=1;
    }
  }

  return retval;
}

/**
 * @brief Get the structure of the Hessian
 * Always call this before the first time using @Eval_Hessian
 */
bool SqpIpoptTNlp::eval_lagrangian_hessian(
    int num_variables, const double* primal_variables,
    bool new_primal_variables, double objective_scaling_factor,
    int num_constraints, const double* constraint_multipliers,
    bool new_constraint_multipliers, int num_nonzeros_hessian, int* row_indices,
    int* column_indices, double* nonzero_values)
{
  // For Ipopt, the constraint multipliers are defined with the opposite sign
  // So we need to create a copy of constraint multipliers with reversed sign
  double* lambda = nullptr;
  if (constraint_multipliers) {
    lambda = new double[num_constraints];
    for (int i = 0; i < num_constraints; ++i) {
      lambda[i] = -constraint_multipliers[i];
    }
  }

  bool retval =
      ipopt_tnlp_->eval_h(num_variables, primal_variables, new_primal_variables,
                          objective_scaling_factor, num_constraints, lambda,
                          new_constraint_multipliers, num_nonzeros_hessian,
                          row_indices, column_indices, nonzero_values);

  delete[] lambda;

  // Adjust numbering if necessary;
  if (retval && row_indices && index_style_ == TNLP::FORTRAN_STYLE) {
    for (int i=0; i<num_nonzeros_hessian; ++i) {
      row_indices[i]-=1;
      column_indices[i]-=1;
    }
  }

  return retval;
}

void SqpIpoptTNlp::finalize_solution(
    SqpSolverExitStatus status, int num_variables,
    const double* primal_solution, const double* bound_multipliers,
    const ActivityStatus* bound_activity_status, int num_constraints,
    const double* constraint_values, const double* constraint_multipliers,
    const ActivityStatus* constraint_activity_status, double objective_value,
    std::shared_ptr<const Statistics> stats)
{
  // Translate SQP exit status to Ipopt exit status
  SolverReturn ipopt_status;
  switch (status) {
    case OPTIMAL:
      ipopt_status = SUCCESS;
      break;
    case EXCEED_MAX_ITERATIONS:
      ipopt_status = MAXITER_EXCEEDED;
      break;
    case TRUST_REGION_TOO_SMALL:
      ipopt_status = STOP_AT_TINY_STEP;
      break;
    case EXCEED_MAX_CPU_TIME:
      ipopt_status = CPUTIME_EXCEEDED;
      break;
    case EXCEED_MAX_WALLCLOCK_TIME:
      ipopt_status = CPUTIME_EXCEEDED;
      break;
    case PENALTY_TOO_LARGE:
      ipopt_status = RESTORATION_FAILURE;
      break;
    default:
      ipopt_status = INTERNAL_ERROR;
      break;
  }

  // Separate the bound multipliers into two parts
  double* z_L = new double[num_variables];
  double* z_U = new double[num_variables];
  for (int i = 0; i < num_variables; ++i) {
    z_L[i] = max(0., bound_multipliers[i]);
    z_U[i] = max(0., -bound_multipliers[i]);
  }

  // Create copy of constraint multipliers with reversed sign
  double* lambda = new double[num_constraints];
  for (int i = 0; i < num_constraints; ++i) {
    lambda[i] = -constraint_multipliers[i];
  }

  // There are not Ipopt-specific objects
  const IpoptData* ip_data = NULL;
  IpoptCalculatedQuantities* ip_cq = NULL;

  // Call the Ipopt method
  ipopt_tnlp_->finalize_solution(ipopt_status, num_variables, primal_solution,
                                 z_L, z_U, num_constraints, constraint_values,
                                 lambda, objective_value, ip_data, ip_cq);

  delete[] z_L;
  delete[] z_U;
  delete[] lambda;
}

}
