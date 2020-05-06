/* Copyright (C) 2019
* All Rights Reserved.
*
* Authors: Xinyi Luo
* Date:2019-07
*/
#include "restartsqp/SqpNlp.hpp"

using namespace std;

namespace RestartSqp {

/** Default constructor*/
SqpNlp::SqpNlp(std::shared_ptr<SqpTNlp> sqp_tnlp)
 : sqp_tnlp_(sqp_tnlp)
{
  sqp_tnlp_->get_nlp_info(num_variables_, num_constraints_,
                          num_nonzeros_jacobian_, num_nonzeros_hessian_,
                          nlp_name_);
}

/** Default constructor*/
SqpNlp::~SqpNlp()
{
}

shared_ptr<const SqpNlpSizeInfo> SqpNlp::get_problem_sizes()
{
  shared_ptr<SqpNlpSizeInfo> retval = make_shared<SqpNlpSizeInfo>(
      num_variables_, num_constraints_, num_nonzeros_jacobian_,
      num_nonzeros_hessian_);
  return retval;
}

/**
 *@brief get the bounds information from the NLP object
 */
bool SqpNlp::get_bounds_info(shared_ptr<Vector> x_l, shared_ptr<Vector> x_u,
                             shared_ptr<Vector> c_l, shared_ptr<Vector> c_u)
{
  // Large value for bounds that should be interpreted as infinity
  const double NLP_INF = 1e18;

  double* x_L = x_l->get_non_const_values();
  double* x_U = x_u->get_non_const_values();
  double* c_L = c_l->get_non_const_values();
  double* c_U = c_u->get_non_const_values();

  bool retval;
  retval = sqp_tnlp_->get_bounds_info(num_variables_, x_L, x_U,
                                      num_constraints_, c_L, c_U);
  if (!retval) {
    return false;
  }

  // We translate very large bounds into true Inf values
  for (int i = 0; i < num_variables_; ++i) {
    if (x_L[i] <= -NLP_INF) {
      x_L[i] = -SqpInf;
    }
    if (x_U[i] >= NLP_INF) {
      x_U[i] = SqpInf;
    }
  }
  for (int i = 0; i < num_constraints_; ++i) {
    if (c_L[i] <= -NLP_INF) {
      c_L[i] = -SqpInf;
    }
    if (c_U[i] >= NLP_INF) {
      c_U[i] = SqpInf;
    }
  }
  return true;
}

bool SqpNlp::use_initial_working_set()
{
  return sqp_tnlp_->use_initial_working_set();
}

bool SqpNlp::get_initial_working_sets(int num_variables,
                                      ActivityStatus* bounds_working_set,
                                      int num_constraints,
                                      ActivityStatus* constraints_working_set)
{
  assert(num_variables == num_variables_);
  assert(num_constraints == num_constraints_);
  return sqp_tnlp_->get_initial_working_sets(num_variables, bounds_working_set,
                                             num_constraints,
                                             constraints_working_set);
}

/*
 * @brief Get the starting point from the NLP object.
 */
bool SqpNlp::get_starting_point(shared_ptr<Vector> primal_point,
                                shared_ptr<Vector> bound_multipliers,
                                shared_ptr<Vector> constraint_multipliers)
{
  return sqp_tnlp_->get_starting_point(
      num_variables_, true, primal_point->get_non_const_values(), true,
      bound_multipliers->get_non_const_values(), num_constraints_, true,
      constraint_multipliers->get_non_const_values());
}

/**
 *@brief Evaluate the objective value
 */
bool SqpNlp::eval_f(shared_ptr<const Vector> x, double& obj_value)
{
  return sqp_tnlp_->eval_objective_value(num_variables_, x->get_values(), true,
                                         obj_value);
}

/**
 * @brief Evaluate the constraints at point x
 *
 */
bool SqpNlp::eval_constraints(shared_ptr<const Vector> x,
                              shared_ptr<Vector> constraints)
{
  return sqp_tnlp_->eval_constraint_values(num_variables_, x->get_values(),
                                           true, num_constraints_,
                                           constraints->get_non_const_values());
}

/**
 *@brief Evaluate gradient at point x
 */
bool SqpNlp::eval_gradient(shared_ptr<const Vector> x,
                           shared_ptr<Vector> gradient)
{
  return sqp_tnlp_->eval_objective_gradient(
      num_variables_, x->get_values(), true, gradient->get_non_const_values());
}

/**
 * @brief Get the matrix structure of the Jacobian
 * Always call this before the first time using @Eval_Jacobian
 */

bool SqpNlp::get_jacobian_structure(shared_ptr<const Vector> x,
                                    shared_ptr<SparseTripletMatrix> Jacobian)
{
  return sqp_tnlp_->eval_constraint_jacobian(
      num_variables_, x->get_values(), true, num_constraints_,
      num_nonzeros_jacobian_, Jacobian->get_nonconst_row_indices(),
      Jacobian->get_nonconst_column_indices(), NULL);
}

/**
 *@brief Evaluate Jacobian at point x
 */

bool SqpNlp::eval_jacobian(shared_ptr<const Vector> x,
                           shared_ptr<SparseTripletMatrix> Jacobian)
{
  return sqp_tnlp_->eval_constraint_jacobian(
      num_variables_, x->get_values(), true, num_constraints_,
      num_nonzeros_jacobian_, NULL, NULL, Jacobian->get_nonconst_values());
}

/**
 * @brief Get the structure of the Hessian
 * Always call this before the first time using @Eval_Hessian
 */
bool SqpNlp::get_hessian_structure(shared_ptr<const Vector> x,
                                   shared_ptr<const Vector> lambda,
                                   shared_ptr<SparseTripletMatrix> Hessian)
{
  return sqp_tnlp_->eval_lagrangian_hessian(
      num_variables_, x->get_values(), true, 1.0, num_constraints_,
      lambda->get_values(), true, num_nonzeros_hessian_,
      Hessian->get_nonconst_row_indices(),
      Hessian->get_nonconst_column_indices(), NULL);
}

/**
 *@brief Evaluate Hessian of Lagragian function at  (x, lambda)
 */
bool SqpNlp::eval_hessian(shared_ptr<const Vector> x,
                          shared_ptr<const Vector> lambda,
                          double objective_scaling_factor,
                          shared_ptr<SparseTripletMatrix> Hessian)
{
  return sqp_tnlp_->eval_lagrangian_hessian(
      num_variables_, x->get_values(), true, objective_scaling_factor, num_constraints_,
      lambda->get_values(), true, num_nonzeros_hessian_, NULL, NULL,
      Hessian->get_nonconst_values());
}

void SqpNlp::finalize_solution(SqpSolverExitStatus status,
                               shared_ptr<const Vector> primal_solution,
                               shared_ptr<const Vector> bound_multipliers,
                               const ActivityStatus* bound_activity_status,
                               shared_ptr<const Vector> constraint_values,
                               shared_ptr<const Vector> constraint_multipliers,
                               const ActivityStatus* constraint_activity_status,
                               double objective_value,
                               shared_ptr<const Statistics> stats)
{
  int n = primal_solution->get_dim();
  const double* x = primal_solution->get_values();
  const double* bound_mults = bound_multipliers->get_values();

  // Get constraint values
  int m = constraint_values->get_dim();
  const double* constr_vals = constraint_values->get_values();
  const double* constr_mults = constraint_multipliers->get_values();

  // Create copy of constraint multipliers with reversed sign
  // Call the Ipopt method
  sqp_tnlp_->finalize_solution(
      status, n, x, bound_mults, bound_activity_status, m, constr_vals,
      constr_mults, constraint_activity_status, objective_value, stats);
}
}
