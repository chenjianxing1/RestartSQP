/* Copyright (C) 2019
* All Rights Reserved.
*
* Authors: Xinyi Luo
* Date:2019-07
*/

#include "restartsqp/CrossoverSqpTNlp.hpp"
#include "restartsqp/Utils.hpp"
#include "restartsqp/SqpIpoptTNlp.hpp"
#include <cassert>

using namespace std;

namespace RestartSqp {

CrossoverSqpTNlp::CrossoverSqpTNlp(std::shared_ptr<SqpTNlp> sqp_tnlp)
 : sqp_tnlp_(sqp_tnlp)
 , bound_activity_status_(nullptr)
 , constraint_activity_status_(nullptr)
 , has_been_solved_before_(false)
 , previous_optimal_solution_(nullptr)
 , previous_optimal_bound_multipliers_(nullptr)
 , previous_optimal_constraint_multipliers_(nullptr)
{
}

CrossoverSqpTNlp::CrossoverSqpTNlp(Ipopt::SmartPtr<Ipopt::TNLP> ipopt_tnlp)
 : sqp_tnlp_(make_shared<SqpIpoptTNlp>(ipopt_tnlp))
 , bound_activity_status_(nullptr)
 , constraint_activity_status_(nullptr)
 , has_been_solved_before_(false)
 , previous_optimal_solution_(nullptr)
 , previous_optimal_bound_multipliers_(nullptr)
 , previous_optimal_constraint_multipliers_(nullptr)
{
}

/** Default constructor*/
CrossoverSqpTNlp::~CrossoverSqpTNlp()
{
  delete[] bound_activity_status_;
  delete[] constraint_activity_status_;
  delete[] previous_optimal_solution_;
  delete[] previous_optimal_bound_multipliers_;
  delete[] previous_optimal_constraint_multipliers_;
}

bool CrossoverSqpTNlp::get_nlp_info(int& num_variables, int& num_constraints,
                               int& num_nonzeros_jacobian,
                               int& num_nonzeros_hessian, std::string& nlp_name)
{
  bool retval;
  retval = sqp_tnlp_->get_nlp_info(num_variables, num_constraints,
                                   num_nonzeros_jacobian,
                                   num_nonzeros_hessian, nlp_name);
  if (!retval) {
    return retval;
  }

  if (has_been_solved_before_) {
    if (num_variables != num_variables_ || num_constraints != num_constraints_) {
      printf("The number of variables or constraints has changed since the last time.");
      assert(num_variables == num_variables_ && num_constraints == num_constraints_);
    }
  }
  else {
    num_variables_ = num_variables;
    num_constraints_ = num_constraints;
  }
  return retval;
}

/**
 *@brief get the bounds information from the NLP object
 */
bool CrossoverSqpTNlp::get_bounds_info(int num_variabes,
                                  double* variable_lower_bounds,
                                  double* variable_upper_bounds,
                                  int num_constraints,
                                  double* constraint_lower_bounds,
                                  double* constraint_upper_bounds)
{
  assert(num_constraints == num_constraints_);

  bool retval;
  retval = sqp_tnlp_->get_bounds_info(
      num_variabes, variable_lower_bounds, variable_upper_bounds,
      num_constraints, constraint_lower_bounds,
      constraint_upper_bounds);

  return retval;
}

/*
 * @brief Get the starting point from the NLP object.
 */
bool CrossoverSqpTNlp::get_starting_point(
    int num_variables, bool init_primal_variables, double* primal_variables,
    bool init_bound_multipliers, double* bound_multipliers, int num_constraints,
    bool init_constraint_multipliers, double* constraint_multipliers)
{
  assert(num_constraints == num_constraints_);
  bool retval;

  if (!has_been_solved_before_) {
    // If this is the first solve and no initial multipliers are available, get
    // the starting point from the user.

    // Get the data from the original NLP
    retval = sqp_tnlp_->get_starting_point(
        num_variables, init_primal_variables, primal_variables,
        init_bound_multipliers, bound_multipliers, num_constraints,
        init_constraint_multipliers, constraint_multipliers);
  } else {
    // If an NLP has already been solved, take the previous optimal solution
    // as a starting point.
    if (init_primal_variables) {
      for (int i = 0; i < num_variables_; ++i) {
        primal_variables[i] = previous_optimal_solution_[i];
      }
    }
    if (init_bound_multipliers) {
      for (int i = 0; i < num_variables_; ++i) {
        bound_multipliers[i] = previous_optimal_bound_multipliers_[i];
      }
    }
    if (init_constraint_multipliers) {
      for (int i = 0; i < num_constraints_; ++i) {
        constraint_multipliers[i] = previous_optimal_constraint_multipliers_[i];
      }
    }
  }

  return retval;
}

/**
 *@brief Evaluate the objective value
 */
bool CrossoverSqpTNlp::eval_objective_value(int num_variables,
                                       const double* primal_variables,
                                       bool new_primal_variables,
                                       double& objective_value)
{
  return sqp_tnlp_->eval_objective_value(num_variables, primal_variables,
                                         new_primal_variables, objective_value);
}

/**
 * @brief Evaluate the constraints at point x
 *
 */
bool CrossoverSqpTNlp::eval_constraint_values(int num_variables,
                                         const double* primal_variables,
                                         bool new_primal_variables,
                                         int num_constraints,
                                         double* constraint_values)
{
  return sqp_tnlp_->eval_constraint_values(
      num_variables, primal_variables, new_primal_variables,
      num_constraints, constraint_values);
}

/**
 *@brief Evaluate gradient at point x
 */
bool CrossoverSqpTNlp::eval_objective_gradient(int num_variables,
                                          const double* primal_variables,
                                          bool new_primal_variables,
                                          double* objective_gradient)
{
  return sqp_tnlp_->eval_objective_gradient(num_variables, primal_variables,
                                            new_primal_variables,
                                            objective_gradient);
}

/**
 * @brief Get the matrix structure of the Jacobian
 * Always call this before the first time using @Eval_Jacobian
 */

bool CrossoverSqpTNlp::eval_constraint_jacobian(
    int num_variables, const double* primal_variables,
    bool new_primal_variables, int num_constraints, int num_nonzeros_jacobian,
    int* row_indices, int* column_indices, double* nonzero_values)
{
  return sqp_tnlp_->eval_constraint_jacobian(
        num_variables, primal_variables, new_primal_variables, num_constraints,
        num_nonzeros_jacobian, row_indices, column_indices,
        nonzero_values);
}

/**
 * @brief Get the structure of the Hessian
 * Always call this before the first time using @Eval_Hessian
 */
bool CrossoverSqpTNlp::eval_lagrangian_hessian(
    int num_variables, const double* primal_variables,
    bool new_primal_variables, double objective_scaling_factor,
    int num_constraints, const double* constraint_multipliers,
    bool new_constraint_multipliers, int num_nonzeros_hessian, int* row_indices,
    int* column_indices, double* nonzero_values)
{
  return sqp_tnlp_->eval_lagrangian_hessian(
      num_variables, primal_variables, new_primal_variables,
      objective_scaling_factor, num_constraints,
      constraint_multipliers, new_constraint_multipliers,
      num_nonzeros_hessian, row_indices, column_indices, nonzero_values);
}

void CrossoverSqpTNlp::finalize_solution(
    SqpSolverExitStatus status, int num_variables,
    const double* primal_solution, const double* bound_multipliers,
    const ActivityStatus* bound_activity_status, int num_constraints,
    const double* constraint_values, const double* constraint_multipliers,
    const ActivityStatus* constraint_activity_status, double objective_value,
    std::shared_ptr<const Statistics> stats)
{
  // Call the original finalization method
  sqp_tnlp_->finalize_solution(
      status, num_variables, primal_solution, bound_multipliers,
      bound_activity_status, num_constraints, constraint_values,
      constraint_multipliers, constraint_activity_status,
      objective_value, stats);

  // Allocate memory for the working sets
  if (!bound_activity_status_) {
    assert(!constraint_activity_status_);
    bound_activity_status_ = new ActivityStatus[num_variables_];
    constraint_activity_status_ = new ActivityStatus[num_constraints_];
  }

  // Store the activity status for the variables and constraints
  for (int i = 0; i < num_variables; i++) {
    bound_activity_status_[i] = bound_activity_status[i];
  }
  for (int i = 0; i < num_constraints; i++) {
    constraint_activity_status_[i] =
        constraint_activity_status[i];
  }

  // If this is the first solved problem, allocate memory for the starting point
  if (!has_been_solved_before_) {
    has_been_solved_before_ = true;
    assert(!previous_optimal_solution_);
    assert(!previous_optimal_bound_multipliers_);
    assert(!previous_optimal_constraint_multipliers_);

    previous_optimal_solution_ = new double[num_variables_];
    previous_optimal_bound_multipliers_ = new double[num_variables_];
    previous_optimal_constraint_multipliers_ =
        new double[num_constraints_];
  }

  // Copy the primal and dual solutions
  for (int i = 0; i < num_variables_; ++i) {
    previous_optimal_solution_[i] = primal_solution[i];
    previous_optimal_bound_multipliers_[i] = bound_multipliers[i];
  }
  for (int i = 0; i < num_constraints; ++i) {
    previous_optimal_constraint_multipliers_[i] =
        constraint_multipliers[i];
  }
}

bool CrossoverSqpTNlp::use_initial_working_set()
{
  if (bound_activity_status_) {
    assert(constraint_activity_status_);
    return true;
  }
  return false;
}

bool CrossoverSqpTNlp::get_initial_working_sets(
    int num_variables, ActivityStatus* bounds_working_set, int num_constraints,
    ActivityStatus* constraints_working_set)
{
  assert(num_constraints == num_constraints_);
  assert(num_variables == num_variables_);

  for (int i = 0; i < num_variables; i++) {
    bounds_working_set[i] = bound_activity_status_[i];
  }
  for (int i = 0; i < num_constraints; i++) {
    constraints_working_set[i] =
        constraint_activity_status_[i];
  }

  return true;
}

}
