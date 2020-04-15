/* Copyright (C) 2019
* All Rights Reserved.
*
* Authors: Xinyi Luo
* Date:2019-07
*/

#include "restartsqp/LazySqpTNlp.hpp"
#include "restartsqp/Utils.hpp"
#include <cassert>

using namespace std;

namespace RestartSqp {

/** Default constructor*/
LazySqpTNlp::LazySqpTNlp(std::shared_ptr<SqpTNlp> sqp_tnlp)
 : sqp_tnlp_(sqp_tnlp)
 , num_constraints_(-1)
 , constraint_indices_(nullptr)
 , sqp_jac_map_(nullptr)
 , bound_activity_status_(nullptr)
 , constraint_activity_status_(nullptr)
 , has_been_solved_before_(false)
 , previous_optimal_solution_(nullptr)
 , previous_optimal_bound_multipliers_(nullptr)
 , previous_optimal_constraint_multipliers_(nullptr)
{
}

/** Default constructor*/
LazySqpTNlp::~LazySqpTNlp()
{
  delete[] constraint_indices_;
  delete[] sqp_jac_map_;
  delete[] bound_activity_status_;
  delete[] constraint_activity_status_;
  delete[] previous_optimal_solution_;
  delete[] previous_optimal_bound_multipliers_;
  delete[] previous_optimal_constraint_multipliers_;
}

void LazySqpTNlp::set_considered_constraints(int num_constraints,
                                             const int* constraint_indices)
{
  num_constraints_ = num_constraints;

  // make a copy of the constraint indices
  delete[] constraint_indices_;
  constraint_indices_ = nullptr;

  constraint_indices_ = new int[num_constraints_];
  for (int i = 0; i < num_constraints_; ++i) {
    constraint_indices_[i] = constraint_indices[i];
  }
}

bool LazySqpTNlp::get_nlp_info(int& num_variables, int& num_constraints,
                               int& num_nonzeros_jacobian,
                               int& num_nonzeros_hessian, std::string& nlp_name)
{
  bool retval;
  retval = sqp_tnlp_->get_nlp_info(num_variables_, num_orig_constraints_,
                                   num_orig_nonzeros_jacobian_,
                                   num_nonzeros_hessian, nlp_name);
  if (!retval) {
    return retval;
  }

  num_variables = num_variables_;
  num_constraints = num_constraints_;

  // Get space to store the jacobian structure.
  int* row_indices = new int[num_orig_nonzeros_jacobian_];
  int* col_indices = new int[num_orig_nonzeros_jacobian_];

  // Ask for the structure of the original Jacobian
  retval = sqp_tnlp_->eval_constraint_jacobian(
      num_variables, NULL, false, num_orig_constraints_,
      num_orig_nonzeros_jacobian_, row_indices, col_indices, NULL);
  if (!retval) {
    delete[] row_indices;
    delete[] col_indices;
  }

  // Create an array that has 1 for all constraints that are in this NLP
  char* constraint_chosen = new char[num_orig_constraints_];
  for (int i = 0; i < num_orig_constraints_; ++i) {
    constraint_chosen[i] = 0;
  }
  for (int i = 0; i < num_constraints_; ++i) {
    constraint_chosen[constraint_indices_[i]] = 1;
  }

  // Now count the number of nonzeros in the shortened Jacobian
  num_nonzeros_jacobian = 0;
  for (int i = 0; i < num_orig_nonzeros_jacobian_; ++i) {
    num_nonzeros_jacobian += constraint_chosen[row_indices[i]];
  }

  num_nonzeros_jacobian_ = num_nonzeros_jacobian;

  delete[] constraint_chosen;
  delete[] row_indices;
  delete[] col_indices;

  // Allocate memory for the working sets
  if (!bound_activity_status_) {
    assert(!constraint_activity_status_);
    bound_activity_status_ = new ActivityStatus[num_variables_];
    constraint_activity_status_ = new ActivityStatus[num_orig_constraints_];
  }

  return retval;
}

/**
 *@brief get the bounds information from the NLP object
 */
bool LazySqpTNlp::get_bounds_info(int num_variabes,
                                  double* variable_lower_bounds,
                                  double* variable_upper_bounds,
                                  int num_constraints,
                                  double* constraint_lower_bounds,
                                  double* constraint_upper_bounds)
{
  assert(num_constraints == num_constraints_);

  // We need to get space for the original bounds
  double* orig_constraint_lower_bounds = new double[num_orig_constraints_];
  double* orig_constraint_upper_bounds = new double[num_orig_constraints_];

  bool retval;
  retval = sqp_tnlp_->get_bounds_info(
      num_variabes, variable_lower_bounds, variable_upper_bounds,
      num_orig_constraints_, orig_constraint_lower_bounds,
      orig_constraint_upper_bounds);

  if (retval) {
    // Now copy the relevant bounds
    for (int i = 0; i < num_constraints_; ++i) {
      constraint_lower_bounds[i] =
          orig_constraint_lower_bounds[constraint_indices_[i]];
      constraint_upper_bounds[i] =
          orig_constraint_upper_bounds[constraint_indices_[i]];
    }
  }

  delete[] orig_constraint_lower_bounds;
  delete[] orig_constraint_upper_bounds;

  return retval;
}

/*
 * @brief Get the starting point from the NLP object.
 */
bool LazySqpTNlp::get_starting_point(
    int num_variables, bool init_primal_variables, double* primal_variables,
    bool init_bound_multipliers, double* bound_multipliers, int num_constraints,
    bool init_constraint_multipliers, double* constraint_multipliers)
{
  assert(num_constraints == num_constraints_);

  bool retval;

  if (!has_been_solved_before_) {
    // If this is the first solve and no initial multipliers are available, get
    // the original
    // starting point from the user.

    // We need to get space for the original constraint multipliers.
    double* orig_constraint_multipliers = nullptr;
    if (init_constraint_multipliers) {
      assert(constraint_multipliers);
      orig_constraint_multipliers = new double[num_orig_constraints_];
    }

    // Get the data from the original NLP
    retval = sqp_tnlp_->get_starting_point(
        num_variables, init_primal_variables, primal_variables,
        init_bound_multipliers, bound_multipliers, num_orig_constraints_,
        init_constraint_multipliers, orig_constraint_multipliers);

    if (retval && init_constraint_multipliers) {
      // Now copy the relevant bounds
      for (int i = 0; i < num_constraints_; ++i) {
        constraint_multipliers[i] =
            orig_constraint_multipliers[constraint_indices_[i]];
      }
    }

    delete[] orig_constraint_multipliers;
  } else {
    // If an NLP has already been solved, take the previous optimal solution
    // as a starting point.
    if (init_primal_variables) {
      for (int i = 0; i < num_variables_; ++i) {
        primal_variables[i] = previous_optimal_solution_[i];
        bound_multipliers[i] = previous_optimal_bound_multipliers_[i];
      }
    }
    if (init_bound_multipliers) {
      for (int i = 0; i < num_variables_; ++i) {
        primal_variables[i] = previous_optimal_solution_[i];
        bound_multipliers[i] = previous_optimal_bound_multipliers_[i];
      }
    }
    if (init_constraint_multipliers) {
      for (int i = 0; i < num_constraints_; ++i) {
        constraint_multipliers[i] =
            previous_optimal_constraint_multipliers_[constraint_indices_[i]];
      }
    }
  }

  return retval;
}

/**
 *@brief Evaluate the objective value
 */
bool LazySqpTNlp::eval_objective_value(int num_variables,
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
bool LazySqpTNlp::eval_constraint_values(int num_variables,
                                         const double* primal_variables,
                                         bool new_primal_variables,
                                         int num_constraints,
                                         double* constraint_values)
{
  assert(num_constraints == num_constraints_);

  // We need to get space for the original constraint values.
  double* orig_constraint_values = new double[num_orig_constraints_];

  bool retval;
  retval = sqp_tnlp_->eval_constraint_values(
      num_variables, primal_variables, new_primal_variables,
      num_orig_constraints_, orig_constraint_values);

  if (retval) {
    for (int i = 0; i < num_constraints_; ++i) {
      constraint_values[i] = orig_constraint_values[constraint_indices_[i]];
    }
  }

  delete[] orig_constraint_values;

  return retval;
}

/**
 *@brief Evaluate gradient at point x
 */
bool LazySqpTNlp::eval_objective_gradient(int num_variables,
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

bool LazySqpTNlp::eval_constraint_jacobian(
    int num_variables, const double* primal_variables,
    bool new_primal_variables, int num_constraints, int num_nonzeros_jacobian,
    int* row_indices, int* column_indices, double* nonzero_values)
{
  assert(num_constraints == num_constraints_);
  assert(num_nonzeros_jacobian == num_nonzeros_jacobian_);

  bool retval;

  if (row_indices) {
    // Now we need to provide the structure
    assert(column_indices);
    assert(!nonzero_values);

    // Get space for the original Jacobian structure
    int* orig_row_indices = new int[num_orig_nonzeros_jacobian_];
    int* orig_col_indices = new int[num_orig_nonzeros_jacobian_];

    retval = sqp_tnlp_->eval_constraint_jacobian(
        num_variables, primal_variables, new_primal_variables, num_constraints,
        num_orig_nonzeros_jacobian_, orig_row_indices, orig_col_indices,
        nonzero_values);

    if (retval) {
      // Compute the inverse mapping for constraint_indices
      int* inv_constraint_indices = new int[num_orig_constraints_];
      for (int i = 0; i < num_orig_constraints_; ++i) {
        inv_constraint_indices[i] = -1;
      }
      for (int i = 0; i < num_constraints_; ++i) {
        inv_constraint_indices[constraint_indices_[i]] = i;
      }

      // Create a map that records the position of the relevant nonzero elements
      // in the original arrays.
      delete[] sqp_jac_map_;
      sqp_jac_map_ = nullptr;

      // Create the map and copy the relevant entries
      sqp_jac_map_ = new int[num_nonzeros_jacobian];
      int nnz = 0;
      for (int i = 0; i < num_orig_nonzeros_jacobian_; ++i) {
        int row = orig_row_indices[i];
        if (inv_constraint_indices[row] >= 0) {
          sqp_jac_map_[nnz] = i;
          row_indices[nnz] = inv_constraint_indices[row];
          column_indices[nnz] = orig_col_indices[i];
          nnz++;
        }
      }
      assert(nnz == num_nonzeros_jacobian_);

      delete[] inv_constraint_indices;
    }

    delete[] orig_row_indices;
    delete[] orig_col_indices;
  } else {
    // In this case we need to get the nonzero values
    assert(nonzero_values);
    assert(!row_indices);
    assert(!column_indices);

    // First let's get space for the original values
    double* orig_nonzero_values = new double[num_orig_nonzeros_jacobian_];

    retval = sqp_tnlp_->eval_constraint_jacobian(
        num_variables, primal_variables, new_primal_variables, num_constraints,
        num_orig_nonzeros_jacobian_, row_indices, column_indices,
        orig_nonzero_values);

    if (retval) {
      // Copy the relevant values
      for (int i = 0; i < num_nonzeros_jacobian_; ++i) {
        nonzero_values[i] = orig_nonzero_values[sqp_jac_map_[i]];
      }
    }

    delete[] orig_nonzero_values;
  }

  return retval;
}

/**
 * @brief Get the structure of the Hessian
 * Always call this before the first time using @Eval_Hessian
 */
bool LazySqpTNlp::eval_lagrangian_hessian(
    int num_variables, const double* primal_variables,
    bool new_primal_variables, double objective_scaling_factor,
    int num_constraints, const double* constraint_multipliers,
    bool new_constraint_multipliers, int num_nonzeros_hessian, int* row_indices,
    int* column_indices, double* nonzero_values)
{
  assert(num_constraints == num_constraints_);

  double* orig_constraint_multipliers = nullptr;
  if (nonzero_values) {
    assert(!row_indices);
    // If we are computing the Hessian values, in order to ignore the unselected
    // constraints,
    // we need to create a multiplier array that has 0 for the constraints that
    // are not included in this NLP.
    orig_constraint_multipliers = new double[num_orig_constraints_];
    for (int i = 0; i < num_orig_constraints_; ++i) {
      orig_constraint_multipliers[i] = 0.;
    }
    // Copy the relevant values
    for (int i = 0; i < num_constraints_; ++i) {
      orig_constraint_multipliers[constraint_indices_[i]] =
          constraint_multipliers[i];
    }
  }

  bool retval = sqp_tnlp_->eval_lagrangian_hessian(
      num_variables, primal_variables, new_primal_variables,
      objective_scaling_factor, num_orig_constraints_,
      orig_constraint_multipliers, new_constraint_multipliers,
      num_nonzeros_hessian, row_indices, column_indices, nonzero_values);

  delete[] orig_constraint_multipliers;

  return retval;
}

void LazySqpTNlp::finalize_solution(
    SqpSolverExitStatus status, int num_variables,
    const double* primal_solution, const double* bound_multipliers,
    const ActivityStatus* bound_activity_status, int num_constraints,
    const double* constraint_values, const double* constraint_multipliers,
    const ActivityStatus* constraint_activity_status, double objective_value,
    std::shared_ptr<const Statistics> stats)
{
  assert(num_constraints == num_constraints_);

  // We pass the constraint information in the original format.
  double* orig_constraint_values = new double[num_orig_constraints_];
  double* orig_constraint_multipliers = new double[num_orig_constraints_];
  ActivityStatus* orig_constraint_activity_status =
      new ActivityStatus[num_orig_constraints_];

  // We set dummy values for those entries that are not included in this NLP.
  for (int i = 0; i < num_orig_constraints_; ++i) {
    orig_constraint_values[i] = SqpNaN;
    orig_constraint_multipliers[i] = SqpNaN;
    orig_constraint_activity_status[i] = ACTIVE_INVALID;
  }

  // Copy the relevant values
  for (int i = 0; i < num_constraints_; ++i) {
    int idx = constraint_indices_[i];
    orig_constraint_values[idx] = constraint_values[i];
    orig_constraint_multipliers[idx] = constraint_multipliers[i];
    orig_constraint_activity_status[idx] = constraint_activity_status[i];
  }

  // Call the finalization method
  sqp_tnlp_->finalize_solution(
      status, num_variables, primal_solution, bound_multipliers,
      bound_activity_status, num_orig_constraints_, orig_constraint_values,
      orig_constraint_multipliers, orig_constraint_activity_status,
      objective_value, stats);
  assert(num_variables == num_variables_);
  assert(num_constraints == num_constraints_);

  delete[] orig_constraint_values;
  delete[] orig_constraint_multipliers;
  delete[] orig_constraint_activity_status;

  // If the solve was not successful, we do not store the solution internally
  // and return here
  if (status != OPTIMAL) {
    return;
  }

  // Store the activity status for the variables and constraints
  for (int i = 0; i < num_variables; i++) {
    bound_activity_status_[i] = bound_activity_status[i];
  }
  for (int i = 0; i < num_constraints; i++) {
    constraint_activity_status_[constraint_indices_[i]] =
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
        new double[num_orig_constraints_];
  }

  // Initialize all constraint multipliers to 0.  When a new constraints is
  // added, this
  // will be the initial value.
  for (int i = 0; i < num_orig_constraints_; ++i) {
    previous_optimal_constraint_multipliers_[i] = 0.;
  }

  // Copy the primal and dual solutions
  for (int i = 0; i < num_variables_; ++i) {
    previous_optimal_solution_[i] = primal_solution[i];
    previous_optimal_bound_multipliers_[i] = bound_multipliers[i];
  }
  for (int i = 0; i < num_constraints; ++i) {
    previous_optimal_constraint_multipliers_[constraint_indices_[i]] =
        constraint_multipliers[i];
  }

  // Copy the optimal solution to
}

bool LazySqpTNlp::get_initial_working_sets(
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
        constraint_activity_status_[constraint_indices_[i]];
  }

  return true;
}

bool LazySqpTNlp::add_new_constraints(
    int num_new_considered_constraints,
    const int* new_considered_constraints_indices)
{
  assert(num_constraints_ + num_new_considered_constraints <=
         num_orig_constraints_);

  // We need to increase the size of the constraint_indices_ array.  Let's put
  // it's current content into a backup.
  int* backup = new int[num_constraints_];
  copy(constraint_indices_, constraint_indices_ + num_constraints_, backup);

  // increase the size of the member array
  delete[] constraint_indices_;
  constraint_indices_ = nullptr;
  constraint_indices_ =
      new int[num_constraints_ + num_new_considered_constraints];

  // copy the original entries back
  copy(backup, backup + num_constraints_, constraint_indices_);

  // Don't need the backup anymore
  delete[] backup;

  // add the new indices at the end, correcting the offer
  for (int i = 0; i < num_new_considered_constraints; ++i) {
    constraint_indices_[num_constraints_ + i] =
        new_considered_constraints_indices[i];
  }

  // Also mark the added constraints as inactive
  for (int i = num_constraints_;
       i < num_constraints_ + num_new_considered_constraints; ++i) {
    constraint_activity_status_[constraint_indices_[i]] = INACTIVE;
  }

  // increase the number of constraints considered
  num_constraints_ += num_new_considered_constraints;

  return true;
}
}
