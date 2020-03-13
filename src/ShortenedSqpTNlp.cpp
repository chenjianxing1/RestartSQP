/* Copyright (C) 2019
* All Rights Reserved.
*
* Authors: Xinyi Luo
* Date:2019-07
*/

#include "sqphot/ShortenedSqpTNlp.hpp"
#include "sqphot/Utils.hpp"
#include <cassert>

using namespace std;

namespace RestartSqp {

/** Default constructor*/
ShortenedSqpTNlp::ShortenedSqpTNlp(std::shared_ptr<SqpTNlp> sqp_tnlp)
 : sqp_tnlp_(sqp_tnlp)
 , num_constraints_(-1)
 , constraint_indices_(nullptr)
 , sqp_jac_map_(nullptr)
{
}

/** Default constructor*/
ShortenedSqpTNlp::~ShortenedSqpTNlp()
{
  delete[] constraint_indices_;
  delete[] sqp_jac_map_;
}

void ShortenedSqpTNlp::set_considered_constraints(int num_constraints,
                                                  int* constraint_indices)
{
  num_constraints_ = num_constraints;

  // make a copy of the constraint indices
  delete[] constraint_indices_;
  constraint_indices_ = nullptr;

  constraint_indices_ = new int[num_constraints];
  for (int i = 0; i < num_constraints; ++i) {
    constraint_indices_[i] = constraint_indices[i] - 1;
  }
}

bool ShortenedSqpTNlp::get_nlp_info(int& num_variables, int& num_constraints,
                                    int& num_nonzeros_jacobian,
                                    int& num_nonzeros_hessian,
                                    std::string& nlp_name)
{
  bool retval;
  retval = sqp_tnlp_->get_nlp_info(num_variables, num_orig_constraints_,
                                   num_orig_nonzeros_jacobian_,
                                   num_nonzeros_hessian, nlp_name);
  if (!retval) {
    return retval;
  }

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
    num_nonzeros_jacobian += constraint_chosen[row_indices[i] - 1];
  }

  num_nonzeros_jacobian_ = num_nonzeros_jacobian;

  delete[] constraint_chosen;
  delete[] row_indices;
  delete[] col_indices;

  return retval;
}

/**
 *@brief get the bounds information from the NLP object
 */
bool ShortenedSqpTNlp::get_bounds_info(int num_variabes,
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
bool ShortenedSqpTNlp::get_starting_point(
    int num_variables, bool init_primal_variables, double* primal_variables,
    bool init_bound_multipliers, double* bound_multipliers, int num_constraints,
    bool init_constraint_multipliers, double* constraint_multipliers)
{
  assert(num_constraints == num_constraints_);

  // We need to get space for the original constraint multipliers.
  double* orig_constraint_multipliers = nullptr;
  if (init_constraint_multipliers) {
    assert(constraint_multipliers);
    orig_constraint_multipliers = new double[num_orig_constraints_];
  }

  // Get the data from the original NLP
  bool retval;
  // Ipopt separates the bound multipliers into two parts
  retval = sqp_tnlp_->get_starting_point(
      num_variables, init_primal_variables, primal_variables,
      init_bound_multipliers, bound_multipliers, num_orig_constraints_,
      init_constraint_multipliers, orig_constraint_multipliers);

  if (retval) {
    // Now copy the relevant bounds
    for (int i = 0; i < num_constraints_; ++i) {
      constraint_multipliers[i] =
          orig_constraint_multipliers[constraint_indices_[i]];
    }
  }

  delete[] orig_constraint_multipliers;

  return retval;
}

/**
 *@brief Evaluate the objective value
 */
bool ShortenedSqpTNlp::eval_objective_value(int num_variables,
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
bool ShortenedSqpTNlp::eval_constraint_values(int num_variables,
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
      num_variables, primal_variables, new_primal_variables, num_constraints,
      orig_constraint_values);

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
bool ShortenedSqpTNlp::eval_objective_gradient(int num_variables,
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

bool ShortenedSqpTNlp::eval_constraint_jacobian(
    int num_variables, const double* primal_variables,
    bool new_primal_variables, int num_constraints, int num_nonzeros_jacobian,
    int* row_indices, int* column_indices, double* nonzero_values)
{
  assert(num_constraints == num_constraints_);
  assert(num_nonzeros_jacobian == num_nonzeros_jacobian_);
  assert(!nonzero_values);

  bool retval;

  if (row_indices) {
    // Now we need to provide the structure
    assert(column_indices);

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
        inv_constraint_indices[i] = 0;
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
bool ShortenedSqpTNlp::eval_lagrangian_hessian(
    int num_variables, const double* primal_variables,
    bool new_primal_variables, double objective_scaling_factor,
    int num_constraints, const double* constraint_multipliers,
    bool new_constraint_multipliers, int num_nonzeros_hessian, int* row_indices,
    int* column_indices, double* nonzero_values)
{
  assert(num_constraints == num_constraints_);

  // We need to create a multiplier array that has 0 for the constraints that
  // are not included in this NLP.
  double* orig_constraint_multipliers = new double[num_orig_constraints_];
  for (int i = 0; i < num_orig_constraints_; ++i) {
    orig_constraint_multipliers[i] = 0.;
  }
  // Copy the relevant values
  for (int i = 0; i < num_constraints_; ++i) {
    orig_constraint_multipliers[constraint_indices_[i]] =
        constraint_multipliers[i];
  }

  bool retval = sqp_tnlp_->eval_lagrangian_hessian(
      num_variables, primal_variables, new_primal_variables,
      objective_scaling_factor, num_constraints, orig_constraint_multipliers,
      new_constraint_multipliers, num_nonzeros_hessian, row_indices,
      column_indices, nonzero_values);

  delete[] orig_constraint_multipliers;

  return retval;
}

void ShortenedSqpTNlp::finalize_solution(
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

  delete[] orig_constraint_values;
  delete[] orig_constraint_multipliers;
  delete[] orig_constraint_activity_status;
}
}
