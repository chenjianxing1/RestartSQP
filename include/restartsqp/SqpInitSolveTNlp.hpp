/* Copyright (C) 2019
 * All Rights Reserved.
 *
 * Authors: Xinyi Luo
 * Date:2019-06
 */
#ifndef SQPHOTSTART_SQPINITSOLVETNLP_HPP
#define SQPHOTSTART_SQPINITSOLVETNLP_HPP

#include "restartsqp/SqpTNlp.hpp"

namespace RestartSqp {

/**
 * This is part of SQPhotstart
 *
 * This class is used in the SqpRestart::initial_solve() method.
 * It is a wrapper around a given SqpTNlp that can overwrite the initial working
 * set and the starting point.
 *
 */
class SqpInitSolveTNlp : public SqpTNlp
{

public:
  /** @brief Default Constructor */
  SqpInitSolveTNlp(std::shared_ptr<SqpTNlp> orig_sqp_tnlp)
   : orig_sqp_tnlp_(orig_sqp_tnlp)
   , initial_bounds_working_set_(nullptr)
   , initial_constraints_working_set_(nullptr)
   , initial_primal_variables_(nullptr)
   , initial_bound_multipliers_(nullptr)
   , initial_constraint_multipliers_(nullptr)
  {}

  /** Default destructor*/
  virtual ~SqpInitSolveTNlp()
  {
    delete[] initial_bounds_working_set_;
    delete[] initial_constraints_working_set_;
    delete[] initial_primal_variables_;
    delete[] initial_bound_multipliers_;
    delete[] initial_constraint_multipliers_;
  }

  /** Get NLP info from original SqpTNlp. */
  bool get_nlp_info(int& num_variables, int& num_constraints,
                    int& num_nonzeros_jacobian, int& num_nonzeros_hessian,
                    std::string& nlp_name)
  {
    return orig_sqp_tnlp_->get_nlp_info(num_variables, num_constraints,
                                        num_nonzeros_jacobian,
                                        num_nonzeros_hessian, nlp_name);
  }

  /** Get bounds from original SqpTNlp */
  bool get_bounds_info(int num_variabes, double* variable_lower_bounds,
                       double* variable_upper_bounds, int num_constraints,
                       double* constraint_lower_bounds,
                       double* constraint_upper_bounds)
  {
    return orig_sqp_tnlp_->get_bounds_info(
        num_variabes, variable_lower_bounds, variable_upper_bounds,
        num_constraints, constraint_lower_bounds, constraint_upper_bounds);
  }

  /** This method returns the starting point that previously have been set with
   * set_starting_point.
   */
  bool get_starting_point(int num_variables, bool init_primal_variables,
                          double* initial_primal_variables,
                          bool init_bound_multipliers,
                          double* initial_bound_multipliers,
                          int num_constraints, bool init_constraint_multipliers,
                          double* initial_constraint_multipliers)
  {
    if (init_primal_variables) {

      if (!initial_primal_variables_) {
        return false;
      } else {
        for (int i = 0; i < num_variables; ++i) {
          initial_primal_variables[i] = initial_primal_variables_[i];
        }
      }
    }

    if (init_bound_multipliers) {
      if (!initial_bound_multipliers_) {
        return false;
      } else {
        for (int i = 0; i < num_variables; ++i) {
          initial_bound_multipliers[i] = initial_bound_multipliers_[i];
        }
      }
    }

    if (init_constraint_multipliers) {
      if (!initial_constraint_multipliers_) {
        return false;
      } else {
        for (int i = 0; i < num_constraints; ++i) {
          initial_constraint_multipliers[i] =
              initial_constraint_multipliers_[i];
        }
      }
    }

    return true;
  }

  /** Return objective function value from original SqpTNlp. */
  bool eval_objective_value(int num_variables, const double* primal_variables,
                            bool new_primal_variables, double& objective_value)
  {
    return orig_sqp_tnlp_->eval_objective_value(
        num_variables, primal_variables, new_primal_variables, objective_value);
  }

  /** Return objective gradient from original SqpTNlp. */
  bool eval_objective_gradient(int num_variables,
                               const double* primal_variables,
                               bool new_primal_variables,
                               double* objective_gradient)
  {
    return orig_sqp_tnlp_->eval_objective_gradient(
        num_variables, primal_variables, new_primal_variables,
        objective_gradient);
  }

  /** Return constraint values from original SqpTNlp. */
  bool eval_constraint_values(int num_variables, const double* primal_variables,
                              bool new_primal_variables, int num_constraints,
                              double* constraint_values)
  {
    return orig_sqp_tnlp_->eval_constraint_values(
        num_variables, primal_variables, new_primal_variables, num_constraints,
        constraint_values);
  }

  /** Return constraint Jaobian structure of values from original SqpTNlp. */
  bool eval_constraint_jacobian(int num_variables,
                                const double* primal_variables,
                                bool new_primal_variables, int num_constraints,
                                int num_nonzeros_jacobian, int* row_indices,
                                int* column_indices, double* nonzero_values)
  {
    return orig_sqp_tnlp_->eval_constraint_jacobian(
        num_variables, primal_variables, new_primal_variables, num_constraints,
        num_nonzeros_jacobian, row_indices, column_indices, nonzero_values);
  }

  /** Return Lagrangian Hessian structure of values from original SqpTNlp.  */
  bool eval_lagrangian_hessian(
      int num_variables, const double* primal_variables,
      bool new_primal_variables, double objective_scaling_factor,
      int num_constraints, const double* constraint_multipliers,
      bool new_constraint_multipliers, int num_nonzeros_hessian,
      int* row_indices, int* column_indices, double* nonzero_values)
  {
    return orig_sqp_tnlp_->eval_lagrangian_hessian(
        num_variables, primal_variables, new_primal_variables,
        objective_scaling_factor, num_constraints, constraint_multipliers,
        new_constraint_multipliers, num_nonzeros_hessian, row_indices,
        column_indices, nonzero_values);
  }

  //@}

  /** @name Solution Methods */
  //@{
  /** Call the finalize_solution method of the original SqpTNlp so that the
   * caller can retrieve the solution. */
  void finalize_solution(SqpSolverExitStatus status, int num_variables,
                         const double* primal_solution,
                         const double* bound_multipliers,
                         const ActivityStatus* bound_activity_status,
                         int num_constraints, const double* constraint_values,
                         const double* constraint_multipliers,
                         const ActivityStatus* constraint_activity_status,
                         double objective_value,
                         std::shared_ptr<const Statistics> stats)
  {
    orig_sqp_tnlp_->finalize_solution(
        status, num_variables, primal_solution, bound_multipliers,
        bound_activity_status, num_constraints, constraint_values,
        constraint_multipliers, constraint_activity_status, objective_value,
        stats);
  }

  //@}

  /** Set the initial working set. */
  bool set_initial_working_sets(int num_variables,
                                const ActivityStatus* bounds_working_set,
                                int num_constraints,
                                const ActivityStatus* constraints_working_set)
  {
    // Delete any previous working sets
    delete[] initial_bounds_working_set_;
    initial_bounds_working_set_ = nullptr;
    delete[] initial_constraints_working_set_;
    initial_constraints_working_set_ = nullptr;

    // Get the memory for the new working sets
    initial_bounds_working_set_ = new ActivityStatus[num_variables];
    initial_constraints_working_set_ = new ActivityStatus[num_constraints];

    // Copy the working sets
    for (int i = 0; i < num_variables; ++i) {
      initial_bounds_working_set_[i] = bounds_working_set[i];
    }
    for (int i = 0; i < num_constraints; ++i) {
      initial_constraints_working_set_[i] = constraints_working_set[i];
    }

    return true;
  }

  /** Set the initial working set. */
  bool get_initial_working_sets(int num_variables,
                                ActivityStatus* bounds_working_set,
                                int num_constraints,
                                ActivityStatus* constraints_working_set)
  {
    for (int i = 0; i < num_variables; ++i) {
      bounds_working_set[i] = initial_bounds_working_set_[i];
    }
    for (int i = 0; i < num_constraints; ++i) {
      constraints_working_set[i] = initial_constraints_working_set_[i];
    }
    return true;
  }

  /** Return true since a working set should have been set. */
  bool use_initial_working_set() override
  {
    return true;
  }

  /** Retrieve the overwriting starting point. */
  void set_starting_point(int num_variables,
                          const double* initial_primal_variables,
                          const double* initial_bound_multipliers,
                          int num_constraints,
                          const double* initial_constraint_multipliers)
  {
    if (initial_primal_variables) {
      initial_primal_variables_ = new double[num_variables];
      for (int i = 0; i < num_variables; ++i) {
        initial_primal_variables_[i] = initial_primal_variables[i];
      }
    }

    if (initial_bound_multipliers) {
      initial_bound_multipliers_ = new double[num_variables];
      for (int i = 0; i < num_variables; ++i) {
        initial_bound_multipliers_[i] = initial_bound_multipliers[i];
      }
    }

    if (initial_constraint_multipliers) {
      initial_constraint_multipliers_ = new double[num_constraints];
      for (int i = 0; i < num_constraints; ++i) {
        initial_constraint_multipliers_[i] = initial_constraint_multipliers[i];
      }
    }
  }

private:
  /** Copy Constructor */
  SqpInitSolveTNlp(const SqpInitSolveTNlp&);

  /** Overloaded Equals Operator */
  void operator=(const SqpInitSolveTNlp&);

  /** The original SqpTNlp. */
  std::shared_ptr<SqpTNlp> orig_sqp_tnlp_;

  /** @name Initial working sets */
  //@{
  /** Working set for the variable bounds. */
  ActivityStatus* initial_bounds_working_set_;
  /** Working set for the contraint bounds. */
  ActivityStatus* initial_constraints_working_set_;
  //@}

  /** @name Starting point. */
  //@{
  /** Initial primal starting point. */
  double* initial_primal_variables_;
  /** Initial bound multipliers. */
  double* initial_bound_multipliers_;
  /** Initial constraint multipliers. */
  double* initial_constraint_multipliers_;
  //@}
};
} // namespace RestartSqp

#endif // SQPHOTSTART_SQPINITSOLVETNLP_HPP
