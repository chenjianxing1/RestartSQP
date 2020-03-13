/* Copyright (C) 2019
 * All Rights Reserved.
 *
 * Authors: Xinyi Luo
 * Date:2019-06
 */
#ifndef SQPHOTSTART_SQPTNLP_HPP
#define SQPHOTSTART_SQPTNLP_HPP

#include "sqphot/Statistics.hpp"
#include "sqphot/Types.hpp"
#include <memory>

namespace RestartSqp {

enum SqpSolverExitStatus
{
  OPTIMAL = 0,
  INVALID_NLP = -1,
  CONVERGE_TO_NONOPTIMAL = 1,
  EXCEED_MAX_ITERATIONS = -2, // exceeds the maximum number of iteration
  PRED_REDUCTION_NEGATIVE = -3,
  TRUST_REGION_TOO_SMALL = -4,
  EXCEED_MAX_CPU_TIME = -6,
  EXCEED_MAX_WALLCLOCK_TIME = -7,
  PENALTY_TOO_LARGE = -8,
  QPERROR_INTERNAL_ERROR = -21, // QP solver internal error
  QPERROR_INFEASIBLE =
      -22, // QP solver error: conclude QP formulation infeasible
  QPERROR_UNBOUNDED = -23,       // QP solver error: unbounded QP
  QPERROR_EXCEED_MAX_ITER = -24, // QP solver error: Exceed maximum iteration,
  QPERROR_NOTINITIALISED = -25,
  QPERROR_PREPARINGAUXILIARYQP = -26,
  QPERROR_AUXILIARYQPSOLVED = -27,
  QPERROR_PERFORMINGHOMOTOPY = -28,
  QPERROR_HOMOTOPYQPSOLVED = -29,
  QPERROR_UNKNOWN = -30,
  UNKNOWN_EXIT_STATUS = -99 // unknown error
};

/**
 * This is part of SQPhotstart
 *
 * This class is the base class for the classes that a user has to provide.
 *
 * Unless the corresponding methods are overloaded, this baseclass can store
 * initial working sets.  This is necessary so that the result of the crossover
 * can be stored.
 *
 */
class SqpTNlp
{

public:
  /** @brief Default Constructor */
  SqpTNlp()
   : initial_bounds_working_set_(nullptr)
   , initial_constraints_working_set_(nullptr)
   , overwrite_starting_point_(false)
   , initial_primal_variables_(nullptr)
   , initial_bound_multipliers_(nullptr)
   , initial_constraint_multipliers_(nullptr)
  {
  }

  /** Default destructor*/
  virtual ~SqpTNlp()
  {
    delete[] initial_bounds_working_set_;
    delete[] initial_constraints_working_set_;
    delete[] initial_primal_variables_;
    delete[] initial_bound_multipliers_;
    delete[] initial_constraint_multipliers_;
  }

  /** Method that returns the dimensions of the problem: double of variables,
   * constraints, nonzeros in the Jacobian and nonzeros on the Hessian. */
  virtual bool get_nlp_info(int& num_variables, int& num_constraints,
                            int& num_nonzeros_jacobian,
                            int& num_nonzeros_hessian,
                            std::string& nlp_name) = 0;

  /** overload this method to return the information about the bound
   *  on the variables and constraints. The value that indicates
   *  that a bound does not exist is specified in the parameters
   *  nlp_lower_bound_inf and nlp_upper_bound_inf.  By default,
   *  nlp_lower_bound_inf is -1e19 and nlp_upper_bound_inf is
   *  1e19. (see TNLPAdapter) */
  virtual bool get_bounds_info(int num_variabes, double* variable_lower_bounds,
                               double* variable_upper_bounds,
                               int num_constraints,
                               double* constraint_lower_bounds,
                               double* constraint_upper_bounds) = 0;

  /** overload this method to return the starting point. The bool
   *  variables indicate whether the algorithm wants you to
   *  initialize x, z_L/z_u, and lambda, respectively.  If, for some
   *  reason, the algorithm wants you to initialize these and you
   *  cannot, return false, which will cause Ipopt to stop.  You
   *  will have to run Ipopt with different options then.
   */
  virtual bool get_starting_point(int num_variables, bool init_primal_variables,
                                  double* primal_variables,
                                  bool init_bound_multipliers,
                                  double* bound_multipliers,
                                  int num_constraints,
                                  bool init_constraint_multipliers,
                                  double* constraint_multipliers) = 0;

  /** overload this method to return the value of the objective function */
  virtual bool eval_objective_value(int num_variables,
                                    const double* primal_variables,
                                    bool new_primal_variables,
                                    double& objective_value) = 0;

  /** overload this method to return the vector of the gradient of
   *  the objective w.r.t. x */
  virtual bool eval_objective_gradient(int num_variables,
                                       const double* primal_variables,
                                       bool new_primal_variables,
                                       double* objective_gradient) = 0;

  /** overload this method to return the vector of constraint values */
  virtual bool eval_constraint_values(int num_variables,
                                      const double* primal_variables,
                                      bool new_primal_variables,
                                      int num_constraints,
                                      double* constraint_values) = 0;

  /** overload this method to return the jacobian of the
   *  constraints. The vectors iRow and jCol only need to be set
   *  once. The first call is used to set the structure only (iRow
   *  and jCol will be non-NULL, and values will be NULL) For
   *  subsequent calls, iRow and jCol will be NULL. */
  virtual bool
  eval_constraint_jacobian(int num_variables, const double* primal_variables,
                           bool new_primal_variables, int num_constraints,
                           int num_nonzeros_jacobian, int* row_indices,
                           int* column_indices, double* nonzero_values) = 0;

  /** overload this method to return the hessian of the
   *  lagrangian. The vectors iRow and jCol only need to be set once
   *  (during the first call). The first call is used to set the
   *  structure only (iRow and jCol will be non-NULL, and values
   *  will be NULL) For subsequent calls, iRow and jCol will be
   *  NULL. This matrix is symmetric - specify the lower diagonal
   *  only.  A default implementation is provided, in case the user
   *  wants to se quasi-Newton approximations to estimate the second
   *  derivatives and doesn't not neet to implement this method. */
  virtual bool eval_lagrangian_hessian(
      int num_variables, const double* primal_variables,
      bool new_primal_variables, double objective_scaling_factor,
      int num_constraints, const double* constraint_multipliers,
      bool new_constraint_multipliers, int num_nonzeros_hessian,
      int* row_indices, int* column_indices, double* nonzero_values) = 0;
  //@}

  /** @name Solution Methods */
  //@{
  /** This method is called when the algorithm is complete so the TNLP can
   * store/write the solution */
  virtual void finalize_solution(
      SqpSolverExitStatus status, int num_variables,
      const double* primal_solution, const double* bound_multipliers,
      const ActivityStatus* bound_activity_status, int num_constraints,
      const double* constraint_values, const double* constraint_multipliers,
      const ActivityStatus* constraint_activity_status, double objective_value,
      std::shared_ptr<const Statistics> stats) = 0;
  //@}

  /** Set the initial working set. */
  virtual bool set_initial_working_sets(
      int num_variables, const ActivityStatus* bounds_working_set,
      int num_constraints, const ActivityStatus* constraints_working_set)
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

    // Store the size of the problem so that we can check later whether the
    // requested working set is consistent.
    num_init_variables_ = num_variables;
    num_init_constraints_ = num_constraints;

    return true;
  }

  /** This method is used to determine if an initial working set is available.
   * If it returns true, the SQP algorithm will ask for the intiial working set.
   */
  bool use_initial_working_set() const
  {
    return initial_bounds_working_set_;
  }

  /** Set the initial working set. */
  bool get_initial_working_sets(int num_variables,
                                ActivityStatus* bounds_working_set,
                                int num_constraints,
                                ActivityStatus* constraints_working_set)
  {
    if (!initial_bounds_working_set_ || !initial_constraints_working_set_ ||
        num_variables != num_init_variables_ ||
        num_constraints != num_init_constraints_) {
      // The initial working sets should have been set earlier and the
      // dimensions must match
      return false;
    }

    for (int i = 0; i < num_variables; ++i) {
      bounds_working_set[i] = initial_bounds_working_set_[i];
    }
    for (int i = 0; i < num_constraints; ++i) {
      constraints_working_set[i] = initial_constraints_working_set_[i];
    }
    return true;
  }

  /** Get the starting point either from the overloaded get_starting_point
   * method or
   *  by the overwriting starting point that has been previously set. */
  bool get_correct_starting_point(int num_variables, bool init_primal_variables,
                                  double* initial_primal_variables,
                                  bool init_bound_multipliers,
                                  double* initial_bound_multipliers,
                                  int num_constraints,
                                  bool init_constraint_multipliers,
                                  double* initial_constraint_multipliers)
  {
    if (overwrite_starting_point_) {
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

    } else {
      return get_starting_point(
          num_variables, init_primal_variables, initial_primal_variables,
          init_bound_multipliers, initial_bound_multipliers, num_constraints,
          init_constraint_multipliers, initial_constraint_multipliers);
    }
  }

  /** Return true if an overwriting starting point has been set. */
  bool use_overwriting_starting_point() const
  {
    return overwrite_starting_point_;
  }

  /** Disable the overwriting starting point. */
  void disable_overwriting_starting_point()
  {
    overwrite_starting_point_ = false;

    delete[] initial_primal_variables_;
    initial_primal_variables_ = nullptr;
    delete[] initial_bound_multipliers_;
    initial_bound_multipliers_ = nullptr;
    delete[] initial_constraint_multipliers_;
    initial_constraint_multipliers_ = nullptr;
  }

  /** Retrieve the overwriting starting point. */
  void set_overwriting_starting_point(
      int num_variables, const double* initial_primal_variables,
      const double* initial_bound_multipliers, int num_constraints,
      const double* initial_constraint_multipliers)
  {
    disable_overwriting_starting_point();

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

    overwrite_starting_point_ = true;
  }

private:
  /** Copy Constructor */
  SqpTNlp(const SqpTNlp&);

  /** Overloaded Equals Operator */
  void operator=(const SqpTNlp&);
  //@}

  /** @name Initial working sets */
  //@{
  /** Working set for the variable bounds. */
  ActivityStatus* initial_bounds_working_set_;
  /** Working set for the contraint bounds. */
  ActivityStatus* initial_constraints_working_set_;
  /** Number of variables when the initial working set was given. */
  int num_init_variables_;
  /** Number of constraints when the initial working set was given. */
  int num_init_constraints_;
  //@}

  /** @name Starting point. */
  //@{
  /** Flag indicating whether a starting point has been given earlier that is
   * overwriting the starting point
   *  obtained from get_starting_point. */
  bool overwrite_starting_point_;
  /** Initial primal starting point. */
  double* initial_primal_variables_;
  /** Initial bound multipliers. */
  double* initial_bound_multipliers_;
  /** Initial constraint multipliers. */
  double* initial_constraint_multipliers_;
  //@}
};
}

#endif // SQPHOTSTART_SqpTnlp_HPP
