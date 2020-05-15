/* Copyright (C) 2019
 * All Rights Reserved.
 *
 * Authors: Xinyi Luo
 * Date:2019-06
 */
#ifndef SQPHOTSTART_SQPTNLP_HPP
#define SQPHOTSTART_SQPTNLP_HPP

#include "restartsqp/Statistics.hpp"
#include "restartsqp/Types.hpp"
#include <memory>
#include <cassert>

namespace RestartSqp {

/**
 * This is part of SQPhotstart
 *
 * This class is the base class for the classes that a user has to provide.
 */
class SqpTNlp
{

public:
  /** @brief Default Constructor */
  SqpTNlp()
  {
  }

  /** Default destructor*/
  virtual ~SqpTNlp()
  {
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

  /** This method is used to determine if an initial working set is available.
   * If it returns true, the SQP algorithm will ask for the intiial working set.
   */
  virtual bool use_initial_working_set() = 0;

  /** Set the initial working set. */
  virtual bool get_initial_working_sets(int num_variables,
                                ActivityStatus* bounds_working_set,
                                int num_constraints,
                                ActivityStatus* constraints_working_set)
  {
    assert("get_initial_working_sets not implemented, but use_initial_working_set returned true.");
    return false;
  }

private:
  /** Copy Constructor */
  SqpTNlp(const SqpTNlp&);

  /** Overloaded Equals Operator */
  void operator=(const SqpTNlp&);
  //@}
};
}

#endif // SQPHOTSTART_SqpTnlp_HPP
