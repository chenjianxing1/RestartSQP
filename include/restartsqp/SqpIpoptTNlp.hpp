/* Copyright (C) 2019
 * All Rights Reserved.
 *
 * Authors: Xinyi Luo
 * Date:2019-06
 */
#ifndef SQPHOTSTART_SQPIPOPTNLP_HPP
#define SQPHOTSTART_SQPIPOPTNLP_HPP

#include "IpTNLP.hpp"
#include "restartsqp/SqpTNlp.hpp"

namespace RestartSqp {
/**
 * This is part of SQPhotstart
 *
 * This makes an Ipopt::TNLP look like an SqpTNlp.
 *
 * IMPORTANT: The Lagrangian function here is defined as L(x,l) = f(x) - sum_i
 * l_ic_j(c)
 *
 *            THIS IS DIFFERENT FROM IPOPT's DEFINITION!!
 */
class SqpIpoptTNlp : public SqpTNlp
{

public:
  /** @brief Constructor that an instance of Ipopt's TNLP */
  SqpIpoptTNlp(Ipopt::SmartPtr<Ipopt::TNLP> ipopt_tnlp,
               const std::string& nlp_name = "Ipopt NLP");

  /** Destructor*/
  virtual ~SqpIpoptTNlp();

  /**
   *@brief Get problem size information
   */
  virtual bool get_nlp_info(int& num_variables, int& num_constraints,
                            int& num_nonzeros_jacobian,
                            int& num_nonzeros_hessian,
                            std::string& nlp_name) override;

  /**
   *@brief get the bounds information from the NLP object
   */
  virtual bool get_bounds_info(int num_variabes, double* variable_lower_bounds,
                               double* variable_upper_bounds,
                               int num_constraints,
                               double* constraint_lower_bounds,
                               double* constraint_upper_bounds) override;

  /*
   * @brief Get the starting point from the NLP object.
   * TODO: add options_ to enable user to choose if to use default input or not
   */
  virtual bool get_starting_point(int num_variables, bool init_primal_variables,
                                  double* primal_variables,
                                  bool init_bound_multipliers,
                                  double* bound_multipliers,
                                  int num_constraints,
                                  bool init_constraint_multipliers,
                                  double* constraint_multipliers) override;

  /**
   *@brief Evaluate the objective value
   */
  virtual bool eval_objective_value(int num_variables,
                                    const double* primal_variables,
                                    bool new_primal_variables,
                                    double& objective_value) override;

  /**
   *@brief Evaluate gradient at point x
   */
  virtual bool eval_objective_gradient(int num_variables,
                                       const double* primal_variables,
                                       bool new_primal_variables,
                                       double* objective_gradient) override;

  /**
   * @brief Evaluate the constraints at point x
   *
   */
  virtual bool eval_constraint_values(int num_variables,
                                      const double* primal_variables,
                                      bool new_primal_variables,
                                      int num_constraints,
                                      double* constraint_values) override;

  /**
   * @brief Get the matrix structure of the Jacobian
   * Always call this before the first time using @Eval_Jacobian
   */
  virtual bool eval_constraint_jacobian(
      int num_variables, const double* primal_variables,
      bool new_primal_variables, int num_constraints, int num_nonzeros_jacobian,
      int* row_indices, int* column_indices, double* nonzero_values) override;

  /**
   *@brief Evaluate Hessian of Lagragian function at  (x, lambda)
   */
  virtual bool eval_lagrangian_hessian(
      int num_variables, const double* primal_variables,
      bool new_primal_variables, double objective_scaling_factor,
      int num_constraints, const double* constraint_multipliers,
      bool new_constraint_multipliers, int num_nonzeros_hessian,
      int* row_indices, int* column_indices, double* nonzero_values) override;

  /**
   * @brief Return the results of the optimization run to the user.
   */
  virtual void finalize_solution(
      SqpSolverExitStatus status, int num_variables,
      const double* primal_solution, const double* bound_multipliers,
      const ActivityStatus* bound_activity_status, int num_constraints,
      const double* constraint_values, const double* constraint_multipliers,
      const ActivityStatus* constraint_activity_status, double objective_value,
      std::shared_ptr<const Statistics> stats) override;

  /** No initial working set is available from an Ipopt TNLP */
  virtual bool use_initial_working_set() override
  {
    return false;
  }

private:
  /** @name Hide unused default methods. */
  //@{
  /** Default constructor*/
  SqpIpoptTNlp();

  /** Copy Constructor */
  SqpIpoptTNlp(const SqpIpoptTNlp&);

  /** Overloaded Equals Operator */
  void operator=(const SqpIpoptTNlp&);
  //@}

protected:
  /** Ipopt's TNLP object that will be called for all evaluations. */
  Ipopt::SmartPtr<Ipopt::TNLP> ipopt_tnlp_;

  /** Numbering style for the matrices from the Ipopt TNLP. */
  Ipopt::TNLP::IndexStyleEnum index_style_;

  /** Name of the NLP that is solved. */
  std::string nlp_name_;
};
} // namespace RestartSqp

#endif // SQPHOTSTART_SqpIpoptNlp_HPP
