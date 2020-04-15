/* Copyright (C) 2019
 * All Rights Reserved.
 *
 * Authors: Xinyi Luo
 * Date:2019-06
 */
#ifndef SQPHOTSTART_SQPNLP_HPP
#define SQPHOTSTART_SQPNLP_HPP

#include "restartsqp/SparseTripletMatrix.hpp"
#include "restartsqp/SqpTNlp.hpp"
#include "restartsqp/Vector.hpp"
#include <memory>

namespace RestartSqp {

/** Simple class for storing the size information about an NLP. */
class SqpNlpSizeInfo
{
public:
  /** constructor/destructor */
  //@{

  /**Default constructor*/
  SqpNlpSizeInfo(int num_variables, int num_constraints,
                 int num_nonzeros_jacobian, int num_nonzeros_hessian)
   : num_variables_(num_variables)
   , num_constraints_(num_constraints)
   , num_nonzeros_jacobian_(num_nonzeros_jacobian)
   , num_nonzeros_hessian_(num_nonzeros_hessian)
  {
  }

  /** Destructor */
  ~SqpNlpSizeInfo()
  {
  }

  /** @name Accessor Methods. */
  //@{
  /** Get number of variables. */
  int get_num_variables() const
  {
    return num_variables_;
  }
  /** Get number of constraints. */
  int get_num_constraints() const
  {
    return num_constraints_;
  }
  /** Get number of nonzeros in constraint Jacobian. */
  int get_num_nonzeros_jacobian() const
  {
    return num_nonzeros_jacobian_;
  }
  /** Get number of nonzeros in Lagrangian Hessian. */
  int get_num_nonzeros_hessian() const
  {
    return num_nonzeros_hessian_;
  }

private:
  /** Default constructor*/
  SqpNlpSizeInfo();

  /** Copy Constructor */
  SqpNlpSizeInfo(const SqpNlpSizeInfo&);

  /** Overloaded Equals Operator */
  void operator=(const SqpNlpSizeInfo&);

  /** @name Data members. */
  //@{
  /** Number of variables. */
  const int num_variables_;

  /** Number of constraints. */
  const int num_constraints_;

  /** Number of nonzeros in constraint Jacobian. */
  const int num_nonzeros_jacobian_;

  /** Number of nonzeros in Lagrangian Hessian. */
  const int num_nonzeros_hessian_;
  //@}
};

/**
 * This is part of SQPhotstart
 *
 * This class enables user to read data from NLP class object with more friendly
 * names and the use of Matrix and Vector objects for data.
 *
 * IMPORTANT: The Lagrangian function here is defined as L(x,l) = f(x) - sum_i
 * l_ic_j(c)
 *
 *            THIS IS DIFFERENT FROM IPOPT's DEFINITION!!
 */
class SqpNlp
{

public:
  /** @brief Constructor that an instance of Ipopt's TNLP */
  SqpNlp(std::shared_ptr<SqpTNlp> sqp_tnlp);

  /** Destructor*/
  ~SqpNlp();

  /** Get problem name. */
  const std::string& get_nlp_name() const
  {
    return nlp_name_;
  }

  /**
   *@brief Get problem size information
   */
  std::shared_ptr<const SqpNlpSizeInfo> get_problem_sizes();

  /**
   *@brief get the bounds information from the NLP object
   */
  bool get_bounds_info(std::shared_ptr<Vector> x_l, std::shared_ptr<Vector> x_u,
                       std::shared_ptr<Vector> c_l,
                       std::shared_ptr<Vector> c_u);

  /** Method indicating whether an initial working set is available.
   *
   *  Returns true if that is the case.
   */
  bool use_initial_working_set();

  /** Get the initial working sets. */
  // TODO: Create a new WorkingSet object
  bool get_initial_working_sets(int num_variables,
                                ActivityStatus* bounds_working_set,
                                int num_constraints,
                                ActivityStatus* constraints_working_set);
  /*
   * @brief Get the starting point from the NLP object.
   * TODO: add options_ to enable user to choose if to use default input or not
   */
  bool get_starting_point(std::shared_ptr<Vector> primal_point,
                          std::shared_ptr<Vector> bound_multipliers,
                          std::shared_ptr<Vector> constraint_multipliers);

  /**
   *@brief Evaluate the objective value
   */
  bool eval_f(std::shared_ptr<const Vector> x, double& obj_value);

  /**
   * @brief Evaluate the constraints at point x
   *
   */
  bool eval_constraints(std::shared_ptr<const Vector> x,
                        std::shared_ptr<Vector> constraints);

  /**
   *@brief Evaluate gradient at point x
   */
  bool eval_gradient(std::shared_ptr<const Vector> x,
                     std::shared_ptr<Vector> gradient);

  /**
   * @brief Get the matrix structure of the Jacobian
   * Always call this before the first time using @Eval_Jacobian
   */
  bool get_jacobian_structure(std::shared_ptr<const Vector> x,
                              std::shared_ptr<SparseTripletMatrix> Jacobian);

  /**
   *@brief Evaluate Jacobian at point x
   */
  bool eval_jacobian(std::shared_ptr<const Vector> x,
                     std::shared_ptr<SparseTripletMatrix> Jacobian);

  /**
   * @brief Get the structure of the Hessian
   * Always call this before the first time using @Eval_Hessian
   */
  bool get_hessian_structure(std::shared_ptr<const Vector> x,
                             std::shared_ptr<const Vector> lambda,
                             std::shared_ptr<SparseTripletMatrix> Hessian);

  /**
   *@brief Evaluate Hessian of Lagragian function at  (x, lambda)
   */
  bool eval_hessian(std::shared_ptr<const Vector> x,
                    std::shared_ptr<const Vector> lambda,
                    double objective_scaling_factor,
                    std::shared_ptr<SparseTripletMatrix> Hessian);

  /**
   * @brief Return the results of the optimization run to the user.
   */
  bool finalize_solution(SqpSolverExitStatus status,
                         std::shared_ptr<const Vector> primal_solution,
                         std::shared_ptr<const Vector> bound_multipliers,
                         const ActivityStatus* bound_activity_status,
                         std::shared_ptr<const Vector> constraint_values,
                         std::shared_ptr<const Vector> constraint_multipliers,
                         const ActivityStatus* constraint_activity_status,
                         double objective_value,
                         std::shared_ptr<const Statistics> stats);

private:
  /** @name Hide unused default methods. */
  //@{
  /** Default constructor*/
  SqpNlp();

  /** Copy Constructor */
  SqpNlp(const SqpNlp&);

  /** Overloaded Equals Operator */
  void operator=(const SqpNlp&);
  //@}

  /** Ipopt's TNLP object that will be called for all evaluations. */
  std::shared_ptr<SqpTNlp> sqp_tnlp_;

  /** Name of the NLP. */
  std::string nlp_name_;

  /** Problem dimensioms. */
  //@{
  /** Number of variables. */
  int num_variables_;

  /** Number of constraints. */
  int num_constraints_;

  /** Number of nonzeros in constraint Jacobian. */
  int num_nonzeros_jacobian_;

  /** Number of nonzeros in Lagrangian Hessian. */
  int num_nonzeros_hessian_;

  //@}
};
}

#endif // SQPHOTSTART_SqpNlp_HPP
