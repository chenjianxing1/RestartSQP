/* Copyright (C) 2019
 * All Rights Reserved.
 *
 * Authors: Xinyi Luo
 * Date:2019-06
 */
#ifndef SQPHOTSTART_SQPNLPBASE_HPP
#define SQPHOTSTART_SQPNLPBASE_HPP

#include "sqphot/SparseTripletMatrix.hpp"
#include "sqphot/Statistics.hpp"
#include "sqphot/Vector.hpp"
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
 * This class is the base class that the optimization code directly interacts
 * with.
 *
 */
class SqpNlpBase
{

public:
  /** @brief Default Constructor */
  SqpNlpBase()
  {
  }

  /** Default destructor*/
  virtual ~SqpNlpBase()
  {
  }

  /** @name Name of the NLP. */
  virtual const std::string& get_nlp_name() const = 0;

  /**
   *@brief get the information about sizes
   */
  virtual std::shared_ptr<const SqpNlpSizeInfo> get_problem_sizes() = 0;

  /**
   *@brief get the bounds information from the NLP object
   */
  virtual bool get_bounds_info(std::shared_ptr<Vector> x_l,
                               std::shared_ptr<Vector> x_u,
                               std::shared_ptr<Vector> c_l,
                               std::shared_ptr<Vector> c_u) = 0;

  /*
   * @brief Get the starting point from the NLP object.
   * TODO: add options_ to enable user to choose if to use default input or not
   */
  virtual bool
  get_starting_point(std::shared_ptr<Vector> primal_point,
                     std::shared_ptr<Vector> bound_multipliers,
                     std::shared_ptr<Vector> constraint_multipliers) = 0;

  /**
   *@brief Evaluate the objective value
   */
  virtual bool eval_f(std::shared_ptr<const Vector> x, double& obj_value) = 0;

  /**
   * @brief Evaluate the constraints at point x
   *
   */
  virtual bool eval_constraints(std::shared_ptr<const Vector> x,
                                std::shared_ptr<Vector> constraints) = 0;

  /**
   *@brief Evaluate gradient at point x
   */
  virtual bool eval_gradient(std::shared_ptr<const Vector> x,
                             std::shared_ptr<Vector> gradient) = 0;

  /**
   * @brief Get the matrix structure of the Jacobian
   * Always call this before the first time using @Eval_Jacobian
   */
  virtual bool
  get_jacobian_structure(std::shared_ptr<const Vector> x,
                         std::shared_ptr<SparseTripletMatrix> jacobian) = 0;

  /**
   *@brief Evaluate Jacobian at point x
   */
  virtual bool eval_jacobian(std::shared_ptr<const Vector> x,
                             std::shared_ptr<SparseTripletMatrix> jacobian) = 0;

  /**
   * @brief Get the structure of the Hessian
   * Always call this before the first time using @Eval_Hessian
   */
  virtual bool
  get_hessian_structure(std::shared_ptr<const Vector> x,
                        std::shared_ptr<const Vector> lambda,
                        std::shared_ptr<SparseTripletMatrix> hessian) = 0;

  /**
   *@brief Evaluate Hessian of Lagragian function at  (x, lambda)
   */
  virtual bool eval_hessian(std::shared_ptr<const Vector> x,
                            std::shared_ptr<const Vector> lambda,
                            std::shared_ptr<SparseTripletMatrix> hessian) = 0;

  /**
   * @brief Return the results of the optimization run to the user.
   */
  virtual bool finalize_solution(
      SqpSolverExitStatus status, std::shared_ptr<const Vector> primal_solution,
      std::shared_ptr<const Vector> bound_multipliers,
      const ActivityStatus* bound_activity_status,
      std::shared_ptr<const Vector> constraint_values,
      std::shared_ptr<const Vector> constraint_multipliers,
      const ActivityStatus* constraint_activity_status, double objective_value,
      std::shared_ptr<const Statistics>) = 0;

private:
  /** Copy Constructor */
  SqpNlpBase(const SqpNlpBase&);

  /** Overloaded Equals Operator */
  void operator=(const SqpNlpBase&);
  //@}
};
}

#endif // SQPHOTSTART_SQPNLPBASE_HPP
