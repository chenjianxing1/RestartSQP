/* Copyright (C) 2019
 * All Rights Reserved.
 *
 * Authors: Xinyi Luo
 * Date:2019-06
 */
#ifndef SQPHOTSTART_SQPTNLP_HPP
#define SQPHOTSTART_SQPTNLP_HPP

#include "IpTNLP.hpp"
#include "sqphot/SqpNlpBase.hpp"

namespace SQPhotstart {
/**
 * This is part of SQPhotstart
 *
 * This class enables user to read data from NLP class object with more friendly
 * names and the use of Matrix and Vector objects for data.
 *
 */
class SqpTNlp : public SqpNlpBase
{

public:
  /** @brief Constructor that an instance of Ipopt's TNLP */
  SqpTNlp(Ipopt::SmartPtr<Ipopt::TNLP> ipopt_tnlp);

  /** Default destructor*/
  ~SqpTNlp();

  /**
   *@brief Get problem size information
   */
  std::shared_ptr<const SqpNlpSizeInfo> get_problem_sizes() override;

  /**
   *@brief get the bounds information from the NLP object
   */
  bool get_bounds_info(std::shared_ptr<Vector> x_l, std::shared_ptr<Vector> x_u,
                       std::shared_ptr<Vector> c_l,
                       std::shared_ptr<Vector> c_u) override;

  /*
   * @brief Get the starting point from the NLP object.
   * TODO: add options_ to enable user to choose if to use default input or not
   */
  bool get_starting_point(std::shared_ptr<Vector> x_0,
                          std::shared_ptr<Vector> lambda_0) override;

  /**
   *@brief Evaluate the objective value
   */
  bool eval_f(std::shared_ptr<const Vector> x, double& obj_value) override;

  /**
   * @brief Evaluate the constraints at point x
   *
   */
  bool eval_constraints(std::shared_ptr<const Vector> x,
                        std::shared_ptr<Vector> constraints) override;

  /**
   *@brief Evaluate gradient at point x
   */
  bool eval_gradient(std::shared_ptr<const Vector> x,
                     std::shared_ptr<Vector> gradient) override;

  /**
   * @brief Get the matrix structure of the Jacobian
   * Always call this before the first time using @Eval_Jacobian
   */
  bool get_jacobian_structure(std::shared_ptr<const Vector> x,
                              std::shared_ptr<SpTripletMat> Jacobian) override;

  /**
   *@brief Evaluate Jacobian at point x
   */
  bool eval_jacobian(std::shared_ptr<const Vector> x,
                     std::shared_ptr<SpTripletMat> Jacobian) override;

  /**
   * @brief Get the structure of the Hessian
   * Always call this before the first time using @Eval_Hessian
   */
  bool get_hessian_structure(std::shared_ptr<const Vector> x,
                             std::shared_ptr<const Vector> lambda,
                             std::shared_ptr<SpTripletMat> Hessian) override;

  /**
   *@brief Evaluate Hessian of Lagragian function at  (x, lambda)
   */
  bool eval_hessian(std::shared_ptr<const Vector> x,
                    std::shared_ptr<const Vector> lambda,
                    std::shared_ptr<SpTripletMat> Hessian) override;

private:
  /** @name Hide unused default methods. */
  //@{
  /** Default constructor*/
  SqpTNlp();

  /** Copy Constructor */
  SqpTNlp(const SqpTNlp&);

  /** Overloaded Equals Operator */
  void operator=(const SqpTNlp&);
  //@}

  /** Ipopt's TNLP object that will be called for all evaluations. */
  Ipopt::SmartPtr<Ipopt::TNLP> ipopt_tnlp_;

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

#endif // SQPHOTSTART_SQPTNLP_HPP
