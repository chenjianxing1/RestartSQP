/* Copyright (C) 2019
 * All Rights Reserved.
 *
 * Authors: Xinyi Luo
 * Date:2019-06
 */
#ifndef SQPHOTSTART_SQPTNLP_HPP
#define SQPHOTSTART_SQPTNLP_HPP

#include "IpTNLP.hpp"
#include "sqphot/SpTripletMat.hpp"
#include "sqphot/Utils.hpp"
#include "sqphot/Vector.hpp"
#include <memory>

namespace SQPhotstart {
/**
 * This is part of SQPhotstart
 *
 * This class enables user to read data from NLP class object with more friendly
 * names and the use of Matrix and Vector objects for data.
 *
 */
class SQPTNLP
{

public:
  /** @brief constructor that copies nlp to _nlp as a local data reader*/
  SQPTNLP(Ipopt::SmartPtr<Ipopt::TNLP> nlp);

  /** Default destructor*/
  virtual ~SQPTNLP();

  /**
   *@brief get the bounds information from the NLP object
   */
  virtual bool get_bounds_info(std::shared_ptr<Vector> x_l,
                               std::shared_ptr<Vector> x_u,
                               std::shared_ptr<Vector> c_l,
                               std::shared_ptr<Vector> c_u);

  /*
   * @brief Get the starting point from the NLP object.
   * TODO: add options_ to enable user to choose if to use default input or not
   */
  virtual bool get_starting_point(std::shared_ptr<Vector> x_0,
                                  std::shared_ptr<Vector> lambda_0);

  /**
   *@brief Evaluate the objective value
   */
  virtual bool eval_f(std::shared_ptr<const Vector> x, double& obj_value);

  /**
   * @brief Evaluate the constraints at point x
   *
   */
  virtual bool eval_constraints(std::shared_ptr<const Vector> x,
                                std::shared_ptr<Vector> constraints);

  /**
   *@brief Evaluate gradient at point x
   */
  virtual bool eval_gradient(std::shared_ptr<const Vector> x,
                             std::shared_ptr<Vector> gradient);

  /**
   * @brief Get the matrix structure of the Jacobian
   * Always call this before the first time using @Eval_Jacobian
   */
  virtual bool get_jacobian_structure(std::shared_ptr<const Vector> x,
                                      std::shared_ptr<SpTripletMat> Jacobian);

  /**
   *@brief Evaluate Jacobian at point x
   */
  virtual bool eval_jacobian(std::shared_ptr<const Vector> x,
                             std::shared_ptr<SpTripletMat> Jacobian);

  /**
   * @brief Get the structure of the Hessian
   * Always call this before the first time using @Eval_Hessian
   */
  virtual bool get_hessian_structure(std::shared_ptr<const Vector> x,
                                     std::shared_ptr<const Vector> lambda,
                                     std::shared_ptr<SpTripletMat> Hessian);

  /**
   *@brief Evaluate Hessian of Lagragian function at  (x, lambda)
   */
  virtual bool eval_hessian(std::shared_ptr<const Vector> x,
                            std::shared_ptr<const Vector> lambda,
                            std::shared_ptr<SpTripletMat> Hessian);

public:
  NLPInfo nlp_info_; /**< the struct record the number of variables, number of
                             constraints, number of nonzeoro entry of Hessian
                        and that of Jacobian
                             Please check Types.hpp for details*/
  Ipopt::SmartPtr<Ipopt::TNLP> nlp_; /**< a local nlp reader */

private:
  /** Default constructor*/
  SQPTNLP();

  /** Copy Constructor */
  SQPTNLP(const SQPTNLP&);

  /** Overloaded Equals Operator */
  void operator=(const SQPTNLP&);
  //@}
};
}

#endif // CPP_SQPTNLP_HPP
