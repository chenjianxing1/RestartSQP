/* Copyright (C) 2019
* All Rights Reserved.
*
* Authors: Xinyi Luo
* Date:2019-07
*/
#include "sqphot/SqpTNlp.hpp"

using namespace std;
using namespace Ipopt;

namespace SQPhotstart {

/** Default constructor*/
SqpTNlp::SqpTNlp(SmartPtr<TNLP> ipopt_tnlp, string nlp_name)
 : ipopt_tnlp_(ipopt_tnlp)
 , nlp_name_(nlp_name)
{
  TNLP::IndexStyleEnum index_style;
  ipopt_tnlp_->get_nlp_info(num_variables_, num_constraints_,
                            num_nonzeros_jacobian_, num_nonzeros_hessian_,
                            index_style);
  assert(index_style == TNLP::FORTRAN_STYLE);
}

/** Default constructor*/
SqpTNlp::~SqpTNlp()
{
}

shared_ptr<const SqpNlpSizeInfo> SqpTNlp::get_problem_sizes()
{
  shared_ptr<SqpNlpSizeInfo> retval = make_shared<SqpNlpSizeInfo>(
      num_variables_, num_constraints_, num_nonzeros_jacobian_,
      num_nonzeros_hessian_);
  return retval;
}

/**
 *@brief get the bounds information from the NLP object
 */
bool SqpTNlp::get_bounds_info(shared_ptr<Vector> x_l, shared_ptr<Vector> x_u,
                              shared_ptr<Vector> c_l, shared_ptr<Vector> c_u)
{

  ipopt_tnlp_->get_bounds_info(num_variables_, x_l->get_values(),
                               x_u->get_values(), num_constraints_,
                               c_l->get_values(), c_u->get_values());
  return true;
}

/*
 * @brief Get the starting point from the NLP object.
 */
bool SqpTNlp::get_starting_point(shared_ptr<Vector> x_0,
                                 shared_ptr<Vector> lambda_0)
{
  ipopt_tnlp_->get_starting_point(num_variables_, true, x_0->get_values(),
                                  false, NULL, NULL, num_constraints_, true,
                                  lambda_0->get_values());

  return true;
}

/**
 *@brief Evaluate the objective value
 */
bool SqpTNlp::eval_f(shared_ptr<const Vector> x, double& obj_value)
{
  ipopt_tnlp_->eval_f(num_variables_, x->get_values(), true, obj_value);
  return true;
}

/**
 * @brief Evaluate the constraints at point x
 *
 */
bool SqpTNlp::eval_constraints(shared_ptr<const Vector> x,
                               shared_ptr<Vector> constraints)
{
  ipopt_tnlp_->eval_g(num_variables_, x->get_values(), true, num_constraints_,
                      constraints->get_values());
  return true;
}

/**
 *@brief Evaluate gradient at point x
 */
bool SqpTNlp::eval_gradient(shared_ptr<const Vector> x,
                            shared_ptr<Vector> gradient)
{
  ipopt_tnlp_->eval_grad_f(num_variables_, x->get_values(), true,
                           gradient->get_values());
  return true;
}

/**
 * @brief Get the matrix structure of the Jacobian
 * Always call this before the first time using @Eval_Jacobian
 */

bool SqpTNlp::get_jacobian_structure(shared_ptr<const Vector> x,
                                     shared_ptr<SpTripletMat> Jacobian)
{
  ipopt_tnlp_->eval_jac_g(num_variables_, x->get_values(), true,
                          num_constraints_, num_nonzeros_jacobian_,
                          Jacobian->get_nonconst_row_indices(),
                          Jacobian->get_nonconst_column_indices(), NULL);
  return true;
}

/**
 *@brief Evaluate Jacobian at point x
 */

bool SqpTNlp::eval_jacobian(shared_ptr<const Vector> x,
                            shared_ptr<SpTripletMat> Jacobian)
{
  ipopt_tnlp_->eval_jac_g(num_variables_, x->get_values(), true,
                          num_constraints_, num_nonzeros_jacobian_, NULL, NULL,
                          Jacobian->get_nonconst_values());
  return true;
}

/**
 * @brief Get the structure of the Hessian
 * Always call this before the first time using @Eval_Hessian
 */
bool SqpTNlp::get_hessian_structure(shared_ptr<const Vector> x,
                                    shared_ptr<const Vector> lambda,
                                    shared_ptr<SpTripletMat> Hessian)
{
  // We need to change the sign of the multipliers
  // TOO: Make consistent
  int dim_lambda = lambda->get_dim();
  double scale_factor = -1.;
  shared_ptr<Vector> negative_lambda = make_shared<Vector>(dim_lambda);
  negative_lambda->copy_vector(lambda, scale_factor);

  ipopt_tnlp_->eval_h(num_variables_, x->get_values(), true, 1.0,
                      num_constraints_, negative_lambda->get_values(), true,
                      num_nonzeros_hessian_,
                      Hessian->get_nonconst_row_indices(),
                      Hessian->get_nonconst_column_indices(), NULL);

  return true;
}

/**
 *@brief Evaluate Hessian of Lagragian function at  (x, lambda)
 */
bool SqpTNlp::eval_hessian(shared_ptr<const Vector> x,
                           shared_ptr<const Vector> lambda,
                           shared_ptr<SpTripletMat> Hessian)
{
  ipopt_tnlp_->eval_h(num_variables_, x->get_values(), true, 1,
                      num_constraints_, lambda->get_values(), true,
                      num_nonzeros_hessian_, NULL, NULL,
                      Hessian->get_nonconst_values());

  return true;
}
}
