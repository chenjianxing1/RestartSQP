/* Copyright (C) 2019
* All Rights Reserved.
*
* Authors: Xinyi Luo
* Date:2019-07
*/
#include "sqphot/SQPTNLP.hpp"

using namespace std;

namespace SQPhotstart {

/** Default constructor*/
SQPTNLP::SQPTNLP(Ipopt::SmartPtr<Ipopt::TNLP> nlp)
{
  nlp_ = nlp;
  Ipopt::TNLP::IndexStyleEnum index_style;
  nlp_->get_nlp_info(nlp_info_.nVar, nlp_info_.nCon, nlp_info_.nnz_jac_g,
                     nlp_info_.nnz_h_lag, index_style);
  assert(index_style == Ipopt::TNLP::FORTRAN_STYLE);
}

/** Default constructor*/
SQPTNLP::~SQPTNLP()
{
}

/**
 *@brief get the bounds information from the NLP object
 */
bool SQPTNLP::get_bounds_info(shared_ptr<Vector> x_l, shared_ptr<Vector> x_u,
                              shared_ptr<Vector> c_l, shared_ptr<Vector> c_u)
{

  nlp_->get_bounds_info(nlp_info_.nVar, x_l->get_values(), x_u->get_values(),
                        nlp_info_.nCon, c_l->get_values(), c_u->get_values());
  return true;
}

/*
 * @brief Get the starting point from the NLP object.
 */
bool SQPTNLP::get_starting_point(shared_ptr<Vector> x_0,
                                 shared_ptr<Vector> lambda_0)
{
  nlp_->get_starting_point(nlp_info_.nVar, true, x_0->get_values(), false, NULL,
                           NULL, nlp_info_.nCon, true, lambda_0->get_values());

  return true;
}

/**
 *@brief Evaluate the objective value
 */
bool SQPTNLP::eval_f(shared_ptr<const Vector> x, double& obj_value)
{
  nlp_->eval_f(nlp_info_.nVar, x->get_values(), true, obj_value);
  return true;
}

/**
 * @brief Evaluate the constraints at point x
 *
 */
bool SQPTNLP::eval_constraints(shared_ptr<const Vector> x,
                               shared_ptr<Vector> constraints)
{
  nlp_->eval_g(nlp_info_.nVar, x->get_values(), true, nlp_info_.nCon,
               constraints->get_values());
  return true;
}

/**
 *@brief Evaluate gradient at point x
 */
bool SQPTNLP::eval_gradient(shared_ptr<const Vector> x,
                            shared_ptr<Vector> gradient)
{
  nlp_->eval_grad_f(nlp_info_.nVar, x->get_values(), true,
                    gradient->get_values());
  return true;
}

/**
 * @brief Get the matrix structure of the Jacobian
 * Always call this before the first time using @Eval_Jacobian
 */

bool SQPTNLP::get_jacobian_structure(shared_ptr<const Vector> x,
                                     shared_ptr<SpTripletMat> Jacobian)
{
  nlp_->eval_jac_g(nlp_info_.nVar, x->get_values(), true, nlp_info_.nCon,
                   nlp_info_.nnz_jac_g, Jacobian->get_nonconst_row_indices(),
                   Jacobian->get_nonconst_column_indices(), NULL);
  return true;
}

/**
 *@brief Evaluate Jacobian at point x
 */

bool SQPTNLP::eval_jacobian(shared_ptr<const Vector> x,
                            shared_ptr<SpTripletMat> Jacobian)
{
  nlp_->eval_jac_g(nlp_info_.nVar, x->get_values(), true, nlp_info_.nCon,
                   nlp_info_.nnz_jac_g, NULL, NULL,
                   Jacobian->get_nonconst_values());
  return true;
}

/**
 * @brief Get the structure of the Hessian
 * Always call this before the first time using @Eval_Hessian
 */
bool SQPTNLP::get_hessian_structure(shared_ptr<const Vector> x,
                                    shared_ptr<const Vector> lambda,
                                    shared_ptr<SpTripletMat> Hessian)
{

  nlp_->eval_h(nlp_info_.nVar, x->get_values(), true, 1.0, nlp_info_.nCon,
               lambda->get_values(), true, nlp_info_.nnz_h_lag,
               Hessian->get_nonconst_row_indices(),
               Hessian->get_nonconst_column_indices(), NULL);

  return true;
}

/**
 *@brief Evaluate Hessian of Lagragian function at  (x, lambda)
 */
bool SQPTNLP::eval_hessian(shared_ptr<const Vector> x,
                           shared_ptr<const Vector> lambda,
                           shared_ptr<SpTripletMat> Hessian)
{
  nlp_->eval_h(nlp_info_.nVar, x->get_values(), true, 1, nlp_info_.nCon,
               lambda->get_values(), true, nlp_info_.nnz_h_lag, NULL, NULL,
               Hessian->get_nonconst_values());

  return true;
}
}
