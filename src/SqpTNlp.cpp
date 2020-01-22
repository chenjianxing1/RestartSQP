/* Copyright (C) 2019
* All Rights Reserved.
*
* Authors: Xinyi Luo
* Date:2019-07
*/
#include "sqphot/SqpTNlp.hpp"

using namespace std;
using namespace Ipopt;

namespace RestartSqp {

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
  // Large value for bounds that should be interpreted as infinity
  const double NLP_INF = 1e18;

  double* x_L = x_l->get_non_const_values();
  double* x_U = x_u->get_non_const_values();
  double* c_L = c_l->get_non_const_values();
  double* c_U = c_u->get_non_const_values();

  bool retval;
  retval = ipopt_tnlp_->get_bounds_info(num_variables_, x_L, x_U, num_constraints_, c_L, c_U);
  if (!retval) {
    return false;
  }
  // TODO: Figure out if this is necessary:
  for (int i=0; i<num_variables_; ++i) {
    if (x_L[i] <= -NLP_INF) {
      x_L[i] = -INF;
    }
    if (x_U[i] >= NLP_INF) {
      x_U[i] = INF;
    }
  }
  for (int i=0; i<num_constraints_; ++i) {
    if (c_L[i] <= -NLP_INF) {
      c_L[i] = -INF;
    }
    if (c_U[i] >= NLP_INF) {
      c_U[i] = INF;
    }
  }
  return true;
}

/*
 * @brief Get the starting point from the NLP object.
 */
bool SqpTNlp::get_starting_point(shared_ptr<Vector> primal_point,
                                 shared_ptr<Vector> bound_multipliers,
                                 shared_ptr<Vector> constraint_multipliers)
{
  //
  // TODO: There is a problem with getting bound multipliers in the Ipopt AMPL
  // interface (somewhat related to the suffix handler).  This might be an Ipopt
  // bug; need to look into this?
  const bool bug_fixed = false;

  bool retval;

  if (bug_fixed) {
    double* z_L = new double[num_variables_];
    double* z_U = new double[num_variables_];
    retval = ipopt_tnlp_->get_starting_point(
        num_variables_, true, primal_point->get_non_const_values(), true, z_L,
        z_U, num_constraints_, true,
        constraint_multipliers->get_non_const_values());
    // We need to collapse the bound multipliers into a single vector
    for (int i = 0; i < num_variables_; ++i) {
      z_L[i] -= z_U[i];
    }
    bound_multipliers->copy_values(z_L);

    delete[] z_L;
    delete[] z_U;
  } else {
    retval = ipopt_tnlp_->get_starting_point(
        num_variables_, true, primal_point->get_non_const_values(), false, NULL,
        NULL, num_constraints_, true,
        constraint_multipliers->get_non_const_values());
    bound_multipliers->set_to_zero();
  }

  return retval;
}

/**
 *@brief Evaluate the objective value
 */
bool SqpTNlp::eval_f(shared_ptr<const Vector> x, double& obj_value)
{
  return ipopt_tnlp_->eval_f(num_variables_, x->get_values(), true, obj_value);
}

/**
 * @brief Evaluate the constraints at point x
 *
 */
bool SqpTNlp::eval_constraints(shared_ptr<const Vector> x,
                               shared_ptr<Vector> constraints)
{
  return ipopt_tnlp_->eval_g(num_variables_, x->get_values(), true, num_constraints_,
                      constraints->get_non_const_values());
}

/**
 *@brief Evaluate gradient at point x
 */
bool SqpTNlp::eval_gradient(shared_ptr<const Vector> x,
                            shared_ptr<Vector> gradient)
{
  return ipopt_tnlp_->eval_grad_f(num_variables_, x->get_values(), true,
                           gradient->get_non_const_values());
}

/**
 * @brief Get the matrix structure of the Jacobian
 * Always call this before the first time using @Eval_Jacobian
 */

bool SqpTNlp::get_jacobian_structure(shared_ptr<const Vector> x,
                                     shared_ptr<SparseTripletMatrix> Jacobian)
{
  return ipopt_tnlp_->eval_jac_g(num_variables_, x->get_values(), true,
                          num_constraints_, num_nonzeros_jacobian_,
                          Jacobian->get_nonconst_row_indices(),
                          Jacobian->get_nonconst_column_indices(), NULL);
}

/**
 *@brief Evaluate Jacobian at point x
 */

bool SqpTNlp::eval_jacobian(shared_ptr<const Vector> x,
                            shared_ptr<SparseTripletMatrix> Jacobian)
{
  return ipopt_tnlp_->eval_jac_g(num_variables_, x->get_values(), true,
                          num_constraints_, num_nonzeros_jacobian_, NULL, NULL,
                          Jacobian->get_nonconst_values());
}

/**
 * @brief Get the structure of the Hessian
 * Always call this before the first time using @Eval_Hessian
 */
bool SqpTNlp::get_hessian_structure(shared_ptr<const Vector> x,
                                    shared_ptr<const Vector> lambda,
                                    shared_ptr<SparseTripletMatrix> Hessian)
{
  // We need to change the sign of the multipliers
  // TOO: Make consistent
  int dim_lambda = lambda->get_dim();
  double scale_factor = -1.;
  shared_ptr<Vector> negative_lambda = make_shared<Vector>(dim_lambda);
  negative_lambda->copy_vector(lambda, scale_factor);

  return ipopt_tnlp_->eval_h(num_variables_, x->get_values(), true, 1.0,
                      num_constraints_, negative_lambda->get_values(), true,
                      num_nonzeros_hessian_,
                      Hessian->get_nonconst_row_indices(),
                      Hessian->get_nonconst_column_indices(), NULL);
}

/**
 *@brief Evaluate Hessian of Lagragian function at  (x, lambda)
 */
bool SqpTNlp::eval_hessian(shared_ptr<const Vector> x,
                           shared_ptr<const Vector> lambda,
                           shared_ptr<SparseTripletMatrix> Hessian)
{
  return ipopt_tnlp_->eval_h(num_variables_, x->get_values(), true, 1,
                      num_constraints_, lambda->get_values(), true,
                      num_nonzeros_hessian_, NULL, NULL,
                      Hessian->get_nonconst_values());
}
}
