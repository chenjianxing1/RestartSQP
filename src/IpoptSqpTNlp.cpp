/* Copyright (C) 2019
 * All Rights Reserved.
 *
 * Authors: Xinyi Luo
 * Date:2019-07
 */
#include "restartsqp/IpoptSqpTNlp.hpp"
#include <cassert>

using namespace std;
using namespace Ipopt;

namespace RestartSqp {

/** Default constructor*/
IpoptSqpTNlp::IpoptSqpTNlp(std::shared_ptr<SqpTNlp> sqp_tnlp)
 : sqp_tnlp_(sqp_tnlp)
 , x_sol_(nullptr)
 , z_L_sol_(nullptr)
 , z_U_sol_(nullptr)
 , lambda_sol_(nullptr)
 , g_sol_(nullptr)
{}

/** Default constructor*/
IpoptSqpTNlp::~IpoptSqpTNlp()
{
  delete[] x_sol_;
  delete[] z_L_sol_;
  delete[] z_U_sol_;
  delete[] lambda_sol_;
  delete[] g_sol_;
}

bool IpoptSqpTNlp::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                                Index& nnz_h_lag,
                                TNLP::IndexStyleEnum& index_style)
{
  string nlp_name = "IpoptTNLP";
  bool retval;
  retval = sqp_tnlp_->get_nlp_info(n, m, nnz_jac_g, nnz_h_lag, nlp_name);

  index_style = TNLP::C_STYLE;

  return retval;
}

/**
 *@brief get the bounds information from the NLP object
 */
bool IpoptSqpTNlp::get_bounds_info(Index n, Number* x_l, Number* x_u, Index m,
                                   Number* g_l, Number* g_u)
{
  bool retval;
  retval = sqp_tnlp_->get_bounds_info(n, x_l, x_u, m, g_l, g_u);
  return retval;
}

/*
 * @brief Get the starting point from the NLP object.
 */
bool IpoptSqpTNlp::get_starting_point(Index n, bool init_x, Number* x,
                                      bool init_z, Number* z_L, Number* z_U,
                                      Index m, bool init_lambda, Number* lambda)
{
  // The SQP solver has only one joint multiplier for the variable lower and
  // upper bounds.  We need to separate it into two
  double* bound_multipliers = nullptr;
  if (init_z) {
    bound_multipliers = new double[n];
  }

  bool retval;
  // Ipopt separates the bound multipliers into two parts
  retval = sqp_tnlp_->get_starting_point(
      n, init_x, x, init_z, bound_multipliers, m, init_lambda, lambda);

  if (!retval) {
    delete[] bound_multipliers;
    return retval;
  }

  if (init_z) {
    for (int i = 0; i < n; i++) {
      z_L[i] = max(bound_multipliers[i], 0.);
      z_U[i] = max(-bound_multipliers[i], 0.);
    }
    delete[] bound_multipliers;
  }

  // We also need to reverse the sign of the constraint multipliers
  if (init_lambda) {
    for (int i = 0; i < m; ++i) {
      lambda[i] = -lambda[i];
    }
  }

  return retval;
}

/**
 *@brief Evaluate the objective value
 */
bool IpoptSqpTNlp::eval_f(Index n, const Number* x, bool new_x,
                          Number& obj_value)
{
  return sqp_tnlp_->eval_objective_value(n, x, new_x, obj_value);
}

/**
 * @brief Evaluate the constraints at point x
 *
 */
bool IpoptSqpTNlp::eval_g(Index n, const Number* x, bool new_x, Index m,
                          Number* g)
{
  return sqp_tnlp_->eval_constraint_values(n, x, new_x, m, g);
}

/**
 *@brief Evaluate gradient at point x
 */
bool IpoptSqpTNlp::eval_grad_f(Index n, const Number* x, bool new_x,
                               Number* grad_f)
{
  return sqp_tnlp_->eval_objective_gradient(n, x, new_x, grad_f);
}

/**
 * @brief Get the matrix structure of the Jacobian
 * Always call this before the first time using @Eval_Jacobian
 */

bool IpoptSqpTNlp::eval_jac_g(Index n, const Number* x, bool new_x, Index m,
                              Index nele_jac, Index* iRow, Index* jCol,
                              Number* values)
{
  bool retval = sqp_tnlp_->eval_constraint_jacobian(n, x, new_x, m, nele_jac,
                                                    iRow, jCol, values);
  return retval;
}

/**
 * @brief Get the structure of the Hessian
 * Always call this before the first time using @Eval_Hessian
 */
bool IpoptSqpTNlp::eval_h(Index n, const Number* x, bool new_x,
                          Number obj_factor, Index m, const Number* lambda,
                          bool new_lambda, Index nele_hess, Index* iRow,
                          Index* jCol, Number* values)
{
  // For Ipopt, the constraint multipliers are defined with the opposite sign
  // So we need to create a copy of constraint multipliers with reversed sign

  double* neg_lambda = nullptr;
  if (lambda) {
    neg_lambda = new double[m];
    for (int i = 0; i < m; ++i) {
      neg_lambda[i] = -lambda[i];
    }
  }

  bool retval = sqp_tnlp_->eval_lagrangian_hessian(
      n, x, new_x, obj_factor, m, neg_lambda, new_lambda, nele_hess, iRow, jCol,
      values);

  delete[] neg_lambda;

  return retval;
}

void IpoptSqpTNlp::finalize_solution(SolverReturn ipopt_status, Index n,
                                     const Number* x, const Number* z_L,
                                     const Number* z_U, Index m,
                                     const Number* g, const Number* lambda,
                                     Number obj_value, const IpoptData* ip_data,
                                     IpoptCalculatedQuantities* ip_cq)
{
  ipopt_status_ = ipopt_status;

  // Delete previous information
  delete[] x_sol_;
  x_sol_ = nullptr;
  delete[] z_L_sol_;
  z_L_sol_ = nullptr;
  delete[] z_U_sol_;
  z_U_sol_ = nullptr;
  delete[] lambda_sol_;
  lambda_sol_ = nullptr;
  delete[] g_sol_;
  g_sol_ = nullptr;

  // Allocate memory for the new solution
  x_sol_ = new double[n];
  z_L_sol_ = new double[n];
  z_U_sol_ = new double[n];
  lambda_sol_ = new double[m];
  g_sol_ = new double[m];

  // Copy the values corresponding to the variables
  for (int i = 0; i < n; ++i) {
    x_sol_[i] = x[i];
    z_L_sol_[i] = z_L[i];
    z_U_sol_[i] = z_U[i];
  }

  // Copy the values corresponding to the constraints.
  // Make sure that the sign convention is that for the SQP algorithm
  for (int i = 0; i < m; ++i) {
    lambda_sol_[i] = -lambda[i];
    g_sol_[i] = g[i];
  }

  // Store objective function value
  obj_value_ = obj_value;
}
} // namespace RestartSqp
