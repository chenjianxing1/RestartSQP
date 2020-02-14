/* Copyright (C) 2019
* All Rights Reserved.
*
* Authors: Xinyi Luo
* Date:2019-07
*/
#include "sqphot/IpoptSqpNlp.hpp"
#include <cassert>

using namespace std;
using namespace Ipopt;

namespace RestartSqp {

/** Default constructor*/
IpoptSqpNlp::IpoptSqpNlp(std::shared_ptr<SqpTNlp> sqp_tnlp)
 : sqp_tnlp_(sqp_tnlp)
{
}

/** Default constructor*/
IpoptSqpNlp::~IpoptSqpNlp()
{
}

bool IpoptSqpNlp::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                               Index& nnz_h_lag,
                               TNLP::IndexStyleEnum& index_style)
{
  string nlp_name = "IpoptTNLP";
  bool retval;
  retval = sqp_tnlp_->get_nlp_info(n, m, nnz_jac_g, nnz_h_lag, nlp_name);

  index_style == TNLP::FORTRAN_STYLE;

  return retval;
}

/**
 *@brief get the bounds information from the NLP object
 */
bool IpoptSqpNlp::get_bounds_info(Index n, Number* x_l, Number* x_u, Index m,
                                  Number* g_l, Number* g_u)
{
  bool retval;
  retval = sqp_tnlp_->get_bounds_info(n, x_l, x_u, m, g_l, g_u);
  return retval;
}

/*
 * @brief Get the starting point from the NLP object.
 */
bool IpoptSqpNlp::get_starting_point(Index n, bool init_x, Number* x,
                                     bool init_z, Number* z_L, Number* z_U,
                                     Index m, bool init_lambda, Number* lambda)
{
  // The SQP solver has only one joint multiplier for the variable lower and
  // upper bounds.  We need to separate it into two
  double* bound_multipliers = nullptr;
  if (init_z) {
    bound_multipliers = new double[n];
    for (int i = 0; i < n; i++) {
      bound_multipliers[i] = z_L[i] - z_U[i];
    }
  }

  bool retval;
  // Ipopt separates the bound multipliers into two parts
  retval = sqp_tnlp_->get_starting_point(
      n, init_x, x, init_z, bound_multipliers, m, init_lambda, lambda);

  delete[] bound_multipliers;

  if (!retval) {
    return retval;
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
bool IpoptSqpNlp::eval_f(Index n, const Number* x, bool new_x,
                         Number& obj_value)
{
  return sqp_tnlp_->eval_objective_value(n, x, new_x, obj_value);
}

/**
 * @brief Evaluate the constraints at point x
 *
 */
bool IpoptSqpNlp::eval_g(Index n, const Number* x, bool new_x, Index m,
                         Number* g)
{
  return sqp_tnlp_->eval_constraint_values(n, x, new_x, m, g);
}

/**
 *@brief Evaluate gradient at point x
 */
bool IpoptSqpNlp::eval_grad_f(Index n, const Number* x, bool new_x,
                              Number* grad_f)
{
  return sqp_tnlp_->eval_objective_gradient(n, x, new_x, grad_f);
}

/**
 * @brief Get the matrix structure of the Jacobian
 * Always call this before the first time using @Eval_Jacobian
 */

bool IpoptSqpNlp::eval_jac_g(Index n, const Number* x, bool new_x, Index m,
                             Index nele_jac, Index* iRow, Index* jCol,
                             Number* values)
{
  return sqp_tnlp_->eval_constraint_jacobian(n, x, new_x, m, nele_jac, iRow,
                                             jCol, values);
}

/**
 * @brief Get the structure of the Hessian
 * Always call this before the first time using @Eval_Hessian
 */
bool IpoptSqpNlp::eval_h(Index n, const Number* x, bool new_x,
                         Number obj_factor, Index m, const Number* lambda,
                         bool new_lambda, Index nele_hess, Index* iRow,
                         Index* jCol, Number* values)
{
  // For Ipopt, the constraint multipliers are defined with the opposite sign
  // So we need to create a copy of constraint multipliers with reversed sign
  double* neg_lambda = new double[m];
  for (int i = 0; i < m; ++i) {
    neg_lambda[i] = -lambda[i];
  }

  bool retval = sqp_tnlp_->eval_lagrangian_hessian(
      n, x, new_x, obj_factor, m, neg_lambda, new_lambda, nele_hess, iRow, jCol,
      values);

  delete[] neg_lambda;

  return retval;
}

void IpoptSqpNlp::finalize_solution(SolverReturn ipopt_status, Index n,
                                    const Number* x, const Number* z_L,
                                    const Number* z_U, Index m, const Number* g,
                                    const Number* lambda, Number obj_value,
                                    const IpoptData* ip_data,
                                    IpoptCalculatedQuantities* ip_cq)
{
  // Translate SQP exit status to Ipopt exit status
  SqpSolverExitStatus sqp_status;

  switch (ipopt_status) {
    case SUCCESS:
      sqp_status = OPTIMAL;
      break;
    case MAXITER_EXCEEDED:
      sqp_status = EXCEED_MAX_ITERATIONS;
      break;
    case STOP_AT_TINY_STEP:
      sqp_status = TRUST_REGION_TOO_SMALL;
      break;
    case CPUTIME_EXCEEDED:
      sqp_status = EXCEED_MAX_CPU_TIME;
      break;
    case RESTORATION_FAILURE:
      sqp_status = PENALTY_TOO_LARGE;
      break;
    default:
      sqp_status = UNKNOWN_EXIT_STATUS;
      break;
  }

  // The SQP solvers combines the bound multipliers into a single array
  double* bound_multipliers = new double[n];
  for (int i = 0; i < n; ++i) {
    bound_multipliers[i] = z_L[i] - z_U[i];
  }

  // We also need to invert the sign of the constraint multipliers
  double* constraint_multipliers = new double[m];

  // Create copy of constraint multipliers with reversed sign
  double* neg_lambda = new double[m];
  for (int i = 0; i < m; ++i) {
    constraint_multipliers[i] = -lambda[i];
  }

  // We don't have information on the activities or the solver statistics
  const ActivityStatus* bound_activity_status = nullptr;
  const ActivityStatus* constraint_activity_status = nullptr;
  shared_ptr<const Statistics> stats;

  // Call the Ipopt method
  sqp_tnlp_->finalize_solution(sqp_status, n, x, bound_multipliers,
                               bound_activity_status, m, g, lambda,
                               constraint_activity_status, obj_value, stats);

  delete[] bound_multipliers;
  delete[] constraint_multipliers;
}
}
