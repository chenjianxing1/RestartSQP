/** Copyright (C) 2019
 * All Rights Reserved.
 *
 * Authors: Xinyi Luo
 * Date:2019-07
 */
#include "sqphot/QPhandler.hpp"
#include "IpOptionsList.hpp"
#include "sqphot/QOREInterface.hpp"
#include "sqphot/SQPDebug.hpp"
#include "sqphot/qpOASESInterface.hpp"

using namespace std;
using namespace Ipopt;

namespace SQPhotstart {

QPhandler::QPhandler(shared_ptr<const SqpNlpSizeInfo> nlp_sizes, QPType qptype,
                     SmartPtr<Journalist> jnlst,
                     SmartPtr<const OptionsList> options)
 : nlp_sizes_(nlp_sizes)
 , jnlst_(jnlst)
{
  // Initialize the matrix structure: where in the Jacobian are multiples of the
  // identity
  int num_vars = nlp_sizes->get_num_variables();
  int num_cons = nlp_sizes->get_num_constraints();
#ifdef NEW_FORMULATION
  num_qp_constraints_ = num_cons + num_vars;
  num_qp_variables_ = num_vars * 3 + 2 * num_cons;
  // Add the different matrices
  //
  // TODO which vars?
  identity_matrix_positions_.add_matrix(1, num_vars + 1, num_cons, 1.);
  // WHICH BLOCK?
  identity_matrix_positions_.add_matrix(1, num_vars + num_cons + 1, num_cons,
                                        -1.);
  // WHICH BLOCK?
  identity_matrix_positions_.add_matrix(num_cons + 1, 1, num_vars, 1.);
  // WHICH BLOCK?
  identity_matrix_positions_.add_matrix(
      num_cons + 1, num_vars + 2 * num_cons + 1, num_vars, 1.);
  // WHICH BLOCK
  identity_matrix_positions_.add_matrix(
      num_cons + 1, 2 * num_vars + 2 * num_cons + 1, num_vars, -1);
#if 0
  identity_matrix_positions_.length = 5;
  identity_matrix_positions_.irow = new int[5];
  identity_matrix_positions_.jcol = new int[5];
  identity_matrix_positions_.size = new int[5];
  identity_matrix_positions_.value = new double[5];
  identity_matrix_positions_.irow[0] = identity_matrix_positions_.irow[1] = 1;
  identity_matrix_positions_.irow[2] = identity_matrix_positions_.irow[3] = identity_matrix_positions_.irow[4] = num_cons + 1;
  identity_matrix_positions_.jcol[0] = num_vars + 1;
  identity_matrix_positions_.jcol[1] = num_vars + num_cons + 1;
  identity_matrix_positions_.jcol[2] = 1;
  identity_matrix_positions_.jcol[3] = num_vars + num_cons * 2 + 1;
  identity_matrix_positions_.jcol[4] = num_vars * 2 + num_cons * 2 + 1;
  identity_matrix_positions_.size[0] = identity_matrix_positions_.size[1] = num_cons;
  identity_matrix_positions_.size[2] = identity_matrix_positions_.size[3] = identity_matrix_positions_.size[4] = num_vars;
  identity_matrix_positions_.value[0] = identity_matrix_positions_.value[2] = identity_matrix_positions_.value[3] = 1.0;
  identity_matrix_positions_.value[1] = identity_matrix_positions_.value[4] = -1.0;
#endif
#else
  // TODO: make IdentityMatrixPosition a class
  num_qp_constraints_ = num_cons;
  num_qp_variables_ = num_vars + 2 * num_cons;

  // WHICH BLOCK
  identity_matrix_positions_.add_matrix(1, num_vars + 1, num_cons, 1.);
  identity_matrix_positions_.add_matrix(1, num_vars + num_cons + 1, num_cons,
                                        -1.);
#if 0
  identity_matrix_positions_.length = 2;
  identity_matrix_positions_.irow = new int[2];
  identity_matrix_positions_.jcol = new int[2];
  identity_matrix_positions_.size = new int[2];
  identity_matrix_positions_.value = new double[2];
  identity_matrix_positions_.irow[0] = identity_matrix_positions_.irow[1] = 1;
  identity_matrix_positions_.jcol[0] = num_vars + 1;
  identity_matrix_positions_.jcol[1] = num_vars + num_cons + 1;
  identity_matrix_positions_.size[0] = identity_matrix_positions_.size[1] = num_cons;
  identity_matrix_positions_.value[0] = 1.0;
  identity_matrix_positions_.value[1] = -1.0;
#endif
#endif

  // Allocate memory for the working sets
  bounds_working_set_ = new ActivityStatus[num_qp_variables_];
  constraints_working_set_ = new ActivityStatus[num_qp_constraints_];

  // Determine which solver should be used
  int enum_int;
  options->GetEnumValue("qp_solver_choice", enum_int, "");
  qp_solver_choice_ = QpSolver(enum_int);

  // Create solver objects
  switch (qp_solver_choice_) {
    case QPOASES:
      solverInterface_ =
          make_shared<qpOASESInterface>(nlp_sizes, qptype, options, jnlst);
      break;
    case QORE:
      solverInterface_ =
          make_shared<QOREInterface>(nlp_sizes, qptype, options, jnlst);
      break;
    case GUROBI:
#ifdef USE_GUROBI
      solverInterface_ =
          make_shared<GurobiInterface>(nlp_sizes, qptype, options, jnlst);
#endif
      break;
    case CPLEX:
#ifdef USE_CPLEX
      solverInterface_ =
          make_shared<CplexInterface>(nlp_sizes, qptype, options, jnlst);
#endif
      break;
  }

#ifdef DEBUG
#ifdef COMPARE_QP_SOLVER
  qpOASESInterface_ =
      make_shared<qpOASESInterface>(nlp_sizes, qptype, options, jnlst);
  QOREInterface_ =
      make_shared<QOREInterface>(nlp_sizes, qptype, options, jnlst);
  W_b_qpOASES_ = new ActivityStatus[num_qp_variables_];
  W_c_qpOASES_ = new ActivityStatus[num_qp_constraints_];
  W_b_qore_ = new ActivityStatus[num_qp_variables_];
  W_c_qore_ = new ActivityStatus[num_qp_constraints_];
#endif
#endif
}

/**
 *Default destructor
 */
QPhandler::~QPhandler()
{
  delete[] bounds_working_set_;
  bounds_working_set_ = NULL;
  delete[] constraints_working_set_;
  constraints_working_set_ = NULL;
#ifdef DEBUG
#ifdef COMPARE_QP_SOLVER
  delete[] W_b_qpOASES_;
  delete[] W_c_qpOASES_;
  delete[] W_b_qore_;
  delete[] W_c_qore_;
#endif
#endif
}

/**
 * Setup the bounds for the QP subproblems according to the information from
 * current
 * iterate. We have
 * 	c_l -c_k <=J_p+ u-v<=c_u-c_k
 * The bound is formulated as
 *   max(-delta, x_l-x_k)<=p<=min(delta,x_u-x_k)
 * and  u,v>=0
 * @param delta      trust region radius
 * @param x_k        current iterate point
 * @param c_k        current constraint value evaluated at x_k
 * @param x_l        the lower bounds for variables
 * @param x_u        the upper bounds for variables
 * @param c_l        the lower bounds for constraints
 * @param c_u        the upper bounds for constraints
 */

void QPhandler::set_bounds(double delta, shared_ptr<const Vector> x_l,
                           shared_ptr<const Vector> x_u,
                           shared_ptr<const Vector> x_k,
                           shared_ptr<const Vector> c_l,
                           shared_ptr<const Vector> c_u,
                           shared_ptr<const Vector> c_k)
{

#ifdef DEBUG
#ifdef COMPARE_QP_SOLVER
  set_bounds_debug(delta, x_l, x_u, x_k, c_l, c_u, c_k);
#endif
#endif

  /*-------------------------------------------------------------*/
  /* Set lbA, ubA as well as lb and ub as qpOASES differentiates */
  /*the bound constraints from the linear constraints            */
  /*-------------------------------------------------------------*/
  int num_vars = nlp_sizes_->get_num_variables();
  int num_cons = nlp_sizes_->get_num_constraints();

  // TODO: There should be no difference between qpOASES and QORE

  if (qp_solver_choice_ != QORE) {
#ifdef NEW_FORMULATION
    const double* c_k_vals = c_k->get_values();
    const double* c_l_vals = c_l->get_values();
    const double* c_u_vals = c_u->get_values();
    for (int i = 0; i < num_cons; i++) {
      solverInterface_->set_lower_constraint_bounds(i, c_l_vals[i] - c_k_vals[i]);
      solverInterface_->set_upper_constraint_bounds(i, c_u_vals[i] - c_k_vals[i]);
    }
    const double* x_k_vals = x_k->get_values();
    const double* x_l_vals = x_l->get_values();
    const double* x_u_vals = x_u->get_values();
    for (int i = 0; i < num_vars; i++) {
      solverInterface_->set_lower_constraint_bounds(num_cons + i, x_l_vals[i] - x_k_vals[i]);
      solverInterface_->set_upper_constraint_bounds(num_cons + i, x_u_vals[i] - x_k_vals[i]);
    }

    for (int i = 0; i < num_vars; i++) {
      solverInterface_->set_lower_variable_bounds(i, -delta);
      solverInterface_->set_upper_variable_bounds(i, delta);
    }
    for (int i = num_vars; i < num_qp_variables_; i++) {
      solverInterface_->set_lower_variable_bounds(i, 0.);
      solverInterface_->set_upper_variable_bounds(i, INF);
    }
#else
    const double* c_k_vals = c_k->get_values();
    const double* c_l_vals = c_l->get_values();
    const double* c_u_vals = c_u->get_values();
    for (int i = 0; i < num_cons; i++) {
      solverInterface_->set_lbA(i, c_l_vals[i] - c_k_vals[i]);
      solverInterface_->set_ubA(i, c_u_vals[i] - c_k_vals[i]);
    }
    const double* x_k_vals = x_k->get_values();
    const double* x_l_vals = x_l->get_values();
    const double* x_u_vals = x_u->get_values();
    for (int i = 0; i < num_vars; i++) {
      solverInterface_->set_lb(i, max(x_l_vals[i] - x_k_vals[i], -delta));
      solverInterface_->set_ub(i, min(x_u_vals[i] - x_k_vals[i], delta));
    }
    /**
     * only set the upper bound for the last half to be infinity(those are slack
     * variables).
     * The lower bounds are initialized as 0
     */
    for (int i = 0; i < num_cons * 2; i++) {
      solverInterface_->set_lb(num_vars + i, 0.);
      solverInterface_->set_ub(num_vars + i, INF);
    }
#endif
  }
  /*-------------------------------------------------------------*/
  /* Only set lb and ub, where lb = [lbx;lbA]; and ub=[ubx; ubA] */
  /*-------------------------------------------------------------*/
  else {
#ifndef NEW_FORMULATION
    // For p_k
    const double* x_k_vals = x_k->get_values();
    const double* x_l_vals = x_l->get_values();
    const double* x_u_vals = x_u->get_values();
    for (int i = 0; i < num_vars; i++) {
      solverInterface_->set_lb(i, max(x_l_vals[i] - x_k_vals[i], -delta));
      solverInterface_->set_ub(i, min(x_u_vals[i] - x_k_vals[i], delta));
    }

    // For u,v
    for (int i = 0; i < num_cons * 2; i++) {
      solverInterface_->set_lb(num_vars + i, 0.);
      solverInterface_->set_ub(num_vars + i, INF);
    }
#endif

    const double* c_k_vals = c_k->get_values();
    const double* c_l_vals = c_l->get_values();
    const double* c_u_vals = c_u->get_values();
    for (int i = 0; i < num_cons; i++) {
      solverInterface_->set_lower_variable_bounds(num_qp_variables_ + i,
                               c_l_vals[i] - c_k_vals[i]);
      solverInterface_->set_upper_variable_bounds(num_qp_variables_ + i,
                               c_u_vals[i] - c_k_vals[i]);
    }

#ifdef NEW_FORMULATION
    const double* x_k_vals = x_k->get_values();
    const double* x_l_vals = x_l->get_values();
    const double* x_u_vals = x_u->get_values();
    for (int i = 0; i < num_vars; i++) {
      solverInterface_->set_lower_variable_bounds(i, -delta);
      solverInterface_->set_upper_variable_bounds(i, delta);

      solverInterface_->set_lower_variable_bounds(num_qp_variables_ + num_cons + i,
                               x_l_vals[i] - x_k_vals[i]);
      solverInterface_->set_upper_variable_bounds(num_qp_variables_ + num_cons + i,
                               x_u_vals[i] - x_k_vals[i]);
    }
    for (int i = num_vars; i < num_qp_variables_; i++) {
      solverInterface_->set_upper_variable_bounds(i, INF);
    }
// DEBUG
//        solverInterface_->getLb()->print("lb");
//        solverInterface_->getUb()->print("ub");
#endif
  }
}

/**
 * This function sets up the object vector g of the QP problem
 * The (2*nCon+nVar) vector g_^T in QP problem will be the same as
 * [grad_f^T, rho* e^T], where the unit vector is of length (2*nCon).
 * @param grad      Gradient vector from nlp class
 * @param rho       Penalty Parameter
 */

void QPhandler::set_gradient(shared_ptr<const Vector> grad, double rho)
{
  int num_vars = nlp_sizes_->get_num_variables();
#ifdef DEBUG
#ifdef COMPARE_QP_SOLVER
  for (int i = 0; i < num_qp_variables_; i++)
    if (i < num_vars) {
      qpOASESInterface_->set_gradient(i, grad->get_value(i));
      QOREInterface_->set_gradient(i, grad->get_value(i));
    } else {
      qpOASESInterface_->set_gradient(i, rho);
      QOREInterface_->set_gradient(i, rho);
    }
#endif
#endif
  for (int i = 0; i < num_vars; i++) {
    solverInterface_->set_linear_objective_coefficients(i, grad->get_value(i));
  }
  for (int i = num_vars; i < num_qp_variables_; i++) {
    solverInterface_->set_linear_objective_coefficients(i, rho);
  }
  // DEBUG
  //    solverInterface_->getG()->print("G");
}

/**
 * @brief Set up the H for the first time in the QP problem.
 * It will be concatenated as [H_k 0]
 *          		         [0   0]
 * where H_k is the Lagragian hessian evaluated at x_k and lambda_k.
 *
 * This method should only be called for once.
 *
 * @param hessian the Lagragian hessian evaluated at x_k and lambda_k from nlp
 * readers.
 */
void QPhandler::set_hessian(shared_ptr<const SpTripletMat> hessian)
{
#ifdef DEBUG
#ifdef COMPARE_QP_SOLVER
  qpOASESInterface_->set_hessian(hessian);
  QOREInterface_->set_hessian(hessian);
#endif
#endif
  solverInterface_->set_objective_hessian(hessian);
}

/**
 * @brief This function sets up the matrix A in the QP subproblem
 * The matrix A in QP problem will be concatenate as [J I -I]
 * @param jacobian  the Matrix object for Jacobian from c(x)
 */
void QPhandler::set_jacobian(shared_ptr<const SpTripletMat> jacobian)
{
#ifdef DEBUG
#ifdef COMPARE_QP_SOLVER
  qpOASESInterface_->set_jacobian(jacobian, identity_matrix_positions_);
  QOREInterface_->set_jacobian(jacobian, identity_matrix_positions_);
#endif
#endif
  solverInterface_->set_constraint_jacobian(jacobian, identity_matrix_positions_);
}

/**
 * @brief This function updates the constraint if there is any changes to
 * the iterates
 */
void QPhandler::update_bounds(double delta, shared_ptr<const Vector> x_l,
                              shared_ptr<const Vector> x_u,
                              shared_ptr<const Vector> x_k,
                              shared_ptr<const Vector> c_l,
                              shared_ptr<const Vector> c_u,
                              shared_ptr<const Vector> c_k)
{
#ifdef DEBUG
#ifdef COMPARE_QP_SOLVER
  set_bounds_debug(delta, x_l, x_u, x_k, c_l, c_u, c_k);
#endif
#endif
  int num_vars = nlp_sizes_->get_num_variables();
  int num_cons = nlp_sizes_->get_num_constraints();
#ifndef NEW_FORMULATION
  if (qp_solver_choice_ != QORE) {
    if (qp_solver_choice_ == GUROBI || qp_solver_choice_ == CPLEX) {
      solverInterface_->reset_constraints();
    }
    for (int i = 0; i < num_cons; i++) {
      solverInterface_->set_lbA(i, c_l->get_value(i) - c_k->get_value(i));
    }

    for (int i = 0; i < num_vars; i++) {
      solverInterface_->set_lb(
          i, max(x_l->get_value(i) - x_k->get_value(i), -delta));
      solverInterface_->set_ub(
          i, min(x_u->get_value(i) - x_k->get_value(i), delta));
    }
  } else {
    for (int i = 0; i < num_vars; i++) {
      solverInterface_->set_lb(
          i, max(x_l->get_value(i) - x_k->get_value(i), -delta));
      solverInterface_->set_ub(
          i, min(x_u->get_value(i) - x_k->get_value(i), delta));
    }
    for (int i = 0; i < num_cons; i++) {
      solverInterface_->set_lb(num_vars + 2 * num_cons + i,
                               c_l->get_value(i) - c_k->get_value(i));
      solverInterface_->set_ub(num_vars + 2 * num_cons + i,
                               c_u->get_value(i) - c_k->get_value(i));
    }
  }
#else
  if (qp_solver_choice_ == QORE) {
    for (int i = 0; i < num_vars; i++) {
      solverInterface_->set_lower_variable_bounds(i, -delta);
      solverInterface_->set_upper_variable_bounds(i, delta);
    }
    for (int i = 0; i < num_cons; i++) {
      solverInterface_->set_lower_variable_bounds(num_qp_variables_ + i,
                               c_l->get_value(i) - c_k->get_value(i));
      solverInterface_->set_upper_variable_bounds(num_qp_variables_ + i,
                               c_u->get_value(i) - c_k->get_value(i));
    }

    for (int i = 0; i < num_vars; i++) {
      solverInterface_->set_lower_variable_bounds(num_qp_variables_ + num_cons + i,
                               x_l->get_value(i) - x_k->get_value(i));
      solverInterface_->set_upper_variable_bounds(num_qp_variables_ + num_cons + i,
                               x_u->get_value(i) - x_k->get_value(i));
    }
  } else {
    for (int i = 0; i < num_vars; i++) {
      solverInterface_->set_lower_variable_bounds(i, -delta);
      solverInterface_->set_upper_variable_bounds(i, delta);
    }

    for (int i = 0; i < num_cons; i++) {
      solverInterface_->set_lower_constraint_bounds(i, c_l->get_value(i) - c_k->get_value(i));
      solverInterface_->set_upper_constraint_bounds(i, c_u->get_value(i) - c_k->get_value(i));
    }
    for (int i = 0; i < num_vars; i++) {
      solverInterface_->set_lower_constraint_bounds(num_cons + i,
                                x_l->get_value(i) - x_k->get_value(i));
      solverInterface_->set_upper_constraint_bounds(num_cons + i,
                                x_u->get_value(i) - x_k->get_value(i));
    }
  }
#endif
}

/**
 * @brief This function updates the vector g in the QP subproblem when there are
 * any
 * change to the values of penalty parameter
 *
 * @param rho               penalty parameter
 * @param nVar              number of variables in NLP
 */

void QPhandler::update_penalty(double rho)
{
  int num_vars = nlp_sizes_->get_num_variables();
#ifdef DEBUG
#ifdef COMPARE_QP_SOLVER
  int num_cons = nlp_sizes_->get_num_constraints();
  for (int i = num_vars; i < num_qp_variables_; i++) {
    qpOASESInterface_->set_gradient(i, rho);
    QOREInterface_->set_gradient(i, rho);
  }
#endif
#endif
  for (int i = num_vars; i < num_qp_variables_; i++)
    solverInterface_->set_linear_objective_coefficients(i, rho);
}

/**
 * @brief This function updates the vector g in the QP subproblem when there are
 * any
 * change to the values of gradient in NLP
 *
 * @param grad              the gradient vector from NLP
 */
void QPhandler::update_grad(shared_ptr<const Vector> grad)
{
  int num_vars = nlp_sizes_->get_num_variables();
#ifdef DEBUG
#ifdef COMPARE_QP_SOLVER

  for (int i = 0; i < num_vars; i++) {
    qpOASESInterface_->set_gradient(i, grad->get_value(i));
    QOREInterface_->set_gradient(i, grad->get_value(i));
  }
#endif
#endif
  for (int i = 0; i < num_vars; i++)
    solverInterface_->set_linear_objective_coefficients(i, grad->get_value(i));
}

/**
 *@brief Solve the QP with objective and constraints defined by its class
 *members
 */
QpSolverExitStatus QPhandler::solve_qp(shared_ptr<Statistics> stats)
{

//    solverInterface_->getA()->print_full("A");
//    solverInterface_->getH()->print_full("H");
//    solverInterface_->getLb()->print("Lb");
//    solverInterface_->getUb()->print("Ub");
//    solverInterface_->getG()->print("G");

#ifdef DEBUG
#ifdef COMPARE_QP_SOLVER
  QOREInterface_->optimizeQP(stats);
  qpOASESInterface_->optimizeQP(stats);
  bool qpOASES_optimal =
      test_optimality(qpOASESInterface_, QPOASES, W_b_qpOASES_, W_c_qpOASES_);
  bool qore_optimal =
      test_optimality(QOREInterface_, QORE, W_b_qore_, W_c_qore_);
  if (!qpOASES_optimal || !qore_optimal) {
    testQPsolverDifference();
  }

#endif
#endif

  QpSolverExitStatus qp_exit_status = solverInterface_->optimize_qp(stats);

  // Check if the QP solver really returned an optimal solution.
  if (qp_exit_status == QPEXIT_OPTIMAL) {
    bool isOptimal =
        test_optimality(solverInterface_, qp_solver_choice_,
                        bounds_working_set_, constraints_working_set_);
    assert(isOptimal);
  }

  return qp_exit_status;
}

double QPhandler::get_objective()
{

  return solverInterface_->get_obj_value();
}

void QPhandler::update_H(shared_ptr<const SpTripletMat> Hessian)
{

#ifdef DEBUG
#ifdef COMPARE_QP_SOLVER
  QOREInterface_->set_hessian(Hessian);
  qpOASESInterface_->set_hessian(Hessian);
#endif
#endif
  solverInterface_->set_objective_hessian(Hessian);
}

void QPhandler::update_A(shared_ptr<const SpTripletMat> Jacobian)
{

#ifdef DEBUG
#ifdef COMPARE_QP_SOLVER
  QOREInterface_->set_jacobian(Jacobian, identity_matrix_positions_);
  qpOASESInterface_->set_jacobian(Jacobian, identity_matrix_positions_);
#endif
#endif
  solverInterface_->set_constraint_jacobian(Jacobian, identity_matrix_positions_);
}

void QPhandler::update_delta(double delta, shared_ptr<const Vector> x_l,
                             shared_ptr<const Vector> x_u,
                             shared_ptr<const Vector> x_k)
{
  int num_vars = nlp_sizes_->get_num_variables();
#ifdef DEBUG
#ifdef COMPARE_QP_SOLVER
  for (int i = 0; i < num_vars; i++) {
    qpOASESInterface_->set_lb(
        i, max(x_l->get_value(i) - x_k->get_value(i), -delta));
    qpOASESInterface_->set_ub(
        i, min(x_u->get_value(i) - x_k->get_value(i), delta));
    QOREInterface_->set_lb(i,
                           max(x_l->get_value(i) - x_k->get_value(i), -delta));
    QOREInterface_->set_ub(i,
                           min(x_u->get_value(i) - x_k->get_value(i), delta));
  }
#endif
#endif

#ifdef NEW_FORMULATION
  for (int i = 0; i < num_vars; i++) {
    solverInterface_->set_lower_variable_bounds(i, -delta);
    solverInterface_->set_upper_variable_bounds(i, delta);
  }
#else

  for (int i = 0; i < num_vars; i++) {
    solverInterface_->set_lb(
        i, max(x_l->get_value(i) - x_k->get_value(i), -delta));
    solverInterface_->set_ub(i,
                             min(x_u->get_value(i) - x_k->get_value(i), delta));
  }

#endif
}

void QPhandler::write_qp_data(const string filename)
{

  solverInterface_->WriteQPDataToFile(J_LAST_LEVEL, J_USER1, filename);
}

Exitflag QPhandler::get_status()
{
  return (solverInterface_->get_status());
}

bool QPhandler::test_optimality(shared_ptr<QPSolverInterface> qpsolverInterface,
                                QpSolver qpSolver,
                                ActivityStatus* bounds_working_set,
                                ActivityStatus* constraints_working_set)
{
  qpOptimalStatus_ = qpsolverInterface->get_optimality_status();
  return (qpsolverInterface->test_optimality(constraints_working_set,
                                             bounds_working_set));
}

double QPhandler::get_model_infeasibility() const
{
  // Here we compute the constraint violation from the penatly values.
  // The QP primal variables are ordered so that the SQP primal variables
  // are first, and the remaining ones are the penalty variables.
  int num_vars = nlp_sizes_->get_num_variables();
  return solverInterface_->get_primal_solution()->calc_subvector_one_norm(
      num_vars, num_qp_variables_ - num_vars);
}

const OptimalityStatus& QPhandler::get_QpOptimalStatus() const
{
  return qpOptimalStatus_;
}

void QPhandler::get_active_set(ActivityStatus* A_c, ActivityStatus* A_b,
                               shared_ptr<Vector> x, shared_ptr<Vector> Ax)
{
  // use the class member to get the qp problem information
  auto lb = solverInterface_->get_lower_variable_bounds();
  auto ub = solverInterface_->get_upper_variable_bounds();
  if (x == nullptr) {
    x = make_shared<Vector>(num_qp_variables_);
    x->copy_vector(get_primal_solution());
  }
  if (Ax == nullptr) {
    Ax = make_shared<Vector>(num_qp_constraints_);
    auto A = solverInterface_->get_constraint_jacobian();
    A->multiply(x, Ax);
  }

  for (int i = 0; i < num_qp_variables_; i++) {
    if (abs(x->get_value(i) - lb->get_value(i)) < sqrt_m_eps) {
      if (abs(ub->get_value(i) - x->get_value(i)) < sqrt_m_eps) {
        A_b[i] = ACTIVE_BOTH_SIDES;
      } else {
        A_b[i] = ACTIVE_BELOW;
      }
    } else if (abs(ub->get_value(i) - x->get_value(i)) < sqrt_m_eps) {
      A_b[i] = ACTIVE_ABOVE;
    } else {
      A_b[i] = INACTIVE;
    }
  }
  if (qp_solver_choice_ == QORE) {
    // if no x and Ax are input
    for (int i = 0; i < num_qp_constraints_; i++) {
      if (abs(Ax->get_value(i) - lb->get_value(i + num_qp_variables_)) <
          sqrt_m_eps) {
        if (abs(ub->get_value(i + num_qp_variables_) - Ax->get_value(i)) <
            sqrt_m_eps) {
          A_c[i] = ACTIVE_BOTH_SIDES;
        } else {
          A_c[i] = ACTIVE_BELOW;
        }
      } else if (abs(ub->get_value(i + num_qp_variables_) - Ax->get_value(i)) <
                 sqrt_m_eps) {
        A_c[i] = ACTIVE_ABOVE;
      } else {
        A_c[i] = INACTIVE;
      }
    }
  } else {
    auto lbA = solverInterface_->get_upper_constraint_bounds();
    auto ubA = solverInterface_->get_upper_constraint_bounds();
    for (int i = 0; i < num_qp_constraints_; i++) {
      if (abs(Ax->get_value(i) - lbA->get_value(i)) < sqrt_m_eps) {
        if (abs(ubA->get_value(i) - Ax->get_value(i)) < sqrt_m_eps) {
          A_c[i] = ACTIVE_BOTH_SIDES;
        } else {
          A_c[i] = ACTIVE_BELOW;
        }
      } else if (abs(ubA->get_value(i) - Ax->get_value(i)) < sqrt_m_eps) {
        A_c[i] = ACTIVE_ABOVE;
      } else {
        A_c[i] = INACTIVE;
      }
    }
  }
}

void QPhandler::set_gradient(double rho)
{
  int num_vars = nlp_sizes_->get_num_variables();
  for (int i = 0; i < num_vars; i++) {
    solverInterface_->set_linear_objective_coefficients(i, rho);
  }
  for (int i = num_vars; i < num_qp_variables_; i++) {
    solverInterface_->set_linear_objective_coefficients(i, rho);
  }
}

#ifdef DEBUG
#ifdef COMPARE_QP_SOLVER
void QPhandler::set_bounds_debug(double delta, shared_ptr<const Vector> x_l,
                                 shared_ptr<const Vector> x_u,
                                 shared_ptr<const Vector> x_k,
                                 shared_ptr<const Vector> c_l,
                                 shared_ptr<const Vector> c_u,
                                 shared_ptr<const Vector> c_k)
{
  int num_vars = nlp_sizes_->get_num_variables();
  int num_cons = nlp_sizes_->get_num_constraints();

  for (int i = 0; i < num_cons; i++) {
    qpOASESInterface_->set_lbA(i, c_l->get_value(i) - c_k->get_value(i));
    qpOASESInterface_->set_ubA(i, c_u->get_value(i) - c_k->get_value(i));
  }

  for (int i = 0; i < num_vars; i++) {
    qpOASESInterface_->set_lb(
        i, max(x_l->get_value(i) - x_k->get_value(i), -delta));
    qpOASESInterface_->set_ub(
        i, min(x_u->get_value(i) - x_k->get_value(i), delta));
  }
  /**
   * only set the upper bound for the last 2*nCon entries (those are slack
   * variables).
   * The lower bounds are initialized as 0
   */
  for (int i = 0; i < num_cons * 2; i++) {
    qpOASESInterface_->set_ub(num_vars + i, INF);
  }

  /*-------------------------------------------------------------*/
  /* Only set lb and ub, where lb = [lbx;lbA]; and ub=[ubx; ubA] */
  /*-------------------------------------------------------------*/
  for (int i = 0; i < num_vars; i++) {
    QOREInterface_->set_lb(i,
                           max(x_l->get_value(i) - x_k->get_value(i), -delta));
    QOREInterface_->set_ub(i,
                           min(x_u->get_value(i) - x_k->get_value(i), delta));
  }

  for (int i = 0; i < num_cons * 2; i++) {
    QOREInterface_->set_ub(num_vars + i, INF);
  }

  for (int i = 0; i < num_cons; i++) {
    QOREInterface_->set_lb(num_qp_variables_ + i,
                           c_l->get_value(i) - c_k->get_value(i));
    QOREInterface_->set_ub(num_qp_variables_ + i,
                           c_u->get_value(i) - c_k->get_value(i));
  }
}

bool QPhandler::testQPsolverDifference()
{
  int num_vars = nlp_sizes_->get_num_variables();
  int num_cons = nlp_sizes_->get_num_constraints();

  shared_ptr<Vector> qpOASESsol = make_shared<Vector>(num_vars);
  shared_ptr<Vector> QOREsol = make_shared<Vector>(num_vars);
  shared_ptr<Vector> difference = make_shared<Vector>(num_vars);
  qpOASESsol->print("qpOASESsol", jnlst_);
  QOREsol->print("QOREsol", jnlst_);
  qpOASESsol->copy_vector(qpOASESInterface_->get_optimal_solution());
  QOREsol->copy_vector(QOREInterface_->get_optimal_solution());
  difference->copy_vector(qpOASESsol);
  difference->add_vector(-1., QOREsol);
  double diff_norm = difference->calc_one_norm();
  if (diff_norm > 1.0e-8) {
    printf("difference is %10e\n", diff_norm);
    qpOASESsol->print("qpOASESsol");
    QOREsol->print("QOREsol");
    // TODO activate this again
    // QOREInterface_->WriteQPDataToFile(jnlst_, J_ALL, J_DBG);
  }
  assert(diff_norm < 1.0e-8);
  return true;
}

#endif
#endif
} // namespace SQPhotstart
