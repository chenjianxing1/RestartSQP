/* Copyright (C) 2019
 * All Rights Reserved.
 *
 * Authors: Xinyi Luo
 * Date:    2019-07
 */
#include "sqphot/qpOASESInterface.hpp"
#include "sqphot/MessageHandling.hpp"

using namespace std;
using namespace Ipopt;

namespace SQPhotstart {

/**
 * @brief Constructor which also initializes the qpOASES SQProblem objects
 *
 * It will either formulate the SQP subproblem as
 * 		\minimize 1/2*p^T*H_k*p + g_k^T*p + rho*(e^T u_1 + e^Tu_2 +e^T
 * v_1
 * e^T v_2)
 * 		\subject  c_l-c_k <= J_k*p + u_1 - u_2 <= c_u-c_k,
 * 			  x_l-x_k <= p + v_1 -v_2 <= x_u -x_k, *
 * 			  ||p||<=\delta,
 * 			  u_1,u_2,v_1,v_2>=0.
 * Or
 * 		\minimize 1/2*p^T*H_k*p + g_k^T*p + rho*(e^T u_1 + e^Tu_2)
 * 		\subject  c_l-c_k <= J_k*p + u_1 - u_2 <= c_u-c_k,
 * 			  x_l-x_k <= p <= x_u -x_k, *
 * 			  ||p||<=\delta,
 * 			  u_1,u_2.
 *
 * depending on user's choice.
 *
 * @param nlp_info the struct that stores simple nlp dimension info
 * @param qptype  is the problem to be solved QP or LP or SOC?
 */
qpOASESInterface::qpOASESInterface(
    std::shared_ptr<const SqpNlpSizeInfo> nlp_sizes, QPType qptype,
    SmartPtr<const OptionsList> options, SmartPtr<Journalist> jnlst)
 : jnlst_(jnlst)
 , first_qp_solved_(false)
 , qp_matrices_changed_(true)
{
  // Get the option values for this object
  get_option_values_(options);

  num_nlp_variables_ = nlp_sizes->get_num_variables();
  num_nlp_constraints_ = nlp_sizes->get_num_constraints();

#ifdef NEW_FORMULATION
  num_qp_constraints_ =
      nlp_sizes->get_num_constraints() + nlp_sizes->get_num_variables();
  num_qp_variables_ =
      nlp_sizes->get_num_variables() * 3 + 2 * nlp_sizes->get_num_constraints();
#else
  num_qp_constraints_ = nlp_sizes->get_num_constraints();
  num_qp_variables_ =
      nlp_sizes->get_num_variables() + 2 * nlp_sizes->get_num_constraints();
#endif
  allocate_memory(nlp_sizes, qptype);
}

void qpOASESInterface::get_option_values_(SmartPtr<const OptionsList> options)
{
  // Get the options from the options list
  options->GetIntegerValue("qp_solver_max_num_iterations",
                           qp_solver_max_num_iterations_, "");
  options->GetIntegerValue("qp_solver_print_level", qp_solver_print_level_, "");
  options->GetIntegerValue("lp_solver_max_num_iterations",
                           lp_solver_max_num_iterations_, "");
}

qpOASESInterface::qpOASESInterface(shared_ptr<SparseHbMatrix> H,
                                   shared_ptr<SparseHbMatrix> A,
                                   shared_ptr<Vector> g, shared_ptr<Vector> lb,
                                   shared_ptr<Vector> ub,
                                   shared_ptr<Vector> lbA,
                                   shared_ptr<Vector> ubA,
                                   SmartPtr<const OptionsList> options)
 : num_qp_variables_(A->get_num_columns())
 , num_qp_constraints_(A->get_num_rows())
 , first_qp_solved_(false)
 , H_(H)
 , A_(A)
 , g_(g)
 , lb_(lb)
 , ub_(ub)
 , lbA_(lbA)
 , ubA_(ubA)
{
  A_->print("bla");
  // Get the option values for this object
  get_option_values_(options);

  x_qp_ = make_shared<Vector>(num_qp_variables_);
  y_qp_ = make_shared<Vector>(num_qp_constraints_ + num_qp_variables_);
  qpoases_solver_ = std::make_shared<qpOASES::SQProblem>(
      (qpOASES::int_t)num_qp_variables_, (qpOASES::int_t)num_qp_constraints_);
  // AW: The following is very bad design in qpOASES: the constructor does not
  // copy
  // the values, but the only time the values
  // are actually changed is if the setVal method is called in
  // qpOASES::SparseMatrix.  It should be const.
  int* H_const_row_indices = const_cast<int*>(H_->get_row_indices());
  int* H_const_col_indices = const_cast<int*>(H_->get_column_indices());
  double* H_const_values = const_cast<double*>(H_->get_values());
  H_qpOASES_ = std::make_shared<qpOASES::SymSparseMat>(
      num_qp_variables_, num_qp_variables_, H_const_row_indices,
      H_const_col_indices, H_const_values);
  H_qpOASES_->createDiagInfo();

  int* A_const_row_indices = const_cast<int*>(A_->get_row_indices());
  int* A_const_col_indices = const_cast<int*>(A_->get_column_indices());
  double* A_const_values = const_cast<double*>(A_->get_values());
  A_qpOASES_ = std::make_shared<qpOASES::SymSparseMat>(
      num_qp_constraints_, num_qp_variables_, A_const_row_indices,
      A_const_col_indices, A_const_values);
}

/**Default destructor*/
qpOASESInterface::~qpOASESInterface()
{
}

/**
 * @brief Allocate memory for the class members
 * @param nlp_info  the struct that stores simple nlp dimension info
 * @param qptype is the problem to be solved QP or LP or SOC?
 * @return
 */
void qpOASESInterface::allocate_memory(
    std::shared_ptr<const SqpNlpSizeInfo> nlp_sizes, QPType qptype)
{

#ifdef NEW_FORMULATION
  int nnz_g_QP = nlp_sizes->get_num_nonzeros_jacobian() +
                 2 * nlp_sizes->get_num_constraints() +
                 3 * nlp_sizes->get_num_variables();
#else
  int nnz_g_QP = nlp_sizes->get_num_nonzeros_jacobian() +
                 2 * nlp_sizes->get_num_constraints();
#endif
  lbA_ = make_shared<Vector>(num_qp_constraints_);
  ubA_ = make_shared<Vector>(num_qp_constraints_);
  lb_ = make_shared<Vector>(num_qp_variables_);
  ub_ = make_shared<Vector>(num_qp_variables_);
  g_ = make_shared<Vector>(num_qp_variables_);
  A_ = make_shared<SparseHbMatrix>(nnz_g_QP, num_qp_constraints_,
                                   num_qp_variables_, false);
  x_qp_ = make_shared<Vector>(num_qp_variables_);
  y_qp_ = make_shared<Vector>(num_qp_constraints_ + num_qp_variables_);

  if (qptype != LP) {
    H_ = make_shared<SparseHbMatrix>(num_qp_variables_, num_qp_variables_,
                                     false);
  }

  qpoases_solver_ = std::make_shared<qpOASES::SQProblem>(
      (qpOASES::int_t)num_qp_variables_, (qpOASES::int_t)num_qp_constraints_);
}

QpSolverExitStatus qpOASESInterface::get_qpoases_exit_status_(
    shared_ptr<qpOASES::SQProblem> qpoases_solver)
{
  QpSolverExitStatus retval = QPEXIT_UNKNOWN_STATUS;
  if (qpoases_solver->isSolved()) {
    retval = QPEXIT_OPTIMAL;
  } else if (qpoases_solver->isInfeasible()) {
    retval = QPEXIT_INFEASIBLE;
  } else if (qpoases_solver->isUnbounded()) {
    retval = QPEXIT_UNBOUNDED;
  }

  return retval;
}

/**
 * @brief This method solves the QP problem specified in the data, with given
 * options.
 * After the QP being solved, it updates the stats, adding the iteration
 * number used to solve the QP to the qp_iter in object stats
 */
QpSolverExitStatus qpOASESInterface::optimize_qp(shared_ptr<Statistics> stats)
{
  QpSolverExitStatus retval = QPEXIT_NOT_SOLVED;

  qpOASES::int_t nWSR = qp_solver_max_num_iterations_; // TODO modify it

  if (!first_qp_solved_) { // if haven't solve any QP before then initialize the
                           // first QP
    set_qp_solver_options_();

#if 0
    //@{
    // for debugging
            H_qpOASES_->print("H_qp_oases");
            A_qpOASES_->print("A_qpoases");
            g_->print("g");
            lbA_->print("LbA");
            ubA_->print("ubA");
            lb_->print("lb");
            ub_->print("ub");
    //@}
#endif

    qpoases_solver_->init(H_qpOASES_.get(), g_->get_values(), A_qpOASES_.get(),
                          lb_->get_values(), ub_->get_values(),
                          lbA_->get_values(), ubA_->get_values(), nWSR);

    retval = get_qpoases_exit_status_(qpoases_solver_);

    if (retval == QPEXIT_OPTIMAL) {
      first_qp_solved_ = true;
    }

    // For now, don't try to fix any QP solver error
    // handle_error(QP, stats);
  } else {
    // for debugging
    //@{
    //                  H_qpOASES_->print("H_qp_oases");
    //                  A_qpOASES_->print("A_qpoases");
    //                  g_->print("g");
    //                  lbA_->print("LbA");
    //                  ubA_->print("ubA");
    //                  lb_->print("lb");
    //                  ub_->print("ub");
    //@}
    get_Matrix_change_status();
    if (new_QP_matrix_status_ == UNDEFINED) {
      assert(old_QP_matrix_status_ != UNDEFINED);
      if (old_QP_matrix_status_ == FIXED)
        qpoases_solver_->hotstart(g_->get_values(), lb_->get_values(),
                                  ub_->get_values(), lbA_->get_values(),
                                  ubA_->get_values(), nWSR);
      else {
        qpoases_solver_->hotstart(H_qpOASES_.get(), g_->get_values(),
                                  A_qpOASES_.get(), lb_->get_values(),
                                  ub_->get_values(), lbA_->get_values(),
                                  ubA_->get_values(), nWSR);
      }
    } else {
      if (new_QP_matrix_status_ == FIXED && old_QP_matrix_status_ == FIXED) {
        qpoases_solver_->hotstart(g_->get_values(), lb_->get_values(),
                                  ub_->get_values(), lbA_->get_values(),
                                  ubA_->get_values(), nWSR);
      } else if (new_QP_matrix_status_ == VARIED &&
                 old_QP_matrix_status_ == VARIED) {
        qpoases_solver_->hotstart(H_qpOASES_.get(), g_->get_values(),
                                  A_qpOASES_.get(), lb_->get_values(),
                                  ub_->get_values(), lbA_->get_values(),
                                  ubA_->get_values(), nWSR);
      } else if (new_QP_matrix_status_ != old_QP_matrix_status_) {
        qpOASES::Bounds tmp_bounds;
        qpoases_solver_->getBounds(tmp_bounds);
        qpoases_solver_->init(H_qpOASES_.get(), g_->get_values(),
                              A_qpOASES_.get(), lb_->get_values(),
                              ub_->get_values(), lbA_->get_values(),
                              ubA_->get_values(), nWSR, 0, x_qp_->get_values(),
                              y_qp_->get_values(), &tmp_bounds);
        new_QP_matrix_status_ = old_QP_matrix_status_ = UNDEFINED;
      }
    }
    retval = get_qpoases_exit_status_(qpoases_solver_);
  }

  if (retval != QPEXIT_OPTIMAL) {
    return retval;
  }

  // reset the flag that remembers whether QP matrices have changed
  qp_matrices_changed_ = false;

  if (stats != nullptr)
    stats->increase_qp_iteration_counter((int)nWSR);

  // For now not trying to fix things if qpOASES fails
  if (!qpoases_solver_->isSolved()) {
    assert(false && "Not trying to fix qpOASES solution yet");
    handle_error(QP, stats);
  }

  // Retrieve the primal and dual solution from qpOASES
  qpoases_solver_->getPrimalSolution(x_qp_->get_non_const_values());
  qpoases_solver_->getDualSolution(y_qp_->get_non_const_values());

  return retval;
}

QpSolverExitStatus qpOASESInterface::optimize_lp(shared_ptr<Statistics> stats)
{
  QpSolverExitStatus retval = QPEXIT_NOT_SOLVED;

  qpOASES::int_t nWSR = lp_solver_max_num_iterations_; // TODO modify it

  // If this is the first solve, we need to call a different method
  if (!first_qp_solved_) {
    set_qp_solver_options_();
    qpoases_solver_->init(0, g_->get_values(), A_qpOASES_.get(),
                          lb_->get_values(), ub_->get_values(),
                          lbA_->get_values(), ubA_->get_values(), nWSR);
    retval = get_qpoases_exit_status_(qpoases_solver_);

    if (retval == QPEXIT_OPTIMAL) {
      first_qp_solved_ = true;
    } else {
      assert("Not trying to fix qpOASES solutions yet" && false);
      handle_error(LP, stats);
    }
  } else {
    get_Matrix_change_status();
    if (new_QP_matrix_status_ == UNDEFINED) {
      if (old_QP_matrix_status_ == FIXED)
        qpoases_solver_->hotstart(g_->get_values(), lb_->get_values(),
                                  ub_->get_values(), lbA_->get_values(),
                                  ubA_->get_values(), nWSR);
      else {
        qpoases_solver_->hotstart(0, g_->get_values(), A_qpOASES_.get(),
                                  lb_->get_values(), ub_->get_values(),
                                  lbA_->get_values(), ubA_->get_values(), nWSR);
      }
    }

    else {
      if (new_QP_matrix_status_ == FIXED && old_QP_matrix_status_ == FIXED) {
        qpoases_solver_->hotstart(g_->get_values(), lb_->get_values(),
                                  ub_->get_values(), lbA_->get_values(),
                                  ubA_->get_values(), nWSR);
      } else if (new_QP_matrix_status_ == VARIED &&
                 old_QP_matrix_status_ == VARIED) {
        qpoases_solver_->hotstart(0, g_->get_values(), A_qpOASES_.get(),
                                  lb_->get_values(), ub_->get_values(),
                                  lbA_->get_values(), ubA_->get_values(), nWSR);
      } else if (new_QP_matrix_status_ != old_QP_matrix_status_) {
        qpoases_solver_->init(0, g_->get_values(), A_qpOASES_.get(),
                              lb_->get_values(), ub_->get_values(),
                              lbA_->get_values(), ubA_->get_values(), nWSR);
        new_QP_matrix_status_ = old_QP_matrix_status_ = UNDEFINED;
      }
      retval = get_qpoases_exit_status_(qpoases_solver_);
    }

    // reset the flag that remembers whether QP matrices have changed
    qp_matrices_changed_ = false;

    if (!qpoases_solver_->isSolved()) {
      assert("Not trying to fix qpOASES solutions yet" && false);
      handle_error(LP, stats);
    }
  }
  // get primal and dual solutions
  if (stats != nullptr)
    stats->increase_qp_iteration_counter((int)nWSR);

  // Retrieve the primal and dual solution from qpOASES
  qpoases_solver_->getPrimalSolution(x_qp_->get_non_const_values());
  qpoases_solver_->getDualSolution(y_qp_->get_non_const_values());

  return retval;
}

shared_ptr<const Vector> qpOASESInterface::get_primal_solution() const
{
  shared_ptr<const Vector> retval = make_shared<Vector>(*x_qp_);
  return retval;
}

/**
 * @brief get the pointer to the multipliers to the bounds constraints.
 */
shared_ptr<const Vector> qpOASESInterface::get_bounds_multipliers() const
{
  // create a new vector with the data
  shared_ptr<Vector> retval = make_shared<Vector>(num_nlp_variables_);

// copy the values from the beginning of the y_qp_ vector
#ifdef NEW_FORMULATION
  // The qpOASES multipliers are organized as follows:
  // bounds for p, u and v slack variables, multpliers for x bounds, multipliers
  // for constraints
  //  retval->copy_values(y_qp_->get_values() + 3 * num_nlp_constraints_);
  retval->copy_values(y_qp_->get_values() + num_nlp_variables_ +
                      2 * num_nlp_constraints_);
#else
  retval->copy_values(y_qp_->get_values());
#endif

  return retval;
}

/**
 * @brief get the pointer to the multipliers to the regular constraints.
 */
shared_ptr<const Vector> qpOASESInterface::get_constraints_multipliers() const
{
  // create a new vector with the data
  shared_ptr<Vector> retval = make_shared<Vector>(num_nlp_constraints_);

  // copy the value from the second part of the y_qp_ vector
  retval->copy_values(y_qp_->get_values() + num_qp_variables_);

  return retval;
}

/**
 *@brief get the objective value from the QP solvers
 *
 * @return the objective function value of the QP problem
 */

inline double qpOASESInterface::get_obj_value()
{

  return (double)(qpoases_solver_->getObjVal());
}

Exitflag qpOASESInterface::get_status()
{
  qpOASES::QProblemStatus finalStatus = qpoases_solver_->getStatus();

  if (qpoases_solver_->isInfeasible()) {
    return QPERROR_INFEASIBLE;
  } else if (qpoases_solver_->isUnbounded()) {
    return QPERROR_UNBOUNDED;
  } else if (qpoases_solver_->isSolved()) {
    return QP_OPTIMAL;
  } else
    switch (finalStatus) {
      case qpOASES::QPS_NOTINITIALISED:
        return QPERROR_NOTINITIALISED;
      case qpOASES::QPS_PREPARINGAUXILIARYQP:
        return QPERROR_PREPARINGAUXILIARYQP;
      case qpOASES::QPS_AUXILIARYQPSOLVED:
        return QPERROR_AUXILIARYQPSOLVED;
      case qpOASES::QPS_PERFORMINGHOMOTOPY:
        return QPERROR_PERFORMINGHOMOTOPY;
      case qpOASES::QPS_HOMOTOPYQPSOLVED:
        return QPERROR_HOMOTOPYQPSOLVED;
    }
}

//@}
void qpOASESInterface::set_lb(int location, double value)
{
  lb_->set_value(location, value);
}

void qpOASESInterface::set_ub(int location, double value)
{
  ub_->set_value(location, value);
}

void qpOASESInterface::set_lbA(int location, double value)
{
  lbA_->set_value(location, value);
}

void qpOASESInterface::set_ubA(int location, double value)
{
  ubA_->set_value(location, value);
}

void qpOASESInterface::set_gradient(int location, double value)
{
  g_->set_value(location, value);
}

void qpOASESInterface::set_hessian(shared_ptr<const SpTripletMat> rhs)
{
  // Remember that a change is made to the Hessian matrix
  qp_matrices_changed_ = true;

  if (!H_->is_initialized()) {
    H_->set_structure(rhs); // TODO: move to somewhere else?
    H_->set_values(rhs);
    int* H_const_row_indices = const_cast<int*>(H_->get_row_indices());
    int* H_const_col_indices = const_cast<int*>(H_->get_column_indices());
    double* H_const_values = const_cast<double*>(H_->get_values());
    H_qpOASES_ = std::make_shared<qpOASES::SymSparseMat>(
        num_qp_variables_, num_qp_variables_, H_const_row_indices,
        H_const_col_indices, H_const_values);
  } else {
    H_->set_values(rhs);
    H_qpOASES_->setVal(H_->get_values());
  }
  H_qpOASES_->createDiagInfo();
}

void qpOASESInterface::set_jacobian(
    shared_ptr<const SQPhotstart::SpTripletMat> rhs,
    IdentityMatrixPositions& identity_matrix_positions)
{
  // Remember that a change is made to the Jacobian matrix
  qp_matrices_changed_ = true;

  if (!A_->is_initialized()) {
    A_->set_structure(rhs, identity_matrix_positions);
    A_->set_values(rhs);
    int* A_const_row_indices = const_cast<int*>(A_->get_row_indices());
    int* A_const_col_indices = const_cast<int*>(A_->get_column_indices());
    double* A_const_values = const_cast<double*>(A_->get_values());
    A_qpOASES_ = std::make_shared<qpOASES::SparseMatrix>(
        num_qp_constraints_, num_qp_variables_, A_const_row_indices,
        A_const_col_indices, A_const_values);
  } else {
    A_->set_values(rhs);
    A_qpOASES_->setVal(A_->get_values());
  }
}

void qpOASESInterface::set_ub(shared_ptr<const Vector> rhs)
{
  ub_->copy_vector(rhs);
}

void qpOASESInterface::set_lb(shared_ptr<const Vector> rhs)
{
  lb_->copy_vector(rhs);
}

void qpOASESInterface::set_lbA(shared_ptr<const Vector> rhs)
{
  lbA_->copy_vector(rhs);
}

void qpOASESInterface::set_ubA(shared_ptr<const Vector> rhs)
{
  ubA_->copy_vector(rhs);
}

void qpOASESInterface::set_gradient(shared_ptr<const Vector> rhs)
{
  g_->copy_vector(rhs);
}

bool qpOASESInterface::test_optimality(ActivityStatus* W_c, ActivityStatus* W_b)
{

  int i;
  // create local variables and set all violation values to be 0
  double primal_violation = 0.0;
  double dual_violation = 0.0;
  double compl_violation = 0.0;
  double statioanrity_violation = 0.0;
  shared_ptr<Vector> Ax = make_shared<Vector>(num_qp_constraints_);

  if (W_c == NULL && W_b == NULL) {
    W_c = new ActivityStatus[num_qp_constraints_];
    W_b = new ActivityStatus[num_qp_variables_];
  }
  get_working_set(W_c, W_b);

  /**-------------------------------------------------------**/
  /**                    primal feasibility                 **/
  /**-------------------------------------------------------**/
  for (i = 0; i < num_qp_variables_; i++) {
    primal_violation += max(0.0, (lb_->get_value(i) - x_qp_->get_value(i)));
    primal_violation += -min(0.0, (ub_->get_value(i) - x_qp_->get_value(i)));
  }
  if (A_ != nullptr) {
    A_->multiply(x_qp_, Ax); // tmp_vec_nCon=A*x
    for (i = 0; i < num_qp_constraints_; i++) {
      primal_violation += max(0.0, (lbA_->get_value(i) - Ax->get_value(i)));
      primal_violation += -min(0.0, (ubA_->get_value(i) - Ax->get_value(i)));
    }
  }

  /**-------------------------------------------------------**/
  /**                    dual feasibility                   **/
  /**-------------------------------------------------------**/
  for (i = 0; i < num_qp_variables_; i++) {
    switch (W_b[i]) {
      case INACTIVE: // the constraint is inactive, then the dual multiplier
        // should be 0
        dual_violation += fabs(y_qp_->get_value(i));
        break;
      case ACTIVE_BELOW: // the constraint is active at the lower bound, so the
        // multiplier should be positive
        dual_violation += -min(0.0, y_qp_->get_value(i));
        break;
      case ACTIVE_ABOVE: // the contraint is active at the upper bounds, so the
        // multiplier should be negavie
        dual_violation += max(0.0, y_qp_->get_value(i));
        break;
      case ACTIVE_BOTH_SIDES:
        break;
      default:
        assert("Invalid working set flag" && false);
    }
  }
  if (A_ != nullptr) {
    for (i = 0; i < num_qp_constraints_; i++) {
      switch (W_c[i]) {
        case INACTIVE: // the constraint is inactive, then the dual multiplier
          // should be 0
          dual_violation += fabs(y_qp_->get_value(i + num_qp_variables_));
          break;
        case ACTIVE_BELOW: // the constraint is active at the lower bound, so
                           // the
          // multiplier should be positive
          dual_violation += -min(0.0, y_qp_->get_value(i + num_qp_variables_));
          break;
        case ACTIVE_ABOVE: // the contraint is active at the upper bounds, so
                           // the
          // multiplier should be negavie
          dual_violation += max(0.0, y_qp_->get_value(i + num_qp_variables_));
          break;
        case ACTIVE_BOTH_SIDES:
          break;
        default:
          assert("Invalid working set flag" && false);
      }
    }
  }

  /**-------------------------------------------------------**/
  /**                   stationarity                        **/
  /**-------------------------------------------------------**/
  // calculate A'*y+lambda-(g+Hx)

  shared_ptr<Vector> stationary_gap = make_shared<Vector>(num_qp_variables_);
  //    A_->print("A");
  //    H_->print("H");
  //    x_qp_->print("x_qp");
  //    lb_->print("lb");
  //    ub_->print("ub");
  //    y_qp_->print("y_qp");

  if (A_ != nullptr) {
    A_->multiply_transpose(y_qp_->get_values() + num_qp_variables_,
                           stationary_gap->get_values());
  }
  shared_ptr<Vector> Hx = make_shared<Vector>(num_qp_variables_);
  H_->multiply(x_qp_, Hx);

  stationary_gap->add_elements(1., y_qp_->get_values());
  stationary_gap->add_vector(-1., g_);
  stationary_gap->add_vector(-1., Hx);
  statioanrity_violation = stationary_gap->calc_one_norm();

  /**-------------------------------------------------------**/
  /**                    Complemtarity                      **/
  /**-------------------------------------------------------**/

  for (i = 0; i < num_qp_variables_; i++) {
    switch (W_b[i]) {
      case INACTIVE: // constraint is inactive, multiplier should be 0
        compl_violation += abs(y_qp_->get_value(i));
        break;
      case ACTIVE_BELOW: // the constraint is active at the lower bound
        compl_violation += abs(y_qp_->get_value(i) *
                               (x_qp_->get_value(i) - lb_->get_value(i)));
        break;
      case ACTIVE_ABOVE: // the contraint is active at the upper bounds, so the
        // multiplier should be negavie
        compl_violation += abs(y_qp_->get_value(i) *
                               (ub_->get_value(i) - x_qp_->get_value(i)));
        break;
      case ACTIVE_BOTH_SIDES:
        break;
      default:
        assert("Invalid working set flag" && false);
    }
  }
  if (A_ != nullptr) {
    for (i = 0; i < num_qp_constraints_; i++) {
      switch (W_c[i]) {
        case INACTIVE: // constraint is inactive, multiplier should be 0
          compl_violation += abs(y_qp_->get_value(i + num_qp_variables_));
          break;
        case ACTIVE_BELOW: // the constraint is active at the lower bound
          compl_violation += abs(y_qp_->get_value(i + num_qp_variables_) *
                                 (Ax->get_value(i) - lbA_->get_value(i)));
          break;
        case ACTIVE_ABOVE: // the contraint is active at the upper bounds, so
                           // the
          // multiplier should be negavie
          compl_violation += abs(y_qp_->get_value(i + num_qp_variables_) *
                                 (ubA_->get_value(i) - Ax->get_value(i)));
          break;
        case ACTIVE_BOTH_SIDES:
          break;
        default:
          assert("Invalid working set flag" && false);
      }
    }
  }

  qpOptimalStatus_.compl_violation = compl_violation;
  qpOptimalStatus_.stationarity_violation = statioanrity_violation;
  qpOptimalStatus_.dual_violation = dual_violation;
  qpOptimalStatus_.primal_violation = primal_violation;
  qpOptimalStatus_.KKT_error = compl_violation + statioanrity_violation +
                               dual_violation + primal_violation;

  if (W_c == NULL && W_b == NULL) {
    delete[] W_c;
    delete[] W_b;
  }

  if (qpOptimalStatus_.KKT_error > 1.0e-6) {
    //        printf("comp_violation %10e\n", compl_violation);
    //        printf("stat_violation %10e\n", statioanrity_violation);
    //        printf("prim_violation %10e\n", primal_violation);
    //        printf("dual_violation %10e\n", dual_violation);
    //        printf("KKT_error %10e\n", qpOptimalStatus_.KKT_error);
    return false;
  }

  return true;
}

void qpOASESInterface::handle_error(QPType qptype, shared_ptr<Statistics> stats)
{
  assert("Not trying to fix qpOASES solutions yet" && false);
  if (qptype == LP) {
    qpOASES::int_t nWSR = lp_solver_max_num_iterations_; // TODO modify it
    if (qpoases_solver_->isInfeasible()) {
      shared_ptr<Vector> x_0 = make_shared<Vector>(num_qp_variables_);
      shared_ptr<Vector> Ax = make_shared<Vector>(num_qp_constraints_);
      x_0->copy_vector(x_qp_);
      A_->multiply(x_0, Ax);

      for (int i = 0; i < num_qp_constraints_; i++) {
        x_0->set_value(i + num_qp_variables_ - 2 * num_qp_constraints_,
                       max(0.0, lbA_->get_value(i)));
        x_0->set_value(i + num_qp_variables_ - num_qp_constraints_,
                       -min(0.0, ubA_->get_value(i)));
      }
      qpoases_solver_->init(H_qpOASES_.get(), g_->get_values(),
                            A_qpOASES_.get(), lb_->get_values(),
                            ub_->get_values(), lbA_->get_values(),
                            ubA_->get_values(), nWSR, NULL, x_0->get_values());
    } else {
      qpoases_solver_->init(0, g_->get_values(), A_qpOASES_.get(),
                            lb_->get_values(), ub_->get_values(),
                            lbA_->get_values(), ubA_->get_values(), nWSR);
    }
    old_QP_matrix_status_ = new_QP_matrix_status_ = UNDEFINED;
    if (stats != nullptr)
      stats->increase_qp_iteration_counter((int)nWSR);

    if (!qpoases_solver_->isSolved()) {
      //  THROW_EXCEPTION(LP_NOT_OPTIMAL, LP_NOT_OPTIMAL_MSG);
    }
  } else {
    qpOASES::int_t nWSR = qp_solver_max_num_iterations_; // TODO modify it
    if (qpoases_solver_->isInfeasible()) {
      shared_ptr<Vector> x_0 = make_shared<Vector>(num_qp_variables_);
      for (int i = 0; i < num_qp_constraints_; i++) {
        x_0->set_value(i + num_qp_variables_ - 2 * num_qp_constraints_,
                       max(0.0, lbA_->get_value(i)));
        x_0->set_value(i + num_qp_variables_ - num_qp_constraints_,
                       -min(0.0, ubA_->get_value(i)));
      }
      qpoases_solver_->init(H_qpOASES_.get(), g_->get_values(),
                            A_qpOASES_.get(), lb_->get_values(),
                            ub_->get_values(), lbA_->get_values(),
                            ubA_->get_values(), nWSR, NULL, x_0->get_values());
      // for debugging
      //@{
      //            H_qpOASES_->print("H_qp_oases");
      //            A_qpOASES_->print("A_qpoases");
      //            g_->print("g");
      //            lbA_->print("LbA");
      //            ubA_->print("ubA");
      //            lb_->print("lb");
      //            ub_->print("ub");
      //            x_0->print("x_0");
      //            shared_ptr<Vector> Ax=
      //            make_shared<Vector>(num_qp_constraints_);
      //
      //            A_->times(x_0, Ax);
      //            Ax->print("Ax");

      //@}
    } else {
      qpoases_solver_->init(H_qpOASES_.get(), g_->get_values(),
                            A_qpOASES_.get(), lb_->get_values(),
                            ub_->get_values(), lbA_->get_values(),
                            ubA_->get_values(), nWSR);
    }
    old_QP_matrix_status_ = new_QP_matrix_status_ = UNDEFINED;
    if (stats != nullptr)
      stats->increase_qp_iteration_counter((int)nWSR);
    if (!qpoases_solver_->isSolved()) {
      // THROW_EXCEPTION(QP_NOT_OPTIMAL, QP_NOT_OPTIMAL_MSG);
    }
  }
}

void qpOASESInterface::set_qp_solver_options_()
{
  qpOASES::Options qp_options;

  qp_options.setToReliable();
  // setup the printlevel of q
  switch (qp_solver_print_level_) {
    case 0:
      qp_options.printLevel = qpOASES::PL_NONE;
      break;
    case 1:
      qp_options.printLevel = qpOASES::PL_TABULAR;
      break;
    case 2:
      qp_options.printLevel = qpOASES::PL_LOW;
      break;
    case 3:
      qp_options.printLevel = qpOASES::PL_MEDIUM;
      break;
    case 4:
      qp_options.printLevel = qpOASES::PL_HIGH;
      break;
    case -2:
      qp_options.printLevel = qpOASES::PL_DEBUG_ITER;
      break;
  }
  qpoases_solver_->setOptions(qp_options);
}

void qpOASESInterface::WriteQPDataToFile(EJournalLevel level,
                                         EJournalCategory category,
                                         const string filename)
{
#ifdef DEBUG
#ifdef PRINT_OUT_QP_WITH_ERROR
  jnlst_->DeleteAllJournals();
  SmartPtr<Journal> QPdata_jrnl =
      jnlst_->AddFileJournal("QPdata", "qpOASES" + filename, J_WARNING);
  QPdata_jrnl->SetAllPrintLevels(level);
  QPdata_jrnl->SetPrintLevel(category, level);

  lb_->write_to_file("lb", jnlst_, level, category, QPOASES);
  lbA_->write_to_file("lbA", jnlst_, level, category, QPOASES);
  ub_->write_to_file("ub", jnlst_, level, category, QPOASES);
  ubA_->write_to_file("ubA", jnlst_, level, category, QPOASES);

  g_->write_to_file("g", jnlst_, level, category, QPOASES);
  A_->write_to_file("A", jnlst_, level, category, QPOASES);
  H_->write_to_file("H", jnlst_, level, category, QPOASES);
  jnlst_->DeleteAllJournals();
#endif
#endif
}

void qpOASESInterface::get_Matrix_change_status()
{

  if (old_QP_matrix_status_ == UNDEFINED) {
    old_QP_matrix_status_ = qp_matrices_changed_ ? VARIED : FIXED;
  } else {
    if (new_QP_matrix_status_ != UNDEFINED)
      old_QP_matrix_status_ = new_QP_matrix_status_;
    new_QP_matrix_status_ = qp_matrices_changed_ ? VARIED : FIXED;
  }
}

void qpOASESInterface::get_working_set(SQPhotstart::ActivityStatus* W_constr,
                                       SQPhotstart::ActivityStatus* W_bounds)
{
  int* tmp_W_c = new int[num_qp_constraints_];
  int* tmp_W_b = new int[num_qp_variables_];

  assert(num_qp_constraints_ == qpoases_solver_->getNC());
  assert(num_qp_variables_ == qpoases_solver_->getNV());
  qpoases_solver_->getWorkingSetConstraints(tmp_W_c);
  qpoases_solver_->getWorkingSetBounds(tmp_W_b);

  for (int i = 0; i < num_qp_variables_; i++) {
    switch ((int)tmp_W_b[i]) {
      case 1:
        if (fabs(x_qp_->get_value(i) - lb_->get_value(i)) < sqrt_m_eps)
          W_bounds[i] = ACTIVE_BOTH_SIDES;
        else
          W_bounds[i] = ACTIVE_ABOVE;
        break;
      case -1:
        if (fabs(x_qp_->get_value(i) - ub_->get_value(i)) < sqrt_m_eps)
          W_bounds[i] = ACTIVE_BOTH_SIDES;
        else
          W_bounds[i] = ACTIVE_BELOW;

        break;
      case 0:
        W_bounds[i] = INACTIVE;
        break;
      default:
        assert("Invalid working set flag" && false);
    }
  }
  auto Ax = make_shared<Vector>(num_qp_constraints_);
  A_->multiply(x_qp_, Ax); // tmp_vec_nCon=A*x
  for (int i = 0; i < num_qp_constraints_; i++) {
    switch ((int)tmp_W_c[i]) {
      case 1:
        if (fabs(Ax->get_value(i) - lbA_->get_value(i) < sqrt_m_eps))
          W_constr[i] = ACTIVE_BOTH_SIDES;
        else
          W_constr[i] = ACTIVE_ABOVE;
        break;
      case -1:
        if (fabs(Ax->get_value(i) - ubA_->get_value(i) < sqrt_m_eps))
          W_constr[i] = ACTIVE_BOTH_SIDES;
        else
          W_constr[i] = ACTIVE_BELOW;
        break;
      case 0:
        W_constr[i] = INACTIVE;
        break;
      default:
        assert("Invalid working set flag" && false);
    }
  }
  delete[] tmp_W_b;
  delete[] tmp_W_c;
}

void qpOASESInterface::reset_constraints()
{
  lb_->set_to_zero();
  ub_->set_to_zero();
  lbA_->set_to_zero();
  ubA_->set_to_zero();
}

} // SQPHOTSTART
