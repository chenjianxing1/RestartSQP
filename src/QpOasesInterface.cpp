/* Copyright (C) 2019
 * All Rights Reserved.
 *
 * Authors: Xinyi Luo
 * Date:    2019-07
 */
#include "sqphot/QpOasesInterface.hpp"
#include "sqphot/MessageHandling.hpp"

using namespace std;
using namespace Ipopt;

namespace RestartSqp {

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
 * @param qp_type  is the problem to be solved QP or LP or SOC?
 */
QpOasesInterface::QpOasesInterface(int num_qp_variables, int num_qp_constraints,
                                   QPType qp_type,
                                   SmartPtr<const OptionsList> options,
                                   SmartPtr<Journalist> jnlst)
 : QpSolverInterface(qp_type, num_qp_variables, num_qp_constraints, jnlst)
 , first_qp_solved_(false)
 , qp_matrices_changed_(true)
{
  // Get the option values for this object
  get_option_values_(options);

  // Create the qpOASES solver object
  qpoases_solver_ = make_shared<qpOASES::SQProblem>(
      (qpOASES::int_t)num_qp_variables_, (qpOASES::int_t)num_qp_constraints_);

  // Set the options for qpOASES
  set_qp_solver_options_();

  // Create matrix objects (without data)
  if (num_qp_constraints_ > 0) {
    bool is_compressed_row = false;
    bool is_symmetric = false;
    jacobian_ =
        make_shared<SparseHbMatrix>(num_qp_constraints_, num_qp_variables_,
                                    is_compressed_row, is_symmetric);
  }

  if (qp_type_ == QP_TYPE_QP) {
    bool is_compressed_row = false;
    bool is_symmetric = true;
    hessian_ = make_shared<SparseHbMatrix>(num_qp_variables_, num_qp_variables_,
                                           is_compressed_row, is_symmetric);
  }
}

void QpOasesInterface::get_option_values_(SmartPtr<const OptionsList> options)
{
  // Get the options from the options list
  options->GetIntegerValue("qp_solver_max_num_iterations",
                           qp_solver_max_num_iterations_, "");
  options->GetIntegerValue("qp_solver_print_level", qp_solver_print_level_, "");
  options->GetIntegerValue("lp_solver_max_num_iterations",
                           lp_solver_max_num_iterations_, "");
}

#if 0
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
 , hessian_(H)
 , jacobian_(A)
 , g_(g)
 , lb_(lb)
 , ub_(ub)
 , lbA_(lbA)
 , ubA_(ubA)
{
  jacobian_->print("bla");
  // Get the option values for this object
  get_option_values_(options);

  x_qp_ = make_shared<Vector>(num_qp_variables_);
  y_qp_ = make_shared<Vector>(num_qp_constraints_ + num_qp_variables_);
  qpoases_solver_ = make_shared<qpOASES::SQProblem>(
      (qpOASES::int_t)num_qp_variables_, (qpOASES::int_t)num_qp_constraints_);
  // AW: The following is very bad design in qpOASES: the constructor does not
  // copy
  // the values, but the only time the values
  // are actually changed is if the setVal method is called in
  // qpOASES::SparseMatrix.  It should be const.
  int* H_const_row_indices = const_cast<int*>(hessian_->get_row_indices());
  int* H_const_col_indices = const_cast<int*>(hessian_->get_column_indices());
  double* H_const_values = const_cast<double*>(hessian_->get_values());
  qpoases_hessian_ = make_shared<qpOASES::SymSparseMat>(
      num_qp_variables_, num_qp_variables_, H_const_row_indices,
      H_const_col_indices, H_const_values);
  qpoases_hessian_->createDiagInfo();

  int* A_const_row_indices = const_cast<int*>(jacobian_->get_row_indices());
  int* A_const_col_indices = const_cast<int*>(jacobian_->get_column_indices());
  double* A_const_values = const_cast<double*>(jacobian_->get_values());
  qpoases_jacobian_ = make_shared<qpOASES::SymSparseMat>(
      num_qp_constraints_, num_qp_variables_, A_const_row_indices,
      A_const_col_indices, A_const_values);
}
#endif

/**Default destructor*/
QpOasesInterface::~QpOasesInterface()
{
}

QpSolverExitStatus QpOasesInterface::get_qpoases_exit_status_()
{
  QpSolverExitStatus retval = QPEXIT_UNKNOWN_STATUS;
  if (qpoases_solver_->isSolved()) {
    retval = QPEXIT_OPTIMAL;
  } else if (qpoases_solver_->isInfeasible()) {
    retval = QPEXIT_INFEASIBLE;
  } else if (qpoases_solver_->isUnbounded()) {
    retval = QPEXIT_UNBOUNDED;
  } else {
    qpOASES::QProblemStatus finalStatus = qpoases_solver_->getStatus();
    retval = QPEXIT_INTERNAL_ERROR;
    jnlst_->Printf(J_ERROR, J_MAIN, "qpOASES Error with finalStatus = %d ",
                   (int)finalStatus);
    string error_code;
    switch (finalStatus) {
      case qpOASES::QPS_NOTINITIALISED:
        error_code = "QPS_NOTINITIALISED";
        break;
      case qpOASES::QPS_PREPARINGAUXILIARYQP:
        error_code = "QPS_PREPARINGAUXILIARYQP";
        break;
      case qpOASES::QPS_AUXILIARYQPSOLVED:
        error_code = "QPS_AUXILIARYQPSOLVED";
        break;
      case qpOASES::QPS_PERFORMINGHOMOTOPY:
        error_code = "QPS_PERFORMINGHOMOTOPY";
        break;
      case qpOASES::QPS_HOMOTOPYQPSOLVED:
        error_code = "QPS_HOMOTOPYQPSOLVED";
        break;
    }
    jnlst_->Printf(J_ERROR, J_MAIN, "(%s)\n", error_code.c_str());
  }

  return retval;
}

/**
 * @brief This method solves the QP problem specified in the data, with given
 * options.
 * After the QP being solved, it updates the stats, adding the iteration
 * number used to solve the QP to the qp_iter in object stats
 */
bool QpOasesInterface::optimize_impl(shared_ptr<Statistics> stats)
{
  bool solve_successful = false;

  qpOASES::int_t num_qp_iterations =
      qp_solver_max_num_iterations_; // TODO modify it

  if (!first_qp_solved_) {
    // This is the first call and we need to use the initial init method
    qpoases_solver_->init(
        qpoases_hessian_.get(), linear_objective_coefficients_->get_values(),
        qpoases_jacobian_.get(), lower_variable_bounds_->get_values(),
        upper_variable_bounds_->get_values(),
        lower_constraint_bounds_->get_values(),
        upper_constraint_bounds_->get_values(), num_qp_iterations);

  } else {
    // We already solved one QP with qpOASES
    if (qp_matrices_changed_) {
      qpoases_solver_->hotstart(
          qpoases_hessian_.get(), linear_objective_coefficients_->get_values(),
          qpoases_jacobian_.get(), lower_variable_bounds_->get_values(),
          upper_variable_bounds_->get_values(),
          lower_constraint_bounds_->get_values(),
          upper_constraint_bounds_->get_values(), num_qp_iterations);
    } else {
      qpoases_solver_->hotstart(linear_objective_coefficients_->get_values(),
                                lower_variable_bounds_->get_values(),
                                upper_variable_bounds_->get_values(),
                                lower_constraint_bounds_->get_values(),
                                upper_constraint_bounds_->get_values(),
                                num_qp_iterations);
    }
#if 0
    get_Matrix_change_status();
    if (new_QP_matrix_status_ == UNDEFINED) {
      assert(old_QP_matrix_status_ != UNDEFINED);
      if (old_QP_matrix_status_ == FIXED)
        qpoases_solver_->hotstart(g_->get_values(), lb_->get_values(),
                                  ub_->get_values(), lbA_->get_values(),
                                  ubA_->get_values(), num_qp_iterations);
      else {
        qpoases_solver_->hotstart(qpoases_hessian_.get(), g_->get_values(),
                                  qpoases_jacobian_.get(), lb_->get_values(),
                                  ub_->get_values(), lbA_->get_values(),
                                  ubA_->get_values(), num_qp_iterations);
      }
    } else {
      if (new_QP_matrix_status_ == FIXED && old_QP_matrix_status_ == FIXED) {
        qpoases_solver_->hotstart(g_->get_values(), lb_->get_values(),
                                  ub_->get_values(), lbA_->get_values(),
                                  ubA_->get_values(), num_qp_iterations);
      } else if (new_QP_matrix_status_ == VARIED &&
                 old_QP_matrix_status_ == VARIED) {
        qpoases_solver_->hotstart(qpoases_hessian_.get(), g_->get_values(),
                                  qpoases_jacobian_.get(), lb_->get_values(),
                                  ub_->get_values(), lbA_->get_values(),
                                  ubA_->get_values(), num_qp_iterations);
      } else if (new_QP_matrix_status_ != old_QP_matrix_status_) {
        qpOASES::Bounds tmp_bounds;
        qpoases_solver_->getBounds(tmp_bounds);
        qpoases_solver_->init(
            qpoases_hessian_.get(), g_->get_values(), qpoases_jacobian_.get(),
            lb_->get_values(), ub_->get_values(), lbA_->get_values(),
            ubA_->get_values(), num_qp_iterations, 0, x_qp_->get_values(),
            y_qp_->get_values(), &tmp_bounds);
        new_QP_matrix_status_ = old_QP_matrix_status_ = UNDEFINED;
      }
    }
    retval = get_qpoases_exit_status_(qpoases_solver_);
#endif
  }

  // Get the solver status from qpOASES
  solver_status_ = get_qpoases_exit_status_();

  // reset the flag that remembers whether QP matrices have changed
  qp_matrices_changed_ = false;

  // Check if the solve was successful
  if (solver_status_ == QPEXIT_OPTIMAL) {
    first_qp_solved_ = true;
    solve_successful = true;

    // Retrieve the primal solution from qpOASES
    qpoases_solver_->getPrimalSolution(
        primal_solution_->get_non_const_values());

    // Retrieve the dual solution.  qpOASES returns bound and constraint
    // multipliers in a single array, so we need to take it apart.
    Vector qpOASES_dual_solution(num_qp_variables_ + num_qp_constraints_);
    qpoases_solver_->getDualSolution(
        qpOASES_dual_solution.get_non_const_values());

    // Get the values for the bound multipliers
    bound_multipliers_->copy_from_subvector(qpOASES_dual_solution, 0);

    // Get the values for the constraint multipliers
    constraint_multipliers_->copy_from_subvector(qpOASES_dual_solution,
                                                 num_qp_variables_);

  } else {
    string qpOASES_error_message;
    switch (solver_status_) {
      case QPEXIT_INFEASIBLE:
        qpOASES_error_message = "INFEASIBLE";
        break;
      case QPEXIT_UNBOUNDED:
        qpOASES_error_message = "UNBOUNDED";
        break;
      case QPEXIT_INTERNAL_ERROR:
        qpOASES_error_message = "INTERNAL ERROR";
        break;
      default:
        qpOASES_error_message = "NEED TO ADD TO OUTPUT";
    }
    jnlst_->Printf(J_ERROR, J_MAIN, "qpOASES error (%d): %s\n",
                   (int)solver_status_, qpOASES_error_message.c_str());
    write_qp_data_to_file("qpoases_failure_qp.txt");
    // assert(false && "Not trying to fix qpOASES solution yet");
    // handle_error(QP_TYPE_QP, stats);
    solve_successful = false;
  }

  // Update solver statistics
  if (stats != nullptr) {
    stats->increase_qp_iteration_counter((int)num_qp_iterations);
  }

  return solve_successful;
}

//@}
void QpOasesInterface::set_objective_hessian(
    shared_ptr<const SparseTripletMatrix> triplet_matrix)
{
  assert(qp_type_ == QP_TYPE_QP);

  // Remember that a change is made to the Hessian matrix
  qp_matrices_changed_ = true;

  if (!hessian_->is_initialized()) {
    // In this case, the Hessian matrix has never been set before and we need to
    // initialize the nonzero structure
    hessian_->set_structure(triplet_matrix);

    // We also need to copy the values
    hessian_->set_values(triplet_matrix);

    // Now we create and initialize the qpOASES version of the Hessian matrix.
    //
    // qpOASES requires those arrays as non-const (even though it should be
    // const, bad design)
    int* H_const_row_indices = const_cast<int*>(hessian_->get_row_indices());
    int* H_const_col_indices = const_cast<int*>(hessian_->get_column_indices());
    double* H_const_values = const_cast<double*>(hessian_->get_values());

    assert(!qpoases_hessian_);
    // IMPORTANT:  This creates a shallow copy of our matrix!
    qpoases_hessian_ = make_shared<qpOASES::SymSparseMat>(
        num_qp_variables_, num_qp_variables_, H_const_row_indices,
        H_const_col_indices, H_const_values);

    // Need to do some additional initialization for the qpOASES matrix
    qpoases_hessian_->createDiagInfo();
  } else {
    // Set the values in our matrix
    hessian_->set_values(triplet_matrix);

    // Since the qpOASES matrix is a shallow copy, we do not need to copy the
    // values again!
    assert(qpoases_hessian_);
  }
}

void QpOasesInterface::set_constraint_jacobian(
    shared_ptr<const RestartSqp::SparseTripletMatrix> triplet_matrix,
    IdentityMatrixPositions& identity_matrix_positions)
{
  if (!jacobian_) {
    // If there are no constraints (and no Jacobian) there is nothing to do
    return;
  }

  // Remember that a change is made to the Jacobian matrix
  qp_matrices_changed_ = true;

  if (!jacobian_->is_initialized()) {
    // In this case, the Jacobian matrix has never been set before and we need
    // to initialize the nonzero structure
    jacobian_->set_structure(triplet_matrix, identity_matrix_positions);

    // We also need to copy the values
    jacobian_->set_values(triplet_matrix);

    // Now we create and initialize the qpOASES version of the Jacobian matrix
    //
    // qpOASES requires those arrays as non-const (even though it should be
    // const, bad design)
    int* A_const_row_indices = const_cast<int*>(jacobian_->get_row_indices());
    int* A_const_col_indices =
        const_cast<int*>(jacobian_->get_column_indices());
    double* A_const_values = const_cast<double*>(jacobian_->get_values());

    qpoases_jacobian_ = make_shared<qpOASES::SparseMatrix>(
        num_qp_constraints_, num_qp_variables_, A_const_row_indices,
        A_const_col_indices, A_const_values);
  } else {
    // Set the values in our matrix
    jacobian_->set_values(triplet_matrix);

    // Since the qpOASES matrix is a shallow copy, we do not need to copy the
    // values again!
    assert(qpoases_jacobian_);
  }
}

#if 0
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
  if (jacobian_ != nullptr) {
    jacobian_->multiply(x_qp_, Ax); // tmp_vec_nCon=A*x
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
  if (jacobian_ != nullptr) {
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

  if (jacobian_ != nullptr) {
    jacobian_->multiply_transpose(y_qp_->get_values() + num_qp_variables_,
                                  stationary_gap->get_values());
  }
  shared_ptr<Vector> Hx = make_shared<Vector>(num_qp_variables_);
  hessian_->multiply(x_qp_, Hx);

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
  if (jacobian_ != nullptr) {
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
#endif

#if 0
void qpOASESInterface::handle_error(QPType qptype, shared_ptr<Statistics> stats)
{
  assert("Not trying to fix qpOASES solutions yet" && false);
  if (qptype == QP_TYPE_LP) {
    qpOASES::int_t num_qp_iterations =
        lp_solver_max_num_iterations_; // TODO modify it
    if (qpoases_solver_->isInfeasible()) {
      shared_ptr<Vector> x_0 = make_shared<Vector>(num_qp_variables_);
      shared_ptr<Vector> Ax = make_shared<Vector>(num_qp_constraints_);
      x_0->copy_vector(x_qp_);
      jacobian_->multiply(x_0, Ax);

      for (int i = 0; i < num_qp_constraints_; i++) {
        x_0->set_value(i + num_qp_variables_ - 2 * num_qp_constraints_,
                       max(0.0, lbA_->get_value(i)));
        x_0->set_value(i + num_qp_variables_ - num_qp_constraints_,
                       -min(0.0, ubA_->get_value(i)));
      }
      qpoases_solver_->init(
          qpoases_hessian_.get(), g_->get_values(), qpoases_jacobian_.get(),
          lb_->get_values(), ub_->get_values(), lbA_->get_values(),
          ubA_->get_values(), num_qp_iterations, NULL, x_0->get_values());
    } else {
      qpoases_solver_->init(0, g_->get_values(), qpoases_jacobian_.get(),
                            lb_->get_values(), ub_->get_values(),
                            lbA_->get_values(), ubA_->get_values(),
                            num_qp_iterations);
    }
    old_QP_matrix_status_ = new_QP_matrix_status_ = UNDEFINED;
    if (stats != nullptr)
      stats->increase_qp_iteration_counter((int)num_qp_iterations);

    if (!qpoases_solver_->isSolved()) {
      //  THROW_EXCEPTION(LP_NOT_OPTIMAL, LP_NOT_OPTIMAL_MSG);
    }
  } else {
    qpOASES::int_t num_qp_iterations =
        qp_solver_max_num_iterations_; // TODO modify it
    if (qpoases_solver_->isInfeasible()) {
      shared_ptr<Vector> x_0 = make_shared<Vector>(num_qp_variables_);
      for (int i = 0; i < num_qp_constraints_; i++) {
        x_0->set_value(i + num_qp_variables_ - 2 * num_qp_constraints_,
                       max(0.0, lbA_->get_value(i)));
        x_0->set_value(i + num_qp_variables_ - num_qp_constraints_,
                       -min(0.0, ubA_->get_value(i)));
      }
      qpoases_solver_->init(
          qpoases_hessian_.get(), g_->get_values(), qpoases_jacobian_.get(),
          lb_->get_values(), ub_->get_values(), lbA_->get_values(),
          ubA_->get_values(), num_qp_iterations, NULL, x_0->get_values());
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
      qpoases_solver_->init(qpoases_hessian_.get(), g_->get_values(),
                            qpoases_jacobian_.get(), lb_->get_values(),
                            ub_->get_values(), lbA_->get_values(),
                            ubA_->get_values(), num_qp_iterations);
    }
    old_QP_matrix_status_ = new_QP_matrix_status_ = UNDEFINED;
    if (stats != nullptr)
      stats->increase_qp_iteration_counter((int)num_qp_iterations);
    if (!qpoases_solver_->isSolved()) {
      // THROW_EXCEPTION(QP_NOT_OPTIMAL, QP_NOT_OPTIMAL_MSG);
    }
  }
}
#endif

void QpOasesInterface::set_qp_solver_options_()
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

void QpOasesInterface::retrieve_working_set_()
{
  assert(solver_status_ == QPEXIT_OPTIMAL);

  assert(num_qp_variables_ == qpoases_solver_->getNV());
  qpOASES::int_t* tmp_W_b = new int[num_qp_variables_];
  qpoases_solver_->getWorkingSetBounds(tmp_W_b);

  for (int i = 0; i < num_qp_variables_; i++) {
    switch ((int)tmp_W_b[i]) {
      case 1:
        if (lower_variable_bounds_->get_value(i) ==
            upper_variable_bounds_->get_value(i)) {
          bounds_working_set_[i] = ACTIVE_EQUALITY;
        } else {
          bounds_working_set_[i] = ACTIVE_ABOVE;
        }
        break;
      case -1:
        if (lower_variable_bounds_->get_value(i) ==
            upper_variable_bounds_->get_value(i)) {
          bounds_working_set_[i] = ACTIVE_EQUALITY;
        } else {
          bounds_working_set_[i] = ACTIVE_BELOW;
        }
        break;
      case 0:
        bounds_working_set_[i] = INACTIVE;
        break;
      default:
        assert("Invalid working set flag" && false);
    }
  }
  delete[] tmp_W_b;

  assert(num_qp_constraints_ == qpoases_solver_->getNC());
  qpOASES::int_t* tmp_W_c = new int[num_qp_constraints_];
  qpoases_solver_->getWorkingSetConstraints(tmp_W_c);

  for (int i = 0; i < num_qp_constraints_; i++) {
    switch ((int)tmp_W_c[i]) {
      case 1:
        if (lower_constraint_bounds_->get_value(i) ==
            upper_constraint_bounds_->get_value(i)) {
          constraints_working_set_[i] = ACTIVE_EQUALITY;
        } else {
          constraints_working_set_[i] = ACTIVE_ABOVE;
        }
        break;
      case -1:
        if (lower_constraint_bounds_->get_value(i) ==
            upper_constraint_bounds_->get_value(i)) {
          constraints_working_set_[i] = ACTIVE_EQUALITY;
        } else {
          constraints_working_set_[i] = ACTIVE_BELOW;
        }
        break;
      case 0:
        constraints_working_set_[i] = INACTIVE;
        break;
      default:
        assert("Invalid working set flag" && false);
    }
  }
  delete[] tmp_W_c;
}

} // SQPHOTSTART
