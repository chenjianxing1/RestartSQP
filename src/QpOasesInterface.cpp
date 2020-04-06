/* Copyright (C) 2019
 * All Rights Reserved.
 *
 * Authors: Xinyi Luo
 * Date:    2019-07
 */
#include "RestartSQP/QpOasesInterface.hpp"
#include "RestartSQP/MessageHandling.hpp"

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
  qpoases_solver_ = make_shared<qpOASES::SQProblemSchur>(
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
QpSolverExitStatus QpOasesInterface::optimize_impl(shared_ptr<Statistics> stats)
{
  qpOASES::int_t num_qp_iterations =
      qp_solver_max_num_iterations_; // TODO modify it

  // First handle the case in which this is the first call
  if (!first_qp_solved_) {

    // If an active set has been provide, copy it into the qpOASES objects
    if (bounds_working_set_) {
      assert(
          constraints_working_set_ ||
          "Need to provide both variable and constraint working sets or none");

      // Create qpOASES bounds object to communicate the active set
      qpOASES::Bounds guessedBounds(num_qp_variables_);
      qpOASES::Constraints guessedConstraints(num_qp_constraints_);

      for (int i = 0; i < num_qp_variables_; ++i) {
        switch (bounds_working_set_[i]) {
          case ACTIVE_ABOVE:
            guessedBounds.setupBound(i, qpOASES::ST_UPPER);
            break;
          case ACTIVE_BELOW:
            guessedBounds.setupBound(i, qpOASES::ST_LOWER);
            break;
          case INACTIVE:
            guessedBounds.setupBound(i, qpOASES::ST_INACTIVE);
            break;
          default:
            assert(false && "Invalue activity value");
        }
      }

      for (int i = 0; i < num_qp_constraints_; ++i) {
        switch (constraints_working_set_[i]) {
          case ACTIVE_ABOVE:
            guessedConstraints.setupConstraint(i, qpOASES::ST_UPPER);
            break;
          case ACTIVE_BELOW:
            guessedConstraints.setupConstraint(i, qpOASES::ST_LOWER);
            break;
          case INACTIVE:
            guessedConstraints.setupConstraint(i, qpOASES::ST_INACTIVE);
            break;
          default:
            assert(false && "Invalue activity value");
        }
      }

      // Call the initial version that takes a working set.  Since no initial
      // starting point is given, qpOASES will initially start from the zero
      // solution.  We might want to change this later.
      double* cputime = nullptr;
      const double* xOpt = nullptr;
      const double* yOpt = nullptr;
      qpoases_solver_->init(
          qpoases_hessian_.get(), linear_objective_coefficients_->get_values(),
          qpoases_jacobian_.get(), lower_variable_bounds_->get_values(),
          upper_variable_bounds_->get_values(),
          lower_constraint_bounds_->get_values(),
          upper_constraint_bounds_->get_values(), num_qp_iterations, cputime,
          xOpt, yOpt, &guessedBounds, &guessedConstraints);

    } else {
      // In this case we call the version that does not use an initial working
      // set.
      qpoases_solver_->init(
          qpoases_hessian_.get(), linear_objective_coefficients_->get_values(),
          qpoases_jacobian_.get(), lower_variable_bounds_->get_values(),
          upper_variable_bounds_->get_values(),
          lower_constraint_bounds_->get_values(),
          upper_constraint_bounds_->get_values(), num_qp_iterations);
    }
  } else {
    // Here we already solved a QP from where we can hotstart
    assert(!constraints_working_set_ ||
           "Need to provide both variable and constraint working sets or none");

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
  }

  // Get the solver status from qpOASES
  QpSolverExitStatus qp_solver_exit_status = get_qpoases_exit_status_();

  // reset the flag that remembers whether QP matrices have changed
  qp_matrices_changed_ = false;

  // Check if the solve was successful
  if (qp_solver_exit_status == QPEXIT_OPTIMAL) {
    first_qp_solved_ = true;

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

    // Determine the number of QP solver iterations. The counter returned from
    // qpOASES does only count the updates, but we want to know how many linear
    // system solves, so we add one.
    qp_solver_iterations_ = num_qp_iterations + 1;
  } else {
    string qpOASES_error_message;
    switch (qp_solver_exit_status) {
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

    qp_solver_iterations_ = 0;
  }

  // Update solver statistics
  if (stats != nullptr) {
    stats->increase_qp_iteration_counter(qp_solver_iterations_);
  }

  return qp_solver_exit_status;
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

  // Allocate memory if that hadn't been done earlier
  if (!bounds_working_set_) {
    assert(!constraints_working_set_);
    bounds_working_set_ = new ActivityStatus[num_qp_variables_];
    constraints_working_set_ = new ActivityStatus[num_qp_constraints_];
  }

  // Get the working set for the bounds from qpOASES
  assert(num_qp_variables_ == qpoases_solver_->getNV());
  qpOASES::int_t* tmp_W_b = new int[num_qp_variables_];
  qpoases_solver_->getWorkingSetBounds(tmp_W_b);

  for (int i = 0; i < num_qp_variables_; i++) {
    switch ((int)tmp_W_b[i]) {
      case qpOASES::ST_UPPER:
#if 0
      if (lower_variable_bounds_->get_value(i) ==
            upper_variable_bounds_->get_value(i)) {
          bounds_working_set_[i] = ACTIVE_EQUALITY;  // we probably should keep information about the direction
        } else {
          bounds_working_set_[i] = ACTIVE_ABOVE;
        }
#endif
        bounds_working_set_[i] = ACTIVE_ABOVE;
        break;
      case qpOASES::ST_LOWER:
#if 0
        if (lower_variable_bounds_->get_value(i) ==
            upper_variable_bounds_->get_value(i)) {
          bounds_working_set_[i] = ACTIVE_EQUALITY;
        } else {
          bounds_working_set_[i] = ACTIVE_BELOW;
        }
#endif
        bounds_working_set_[i] = ACTIVE_BELOW;
        break;
      case qpOASES::ST_INACTIVE:
        bounds_working_set_[i] = INACTIVE;
        break;
      default:
        assert("Invalid working set flag" && false);
    }
  }
  delete[] tmp_W_b;

  // Get the constraint working set from qpOASES
  assert(num_qp_constraints_ == qpoases_solver_->getNC());
  qpOASES::int_t* tmp_W_c = new int[num_qp_constraints_];
  qpoases_solver_->getWorkingSetConstraints(tmp_W_c);

  for (int i = 0; i < num_qp_constraints_; i++) {
    switch ((int)tmp_W_c[i]) {
      case qpOASES::ST_UPPER:
#if 0
        if (lower_constraint_bounds_->get_value(i) ==
            upper_constraint_bounds_->get_value(i)) {
          constraints_working_set_[i] = ACTIVE_EQUALITY;
        } else {
          constraints_working_set_[i] = ACTIVE_ABOVE;
        }
#endif
        constraints_working_set_[i] = ACTIVE_ABOVE;
        break;
      case qpOASES::ST_LOWER:
#if 0
        if (lower_constraint_bounds_->get_value(i) ==
            upper_constraint_bounds_->get_value(i)) {
          constraints_working_set_[i] = ACTIVE_EQUALITY;
        } else {
          constraints_working_set_[i] = ACTIVE_BELOW;
        }
#endif
        constraints_working_set_[i] = ACTIVE_BELOW;
        break;
      case qpOASES::ST_INACTIVE:
        constraints_working_set_[i] = INACTIVE;
        break;
      default:
        assert("Invalid working set flag" && false);
    }
  }
  delete[] tmp_W_c;
}

} // RestartSQPSTART
