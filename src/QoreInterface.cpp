/* Copyright (C) 2019
 *
 * Authors: Xinyi Luo
 * Date:    2019-08-15
 */

#include "restartsqp/QoreInterface.hpp"

using namespace std;
using namespace Ipopt;

namespace RestartSqp {

QoreInterface::QoreInterface(int num_qp_variables, int num_qp_constraints,
                             QPType qp_type,
                             SmartPtr<const OptionsList> options,
                             SmartPtr<Journalist> jnlst)
 : QpSolverInterface(qp_type, num_qp_variables, num_qp_constraints, jnlst)
 , qore_solver_(nullptr)
 , first_qp_solved_(false)
 , qp_matrices_changed_(true)
 , qore_primal_solution_(NULL)
 , qore_dual_solution_(NULL)
{
  /** Get the options for this object. */
  get_option_values_(options);

  // We need to delay the creation of the QORE object since we need to know the
  // number of nonzeros in the matrices for that

  // Create matrix objects (without data)
  if (num_qp_constraints_ > 0) {
    bool is_compressed_row = true;
    bool is_symmetric = false;
    jacobian_ =
        make_shared<SparseHbMatrix>(num_qp_constraints_, num_qp_variables_,
                                    is_compressed_row, is_symmetric);
  }

  if (qp_type_ == QP_TYPE_QP) {
    bool is_compressed_row = true;
    bool is_symmetric = true;
    hessian_ = make_shared<SparseHbMatrix>(num_qp_variables_, num_qp_variables_,
                                           is_compressed_row, is_symmetric);
  }

  // Allocate memory for QORE's optimal solutions. */
  qore_primal_solution_ = new double[num_qp_variables + num_qp_constraints];
  qore_dual_solution_ = new double[num_qp_variables + num_qp_constraints];
}

void QoreInterface::create_qore_solver_(int num_nnz_jacobian,
                                        int num_nnz_hessian)
{
  assert(!qore_solver_);

  // Create the new qore_solver object
  int rv = QPNew(&qore_solver_, num_qp_variables_, num_qp_constraints_,
                 num_nnz_jacobian, num_nnz_hessian);
  if (rv != QPSOLVER_OK) {
    jnlst_->Printf(J_ERROR, J_MAIN, "QORE QPNew fails and returns rv = %d\n", rv);
    THROW_EXCEPTION(SQP_EXCEPTION_QP_SOLVER_FAILS, "QORE fails in initialization");
  }

  // Set the options in the QORE object
  set_qp_solver_options_();
}

void QoreInterface::get_option_values_(SmartPtr<const OptionsList> options)
{
  // Get the options from the options list
  options->GetIntegerValue("qp_solver_max_num_iterations",
                           qp_solver_max_num_iterations_, "");
  options->GetIntegerValue("qp_solver_print_level", qp_solver_print_level_, "");
  options->GetIntegerValue("lp_solver_max_num_iterations",
                           lp_solver_max_num_iterations_, "");
  options->GetBoolValue("qore_init_primal_variables", qore_init_primal_variables_,
                           "");
}

/**
 * @brief Destructor
 */
QoreInterface::~QoreInterface()
{
  if (qore_solver_) {
    QPFree(&qore_solver_);
  }

  delete[] qore_primal_solution_;
  delete[] qore_dual_solution_;
}

// Helper function to convert working set
static qp_int woringset_sqp_to_qore(ActivityStatus status_sqp)
{
  switch (status_sqp) {
    case ACTIVE_ABOVE:
      return -1;
      break;
    case ACTIVE_BELOW:
      return 1;
      break;
    case INACTIVE:
      return 0;
      break;
    default:
      assert("Invalid working set flag" && false);
  }
}

/**
 * @brief Solve a regular QP with given data and options.
 */

QpSolverExitStatus QoreInterface::optimize_impl(shared_ptr<Statistics> stats)
{
  // If not done earlier, allocate the QORE solver object
  if (!qore_solver_) {
    int num_nnz_jacobian = 0;
    if (jacobian_) {
      num_nnz_jacobian = jacobian_->get_num_entries();
    }
    int num_nnz_hessian = 0;
    if (hessian_) {
      num_nnz_hessian = hessian_->get_num_entries();
    }
    create_qore_solver_(num_nnz_jacobian, num_nnz_hessian);
  }

  // If the matrices have changed (or this is the first call), we need to set
  // the matrix data
  if (qp_matrices_changed_) {
    const int* A_row_indices = NULL;
    const int* A_col_indices = NULL;
    const double* A_values = NULL;
    if (jacobian_) {
      A_row_indices = jacobian_->get_row_indices();
      A_col_indices = jacobian_->get_column_indices();
      A_values = jacobian_->get_values();
    }
    const int* H_row_indices = NULL;
    const int* H_col_indices = NULL;
    const double* H_values = NULL;
    if (hessian_) {
      H_row_indices = hessian_->get_row_indices();
      H_col_indices = hessian_->get_column_indices();
      H_values = hessian_->get_values();
    }
    qp_int pflags = 0;
    int rv = QPSetData(qore_solver_, num_qp_variables_, num_qp_constraints_,
                       A_row_indices, A_col_indices, A_values, H_row_indices,
                       H_col_indices, H_values, pflags);
  }

  // Create arrays that include both variable and constraint bounds
  Vector qore_lb(num_qp_variables_ + num_qp_variables_);
  qore_lb.copy_into_subvector(lower_variable_bounds_, 0);
  qore_lb.copy_into_subvector(lower_constraint_bounds_, num_qp_variables_);

  Vector qore_ub(num_qp_variables_ + num_qp_variables_);
  qore_ub.copy_into_subvector(upper_variable_bounds_, 0);
  qore_ub.copy_into_subvector(upper_constraint_bounds_, num_qp_variables_);

  // Call the solver
  // TODO: Do we need to provide a starting point?
  // YES: for the primal variables, set them to zero.  Not for the dual
  // variables
  int qore_retval = 0;
  if (first_qp_solved_) {
    double* x_qore = NULL;
    if (qore_init_primal_variables_) {
    // We set the starting point for the primal variables (search direction and
    // penalty variables) to zero.
    // TODO QP: Should we set the penalty variables differently?
    x_qore = new double[num_qp_variables_ + num_qp_constraints_];
    for (int i = 0; i < num_qp_variables_ + num_qp_constraints_; ++i) {
      x_qore[i] = 0.;
    }
    }
    // We ask QORE to use the last dual solution as starting point
    const double* y_qore = NULL;
    // We ask QORE to use the working set from the last iteration
    const qp_int* ws_qore = NULL;
    qore_retval =
        QPOptimize(qore_solver_, qore_lb.get_values(), qore_ub.get_values(),
                   linear_objective_coefficients_->get_values(), x_qore, y_qore,
                   ws_qore, QPSOLVER_WARM);
    // Free memory
    delete [] x_qore;
  } else {
    // This is the first SQP iteration.
    //
    // We do not provide a starting point.
    const double* x_qore = NULL;
    const double* y_qore = NULL;

    // If the SQP user did not provide a starting point, we ask QORE to do a
    // cold start.
    // In that case
    // In COLD: Do not give starting point for primal and dual
    // In case SQP user provided initial working set, need to call WARM start
    // here
    //
    // TODO QP: Maybe we should not set the variables to zero if only the trust
    // region changed. Similar for second-order correction?  So, only if
    // matrices
    // changed?
    if (!bounds_working_set_) {
      assert(
          !constraints_working_set_ ||
          "Need to provide both variable and constraint working sets or none");
      const qp_int* ws_qore = NULL;
      qore_retval =
          QPOptimize(qore_solver_, qore_lb.get_values(), qore_ub.get_values(),
                     linear_objective_coefficients_->get_values(), x_qore,
                     y_qore, ws_qore, QPSOLVER_COLD);
    } else {
      assert(
          constraints_working_set_ ||
          "Need to provide both variable and constraint working sets or none");
      // If the user provided an initial working set, we translate that into
      // QORE's
      // structure and call it with a warm start.
      qp_int* ws_qore = new qp_int[num_qp_variables_ + num_qp_constraints_];
      // First set the flags for the bounds
      for (int i = 0; i < num_qp_variables_; ++i) {
        ws_qore[i] = woringset_sqp_to_qore(bounds_working_set_[i]);
      }
      // Now set the flags for the constraints
      for (int i = 0; i < num_qp_constraints_; ++i) {
        ws_qore[num_qp_variables_ + i] =
            woringset_sqp_to_qore(constraints_working_set_[i]);
      }
      // Now call QORE with WARM start
      qore_retval =
          QPOptimize(qore_solver_, qore_lb.get_values(), qore_ub.get_values(),
                     linear_objective_coefficients_->get_values(), x_qore,
                     y_qore, ws_qore, QPSOLVER_WARM);
      // Free memory
      delete [] ws_qore;
    }
  }
  // Get the solver status
  QpSolverExitStatus qp_solver_exit_status = get_qore_exit_status_();
  solver_status_ = qp_solver_exit_status;

  // reset the flag that remembers whether QP matrices have changed
  qp_matrices_changed_ = false;

  // Check if the solve was successful
  if (qp_solver_exit_status == QPEXIT_OPTIMAL) {
    first_qp_solved_ = true;

    // Retrieve the primal solution.  QORE's primal solution includes something
    // related to the constraints and is long.  We need to extract the first
    // part.
    int size = num_qp_variables_ + num_qp_constraints_;
    int rv = QPGetDbl(qore_solver_, QPSOLVER_DBL_PRIMALSOL,
                      qore_primal_solution_, &size);
    assert(rv == QPSOLVER_OK);

    // Get the values for the variables
    primal_solution_->copy_values(qore_primal_solution_);

    // Retrieve the dual solution.  They come in a single array and we need to
    // take them apart
    size = num_qp_variables_ + num_qp_constraints_;
    rv = QPGetDbl(qore_solver_, QPSOLVER_DBL_DUALSOL, qore_dual_solution_,
                  &size);
    assert(rv == QPSOLVER_OK);

    // Get the values for the bound multipliers
    bound_multipliers_->copy_values(qore_dual_solution_);

    // Get the values for the constraint multipliers
    constraint_multipliers_->copy_values(
        &qore_dual_solution_[num_qp_variables_]);

    // Determine the number of QP solver iterations
    int itercount;
    size = 1;
    rv = QPGetInt(qore_solver_, QPSOLVER_INT_ITERCOUNT, &itercount, &size);
    assert(rv == QPSOLVER_OK);
    qp_solver_iterations_ = itercount;
  } else {
    string qore_error_message;
    switch (qp_solver_exit_status) {
      case QPEXIT_INFEASIBLE:
        qore_error_message = "INFEASIBLE";
        break;
      case QPEXIT_UNBOUNDED:
        qore_error_message = "UNBOUNDED";
        break;
      case QPEXIT_INTERNAL_ERROR:
        qore_error_message = "INTERNAL ERROR";
        break;
      default:
        qore_error_message = "NEED TO ADD TO OUTPUT";
    }
    jnlst_->Printf(J_ERROR, J_MAIN, "qore error (%d): %s\n",
                   (int)qp_solver_exit_status, qore_error_message.c_str());
    write_qp_data_to_file("qore_failure_qp.txt");
  }

  // Update solver statistics
  if (stats != nullptr) {
    int num_qp_iterations;
    int size = 1;
    QPGetInt(qore_solver_, QPSOLVER_INT_ITERCOUNT, &num_qp_iterations, &size);
    stats->increase_qp_iteration_counter(num_qp_iterations);
  }

  return qp_solver_exit_status;
}

QpSolverExitStatus QoreInterface::get_qore_exit_status_()
{
  int qore_status, size = 1;
  QPGetInt(qore_solver_, QPSOLVER_INT_STATUS, &qore_status, &size);
  switch (qore_status) {
    case QPSOLVER_ITER_LIMIT:
      return QPEXIT_ITERLIMIT;
    case QPSOLVER_INFEASIBLE:
      return QPEXIT_INFEASIBLE;
    case QPSOLVER_UNBOUNDED:
      return QPEXIT_UNBOUNDED;
    case QPSOLVER_OPTIMAL:
      return QPEXIT_OPTIMAL;
    default:
      jnlst_->Printf(J_ERROR, J_MAIN, "QORE Error with return value = %d\n",
                     qore_status);
      return QPEXIT_INTERNAL_ERROR;
  }
}

void QoreInterface::retrieve_working_set_()
{
  assert(solver_status_ == QPEXIT_OPTIMAL);

  // Allocate memory to store the working set from QORE
  int size = num_qp_variables_ + num_qp_constraints_;
  int* working_set = new int[size];

  // Get the vector with the working set information from QORE
  QPGetInt(qore_solver_, QPSOLVER_INT_WORKINGSET, working_set, &size);

  // Allocate memory if that hadn't been done earlier
  if (!bounds_working_set_) {
    bounds_working_set_ = new ActivityStatus[num_qp_variables_];
  }

  for (int i = 0; i < num_qp_variables_; i++) {
    switch (working_set[i]) {
      case -1:
#if 0
        if (lower_variable_bounds_->get_value(i) ==
            upper_variable_bounds_->get_value(i)) {
          bounds_working_set_[i] = ACTIVE_EQUALITY;
        } else {
          bounds_working_set_[i] = ACTIVE_ABOVE;
        }
#endif
        bounds_working_set_[i] = ACTIVE_ABOVE;
        break;
      case 1:
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
      case 0:
        bounds_working_set_[i] = INACTIVE;
        break;
      default:
        assert("Invalid working set flag" && false);
    }
  }

  // Allocate memory if that hadn't been done earlier
  if (!constraints_working_set_) {
    constraints_working_set_ = new ActivityStatus[num_qp_constraints_];
  }

  for (int i = 0; i < num_qp_constraints_; i++) {
    switch (working_set[num_qp_variables_ + i]) {
      case -1:
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
      case 1:
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
      case 0:
        constraints_working_set_[i] = INACTIVE;
        break;
      default:
        assert("Invalid working set flag" && false);
    }
  }

  delete[] working_set;
}

void QoreInterface::set_qp_solver_options_()
{
  int value;
  if (qp_solver_print_level_ == 0) {
    // does not print anything
    value = -1;
    QPSetInt(qore_solver_, QPSOLVER_INT_PRTFREQ, &value, 1);
  }
  value = qp_solver_print_level_;
  QPSetInt(qore_solver_, QPSOLVER_INT_LOGLEVEL, &value, 1);
  QPSetInt(qore_solver_, QPSOLVER_INT_MAXITER, &qp_solver_max_num_iterations_,
           1);
}

void QoreInterface::set_constraint_jacobian(
    shared_ptr<const SparseTripletMatrix> triplet_matrix,
    IdentityMatrixPositions& identity_matrix_positions)
{
  // If there are no constraints (and no Jacobian), there is nothing to do
  if (!jacobian_) {
    return;
  }

  // Remember that a change is made to the Jacobian matrix
  qp_matrices_changed_ = true;

  if (!jacobian_->is_initialized()) {
    // In this case, the Jacobian matrix has never been set before and we need
    // to initialize the nonzero structure
    jacobian_->set_structure(triplet_matrix, identity_matrix_positions);
  }

  // Set the values in our matrix
  jacobian_->set_values(triplet_matrix);
}

void QoreInterface::set_objective_hessian(
    shared_ptr<const SparseTripletMatrix> triplet_matrix)
{
  assert(qp_type_ == QP_TYPE_QP);

  // Remember that a change is made to the Jacobian matrix
  qp_matrices_changed_ = true;

  if (!hessian_->is_initialized()) {
    // In this case, the Jacobian matrix has never been set before and we need
    // to initialize the nonzero structure
    hessian_->set_structure(triplet_matrix);
  }

  // Set the values in our matrix
  hessian_->set_values(triplet_matrix);
}

} // SQP_HOTSTART
