/* Copyright (C) 2019
 *
 * Authors: Xinyi Luo
 * Date:    2019-08-15
 */

#include "sqphot/QoreInterface.hpp"

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
  assert(rv == QPSOLVER_OK);

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
    int rv = QPSetData(qore_solver_, num_qp_variables_, num_qp_constraints_,
                       A_row_indices, A_col_indices, A_values, H_row_indices,
                       H_col_indices, H_values);
    assert(rv == QPSOLVER_OK);
  }

  // Pointers specifying the QORE starting point (NULL if none)
  double* qore_x;
  double* qore_y;

  if (first_qp_solved_) {
    // If this is not the first QP, we need to tell the solver whether the
    // matrices have changed.
    double strange_factor;
    if (qp_matrices_changed_) {
      strange_factor = 1.;
    } else {
      strange_factor = -1.;
    }
    int rv = QPAdjust(qore_solver_, strange_factor);
    assert(rv == QPSOLVER_OK);

    // Give QORE the starting point from the most recent solve
    qore_x = qore_primal_solution_;
    qore_y = qore_dual_solution_;
  } else {
    // In this case, we do not have a starting point yet
    qore_x = NULL;
    qore_y = NULL;
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
  int qore_retval =
      QPOptimize(qore_solver_, qore_lb.get_values(), qore_ub.get_values(),
                 linear_objective_coefficients_->get_values(), qore_x, qore_y);

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
    int rv = QPGetDblVector(qore_solver_, "primalsol", qore_primal_solution_);
    assert(rv == QPSOLVER_OK);

    // Get the values for the variables
    primal_solution_->copy_values(qore_primal_solution_);

    // Retrieve the dual solution.  They come in a single array and we need to
    // take them apart
    rv = QPGetDblVector(qore_solver_, "dualsol", qore_dual_solution_);
    assert(rv == QPSOLVER_OK);

    // Get the values for the bound multipliers
    bound_multipliers_->copy_values(qore_dual_solution_);

    // Get the values for the constraint multipliers
    constraint_multipliers_->copy_values(
        &qore_dual_solution_[num_qp_variables_]);

    // Determine the number of QP solver iterations
    int itercount;
    rv = QPGetInt(qore_solver_, "itercount", &itercount);
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
    QPGetInt(qore_solver_, "itercount", &num_qp_iterations);
    stats->increase_qp_iteration_counter(num_qp_iterations);
  }

  return qp_solver_exit_status;
}

QpSolverExitStatus QoreInterface::get_qore_exit_status_()
{
  int qore_status;
  QPGetInt(qore_solver_, "status", &qore_status);
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
  int* working_set = new int[num_qp_variables_ + num_qp_constraints_];

  // Get the vector with the working set information from QORE
  QPGetIntVector(qore_solver_, "workingset", working_set);

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

  if (qp_solver_print_level_ == 0) {
    // does not print anything
    QPSetInt(qore_solver_, "prtfreq", -1);
  }
  QPSetInt(qore_solver_, "maxiter", qp_solver_max_num_iterations_);
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
