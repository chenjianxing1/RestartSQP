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

#if 0
QoreInterface::QoreInterface(shared_ptr<SparseHbMatrix> H,
                             shared_ptr<SparseHbMatrix> A, shared_ptr<Vector> g,
                             shared_ptr<Vector> lb, shared_ptr<Vector> ub,
                             SmartPtr<const OptionsList> options)
 : A_(A)
 , H_(H)
 , g_(g)
 , lb_(lb)
 , ub_(ub)
 , num_qp_variables_(A->get_num_columns())
 , num_qp_constraints_(A->get_num_rows())
 , first_qp_solved_(false)
 , qore_solver_(0)
{
  get_option_values_(options);

  x_qp_ = make_shared<Vector>(num_qp_variables_ + num_qp_constraints_);
  y_qp_ = make_shared<Vector>(num_qp_variables_ + num_qp_constraints_);
  working_set_ = new int[num_qp_variables_ + num_qp_constraints_]();
  rv_ = QPNew(&qore_solver_, num_qp_variables_, num_qp_constraints_, A->get_num_entries(),
              H->get_num_entries());
  assert(rv_ == QPSOLVER_OK);
  set_qp_solver_options_();
}
#endif

/**
 * @brief Destructor
 */
QoreInterface::~QoreInterface()
{
  if (qore_solver_) {
    QPFree(&qore_solver_);
  }
}

/**
 * @brief Solve a regular QP with given data and options.
 */

bool QoreInterface::optimize_impl(shared_ptr<Statistics> stats)
{
  bool solve_successful = false;

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
      jacobian_->print("Jac");
      A_row_indices = jacobian_->get_row_indices();
      A_col_indices = jacobian_->get_column_indices();
      A_values = jacobian_->get_values();
    }
    const int* H_row_indices = NULL;
    const int* H_col_indices = NULL;
    const double* H_values = NULL;
    if (hessian_) {
      hessian_->print("Hess");
      H_row_indices = hessian_->get_row_indices();
      H_col_indices = hessian_->get_column_indices();
      H_values = hessian_->get_values();
    }
    int rv = QPSetData(qore_solver_, num_qp_variables_, num_qp_constraints_,
                       A_row_indices, A_col_indices, A_values, H_row_indices,
                       H_col_indices, H_values);
    assert(rv == QPSOLVER_OK);
  }

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
                 linear_objective_coefficients_->get_values(), 0, 0);

  // Get the solver status
  solver_status_ = get_qore_exit_status_(qore_retval);

  // reset the flag that remembers whether QP matrices have changed
  qp_matrices_changed_ = false;

  // Check if the solve was successful
  if (solver_status_ == QPEXIT_OPTIMAL) {
    first_qp_solved_ = true;
    solve_successful = true;

    // Retrieve the primal solution
    int rv = QPGetDblVector(qore_solver_, "primalsol",
                            primal_solution_->get_non_const_values());
    assert(rv == QPSOLVER_OK);

    // Retrieve the dual solution.  They come in a single array and we need to
    // take them apart
    Vector qore_dual_solution(num_qp_variables_ + num_qp_constraints_);
    rv = QPGetDblVector(qore_solver_, "dualsol",
                        qore_dual_solution.get_non_const_values());
    assert(rv == QPSOLVER_OK);

    // Get the values for the bound multipliers
    bound_multipliers_->copy_from_subvector(qore_dual_solution, 0);

    // Get the values for the constraint multipliers
    constraint_multipliers_->copy_from_subvector(qore_dual_solution,
                                                 num_qp_variables_);

  } else {
    string qore_error_message;
    switch (solver_status_) {
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
                   (int)solver_status_, qore_error_message.c_str());
    write_qp_data_to_file("qore_failure_qp.txt");
    // assert(false && "Not trying to fix qore solution yet");
    // handle_error(QP_TYPE_QP, stats);
    solve_successful = false;
  }

  // Update solver statistics
  if (stats != nullptr) {
    int num_qp_iterations;
    QPGetInt(qore_solver_, "itercount", &num_qp_iterations);
    stats->increase_qp_iteration_counter(num_qp_iterations);
  }

  return solve_successful;
}

QpSolverExitStatus QoreInterface::get_qore_exit_status_(int qore_retval)
{
  switch (qore_retval) {
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
                     qore_retval);
      return QPEXIT_INTERNAL_ERROR;
  }
}

void QoreInterface::retrieve_working_set_()
{
  assert(solver_status_ == QPEXIT_OPTIMAL);

  // Allocate memory to store the working set from QORE
  int* working_set = new int[num_qp_variables_ + num_qp_constraints_];

  QPGetIntVector(qore_solver_, "workingset", working_set);

  for (int i = 0; i < num_qp_variables_; i++) {
    switch (working_set[i]) {
      case -1:
        if (lower_variable_bounds_->get_value(i) ==
            upper_variable_bounds_->get_value(i)) {
          bounds_working_set_[i] = ACTIVE_EQUALITY;
        } else {
          bounds_working_set_[i] = ACTIVE_ABOVE;
        }
        break;
      case 1:
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

  for (int i = 0; i < num_qp_constraints_; i++) {
    switch (working_set[i]) {
      case -1:
        if (lower_constraint_bounds_->get_value(i) ==
            upper_constraint_bounds_->get_value(i)) {
          constraints_working_set_[i] = ACTIVE_EQUALITY;
        } else {
          constraints_working_set_[i] = ACTIVE_ABOVE;
        }
        break;
      case 1:
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
}

#if 0
void QoreInterface::handle_error(QPType qptype, shared_ptr<Statistics> stats)
{
  switch (status_) {
    case QPSOLVER_OPTIMAL:
      // do nothing here
      break;
    case QPSOLVER_INFEASIBLE:
      assert(false && "We don't handle QP solver yet");
      shared_ptr<Vector> x_0 = make_shared<Vector>(num_qp_variables_);
      // setup the slack variables to satisfy the bound constraints
      for (int i = 0; i < num_qp_constraints_; i++) {
        x_0->set_value(i + num_qp_variables_ - 2 * num_qp_constraints_,
                       max(0.0, lb_->get_value(num_qp_variables_ + i)));
        x_0->set_value(i + num_qp_variables_ - num_qp_constraints_,
                       -min(0.0, ub_->get_value(num_qp_variables_ + i)));
      }

      rv_ = QPOptimize(qore_solver_, lb_->get_values(), ub_->get_values(),
                       g_->get_values(), x_0->get_values(), NULL);

      int num_qp_iterations;
      QPGetInt(qore_solver_, "itercount", &num_qp_iterations);
      if (stats != nullptr)
        stats->increase_qp_iteration_counter(num_qp_iterations);

      assert(rv_ == QPSOLVER_OK);
      break;
  }
}
#endif

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
