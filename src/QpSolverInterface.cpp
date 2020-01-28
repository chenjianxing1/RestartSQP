/* Copyright (C) 2019
 * All Rights Reserved.
 *
 * Authors: Xinyi Luo
 * Date:    2019-07
 */
#include "sqphot/QpSolverInterface.hpp"

using namespace std;
using namespace Ipopt;

namespace RestartSqp {

QpSolverInterface::QpSolverInterface(QPType qp_type, int num_qp_variables,
                                     int num_qp_constraints,
                                     Ipopt::SmartPtr<Ipopt::Journalist> jnlst)
 : qp_type_(qp_type)
 , num_qp_variables_(num_qp_variables)
 , num_qp_constraints_(num_qp_constraints)
 , linear_objective_coefficients_(std::make_shared<Vector>(num_qp_variables))
 , lower_constraint_bounds_(std::make_shared<Vector>(num_qp_constraints))
 , upper_constraint_bounds_(std::make_shared<Vector>(num_qp_constraints))
 , lower_variable_bounds_(std::make_shared<Vector>(num_qp_variables))
 , upper_variable_bounds_(std::make_shared<Vector>(num_qp_variables))
 , solver_status_(QPEXIT_UNINITIALIZED)
 , primal_solution_(std::make_shared<Vector>(num_qp_variables))
 , constraint_multipliers_(std::make_shared<Vector>(num_qp_constraints))
 , bound_multipliers_(std::make_shared<Vector>(num_qp_variables))
 , jnlst_(jnlst)
 , file_counter_(1)
{
  working_set_up_to_date_ = false;
  bounds_working_set_ = new ActivityStatus[num_qp_variables];
  constraints_working_set_ = new ActivityStatus[num_qp_constraints];
}

QpSolverInterface::~QpSolverInterface()
{
  delete[] bounds_working_set_;
  bounds_working_set_ = NULL;
  delete[] constraints_working_set_;
  constraints_working_set_ = NULL;
}

QpSolverExitStatus QpSolverInterface::optimize(shared_ptr<Statistics> stats)
{
  bool write_all_qps = false;
  string qp_type_str;
  if (qp_type_ == QP_TYPE_QP) {
    qp_type_str = "qp";
  } else {
    qp_type_str = "lp";
  }
  if (write_all_qps) {
    write_qp_data_to_file("subprob_" + qp_type_str + to_string(file_counter_) +
                          ".txt");
    file_counter_++;
  }

  working_set_up_to_date_ = false;
  QpSolverExitStatus qp_solver_exit_status = optimize_impl(stats);
  solver_status_ = qp_solver_exit_status;
  if (qp_solver_exit_status == QPEXIT_OPTIMAL) {
    // TODO: Could postpone computation of working set until it is actually
    // needed, but this will require that some data members are declared
    // mutable
    retrieve_working_set_();
    working_set_up_to_date_ = true;
  }

  return qp_solver_exit_status;
}

KktError QpSolverInterface::calc_kkt_error(EJournalLevel level) const
{
  assert(solver_status_ == QPEXIT_OPTIMAL);
  assert(working_set_up_to_date_);

  bool is_nlp = false;
  KktError retval = calc_kkt_error_(
      jnlst_, level, is_nlp, lower_variable_bounds_, upper_variable_bounds_,
      lower_constraint_bounds_, upper_constraint_bounds_,
      linear_objective_coefficients_, nullptr, jacobian_, hessian_,
      primal_solution_, bound_multipliers_, constraint_multipliers_,
      bounds_working_set_, constraints_working_set_);

  return retval;
}

KktError calc_kkt_error_(SmartPtr<Journalist> jnlst, EJournalLevel level,
                         bool is_nlp,
                         shared_ptr<const Vector> lower_variable_bounds,
                         shared_ptr<const Vector> upper_variable_bounds,
                         shared_ptr<const Vector> lower_constraint_bounds,
                         shared_ptr<const Vector> upper_constraint_bounds,
                         shared_ptr<const Vector> linear_objective_coefficients,
                         shared_ptr<const Vector> constraint_values,
                         shared_ptr<const Matrix> jacobian,
                         shared_ptr<const Matrix> hessian,
                         shared_ptr<const Vector> primal_solution,
                         shared_ptr<const Vector> bound_multipliers,
                         shared_ptr<const Vector> constraint_multipliers,
                         const ActivityStatus* bounds_working_set,
                         const ActivityStatus* constraints_working_set)
{
  int num_variables = primal_solution->get_dim();
  int num_constraints = constraint_multipliers->get_dim();

  /*---------------------------------------------------------*/
  /*                    Primal Infeasibility                 */
  /*---------------------------------------------------------*/

  double primal_infeasibility = 0.;

  // First look at the bounds
  for (int i = 0; i < num_variables; i++) {
    primal_infeasibility =
        max(primal_infeasibility, max(0.0, lower_variable_bounds->get_value(i) -
                                               primal_solution->get_value(i)));
    primal_infeasibility =
        max(primal_infeasibility,
            max(0.0, primal_solution->get_value(i) -
                         upper_variable_bounds->get_value(i)));
  }

  shared_ptr<const Vector> constraint_body;

  // Compute the constraint body
  if (!is_nlp && jacobian) {
    // Compute the jacobian-vector project (this is for the QP)
    shared_ptr<Vector> Ax = make_shared<Vector>(num_constraints);
    Ax->set_to_zero();
    jacobian->multiply(primal_solution, Ax);
    constraint_body = Ax;
  } else {
    constraint_body = constraint_values;
  }

  for (int i = 0; i < num_constraints; i++) {
    primal_infeasibility = max(primal_infeasibility,
                               max(0.0, lower_constraint_bounds->get_value(i) -
                                            constraint_body->get_value(i)));
    primal_infeasibility =
        max(primal_infeasibility,
            max(0.0, constraint_body->get_value(i) -
                         upper_constraint_bounds->get_value(i)));
  }

  /*-------------------------------------------------------*/
  /*                    Dual Infeasibility                 */
  /*-------------------------------------------------------*/

  double dual_infeasibility = 0.;

  shared_ptr<Vector> lagrangian_gradient = make_shared<Vector>(num_variables);

  // Start with the linear part of the objective
  lagrangian_gradient->copy_vector(linear_objective_coefficients);

  // Add the Hessian part
  if (!is_nlp && hessian) {
    hessian->multiply(primal_solution, lagrangian_gradient);
  }

  // Subtract the bound mulipliers
  lagrangian_gradient->add_vector(-1., bound_multipliers);

  // Now the constraint multipliers
  if (jacobian) {
    jacobian->multiply_transpose(constraint_multipliers, lagrangian_gradient,
                                 -1.);
  }

  dual_infeasibility = lagrangian_gradient->calc_inf_norm();

  /*-------------------------------------------------------*/
  /*             Complementarity violation                 */
  /*-------------------------------------------------------*/

  double complementarity_violation = 0.;

  // Let's start with the variable bounds
  for (int i = 0; i < num_variables; ++i) {
    complementarity_violation =
        max(complementarity_violation,
            min(max(0., bound_multipliers->get_value(i)),
                primal_solution->get_value(i) -
                    lower_variable_bounds->get_value(i)));
    complementarity_violation =
        max(complementarity_violation,
            min(max(0., -bound_multipliers->get_value(i)),
                -primal_solution->get_value(i) +
                    upper_variable_bounds->get_value(i)));
  }

  // And now look at constraint bounds
  for (int i = 0; i < num_constraints; ++i) {
    complementarity_violation =
        max(complementarity_violation,
            min(max(0., constraint_multipliers->get_value(i)),
                constraint_body->get_value(i) -
                    lower_constraint_bounds->get_value(i)));
    complementarity_violation =
        max(complementarity_violation,
            min(max(0., -constraint_multipliers->get_value(i)),
                -constraint_body->get_value(i) +
                    upper_constraint_bounds->get_value(i)));
  }

  /*-------------------------------------------------------*/
  /*                       Working set                     */
  /*-------------------------------------------------------*/

  double working_set_error = 0.;

  if (bounds_working_set) {
    assert(constraints_working_set);

    // For bound constraints
    for (int i = 0; i < num_variables; ++i) {
      switch (bounds_working_set[i]) {
        case ACTIVE_ABOVE:
          working_set_error =
              max(working_set_error, fabs(primal_solution->get_value(i) -
                                          upper_variable_bounds->get_value(i)));
          break;
        case ACTIVE_BELOW:
          working_set_error =
              max(working_set_error, fabs(primal_solution->get_value(i) -
                                          lower_variable_bounds->get_value(i)));
          break;
        case ACTIVE_EQUALITY:
          working_set_error =
              max(working_set_error, fabs(primal_solution->get_value(i) -
                                          upper_variable_bounds->get_value(i)));
          working_set_error =
              max(working_set_error, fabs(primal_solution->get_value(i) -
                                          lower_variable_bounds->get_value(i)));
          break;
        case INACTIVE:
          // Nothing to check here
          break;
      }
    }

    // For regular constraints
    for (int i = 0; i < num_constraints; ++i) {
      switch (constraints_working_set[i]) {
        case ACTIVE_ABOVE:
          working_set_error = max(working_set_error,
                                  fabs(constraint_body->get_value(i) -
                                       upper_constraint_bounds->get_value(i)));
          break;
        case ACTIVE_BELOW:
          working_set_error = max(working_set_error,
                                  fabs(constraint_body->get_value(i) -
                                       lower_constraint_bounds->get_value(i)));
          break;
        case ACTIVE_EQUALITY:
          working_set_error = max(working_set_error,
                                  fabs(constraint_body->get_value(i) -
                                       upper_constraint_bounds->get_value(i)));
          working_set_error = max(working_set_error,
                                  fabs(constraint_body->get_value(i) -
                                       lower_constraint_bounds->get_value(i)));
          break;
        case INACTIVE:
          // Nothing to check here
          break;
      }
    }
  }

  // Compute the return value as the worst violation
  double worst_violation = primal_infeasibility;
  worst_violation = max(worst_violation, dual_infeasibility);
  worst_violation = max(worst_violation, complementarity_violation);
  worst_violation = max(worst_violation, working_set_error);

  // Write output
  string problem_name;
  if (is_nlp) {
    problem_name = "NLP";
  } else {
    problem_name = "QP";
  }
  if (jnlst->ProduceOutput(level, J_MAIN)) {
    jnlst->Printf(level, J_MAIN, "\nOptimality error in %s solution:\n",
                  problem_name.c_str());
    jnlst->Printf(level, J_MAIN,
                  "  Worst violation......................: %9.2e\n",
                  worst_violation);
    jnlst->Printf(level, J_MAIN,
                  "    Primal constraint violation........: %9.2e\n",
                  primal_infeasibility);
    jnlst->Printf(level, J_MAIN,
                  "    Dual infeasibility.................: %9.2e\n",
                  dual_infeasibility);
    jnlst->Printf(level, J_MAIN,
                  "    Complementarity vaiolation.........: %9.2e\n",
                  complementarity_violation);
    jnlst->Printf(level, J_MAIN,
                  "    Working set error..................: %9.2e\n\n",
                  working_set_error);
  }

  // Construct return value
  KktError retval;

  retval.primal_infeasibility = primal_infeasibility;
  retval.dual_infeasibility = dual_infeasibility;
  retval.complementarity_violation = complementarity_violation;
  retval.working_set_error = working_set_error;
  retval.worst_violation = worst_violation;

  return retval;
}

void QpSolverInterface::write_qp_data_to_file(const string& filename)
{
  FILE* file = fopen(filename.c_str(), "w");

  linear_objective_coefficients_->write_to_file(file, "g");
  lower_variable_bounds_->write_to_file(file, "lb");
  upper_variable_bounds_->write_to_file(file, "ub");
  lower_constraint_bounds_->write_to_file(file, "lbA");
  upper_constraint_bounds_->write_to_file(file, "ubA");
  if (hessian_) {
    hessian_->write_to_file(file, "H");
  }
  if (jacobian_) {
    jacobian_->write_to_file(file, "A");
  }

  fclose(file);
}

} // SQPHOTSTART
