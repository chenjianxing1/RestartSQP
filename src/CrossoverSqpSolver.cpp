/* Copyright (C) 2019
* All Rights Reserved.
*
* Authors: Xinyi Luo
* Date:2019-06
*/

#include "restartsqp/CrossoverSqpSolver.hpp"
#include "restartsqp/SqpInitSolveTNlp.hpp"
#include "IpIpoptApplication.hpp"

using namespace std;
using namespace Ipopt;

namespace RestartSqp {

/**
 * Default Constructor
 */
CrossoverSqpSolver::CrossoverSqpSolver()
{
  // First create the Ipopt application.  We will use the options and journalist
  // created in there
  ipopt_app_ = IpoptApplicationFactory();

  // Get the options and journalist
  SmartPtr<RegisteredOptions> reg_options = ipopt_app_->RegOptions();
  SmartPtr<OptionsList> options = ipopt_app_->Options();
  jnlst_ = ipopt_app_->Jnlst();

  // Register any extra options
  register_options_(reg_options);

  // Create an SqpSolver object that will be used for all active-set calls and
  // that is responsible for the journalist and options.
  sqp_solver_ = make_shared<SqpSolver>(reg_options, options, jnlst_);
}

/**
 * Destructor
 */
CrossoverSqpSolver::~CrossoverSqpSolver()
{
}

void CrossoverSqpSolver::determine_activities_(
    ActivityStatus* bound_activity_status,
    ActivityStatus* constraint_activity_status, IpoptSqpNlp& ipopt_tnlp)
{
  // Get the bounds
  double* variable_lower_bounds = new double[num_variables_];
  double* variable_upper_bounds = new double[num_variables_];
  double* constraint_lower_bounds = new double[num_constraints_];
  double* constraint_upper_bounds = new double[num_constraints_];
  bool retval = ipopt_tnlp.get_bounds_info(
      num_variables_, variable_lower_bounds, variable_upper_bounds,
      num_constraints_, constraint_lower_bounds, constraint_upper_bounds);
  assert(retval);

  // Get Ipopt solution
  const double* x_sol = ipopt_tnlp.x_sol();
  const double* z_L_sol = ipopt_tnlp.z_L_sol();
  const double* z_U_sol = ipopt_tnlp.z_U_sol();
  const double* g_sol = ipopt_tnlp.g_sol();
  const double* lambda_sol = ipopt_tnlp.lambda_sol();

  // Determine activity status for the variables based on slacks and multipliers
  const double active_bound_tol = 1e-4; // TODO, find more reliable tolerances
  const double active_mult_tol = 1e-8;

  jnlst_->Printf(J_DETAILED, J_MAIN,
                 "\nDetermine active variable bounds with tolderance %e:\n",
                 active_bound_tol);

  int num_variable_lower_bounds_active = 0;
  int num_variable_upper_bounds_active = 0;
  for (int i = 0; i < num_variables_; ++i) {
    // Check if this is a fixed variable
    if (variable_lower_bounds[i] == variable_upper_bounds[i]) {
      double bound_mult = z_L_sol[i] - z_U_sol[i];
      jnlst_->Printf(J_MOREDETAILED, J_MAIN,
                     " Variable %5d: FIXED: z_L = %e z_U = %e ", z_L_sol[i],
                     z_U_sol[i]);
      // Choose activity based on the sign of the multipliers
      if (bound_mult > active_mult_tol) {
        bound_activity_status[i] = ACTIVE_BELOW;
        jnlst_->Printf(J_MOREDETAILED, J_MAIN, " Active below\n");
        num_variable_lower_bounds_active++;
      } else if (bound_mult > active_mult_tol) {
        bound_activity_status[i] = ACTIVE_ABOVE;
        jnlst_->Printf(J_MOREDETAILED, J_MAIN, " Active above\n");
        num_variable_upper_bounds_active++;
      } else {
        bound_activity_status[i] = INACTIVE;
        jnlst_->Printf(J_MOREDETAILED, J_MAIN, " Inactive\n");
      }
      continue;
    }

    // Otherwise determine which side is active by looking at the ratio of the
    // slack and multipliers
    double slack_lower = x_sol[i] - variable_lower_bounds[i];
    double slack_upper = variable_upper_bounds[i] - x_sol[i];
    double mult_lower = max(1e-16, z_L_sol[i]);
    double mult_upper = max(1e-16, z_U_sol[i]);
    double ratio_lower = slack_lower / mult_lower;
    double ratio_upper = slack_upper / mult_upper;
    jnlst_->Printf(J_MOREDETAILED, J_MAIN,
                   " Variable %5d: LOWER: slack/mult = %e/%e=%e", i,
                   slack_lower, mult_lower, ratio_lower);
    jnlst_->Printf(J_MOREDETAILED, J_MAIN, "  UPPER: slack/mult = %e/%e=%e: ",
                   i, slack_upper, mult_upper, ratio_upper);

    if (ratio_lower > active_bound_tol && ratio_upper > active_bound_tol) {
      bound_activity_status[i] = INACTIVE;
      jnlst_->Printf(J_MOREDETAILED, J_MAIN, "Inactive\n");
    } else if (ratio_lower < ratio_upper) {
      bound_activity_status[i] = ACTIVE_BELOW;
      jnlst_->Printf(J_MOREDETAILED, J_MAIN, "Lower active\n");
      num_variable_lower_bounds_active++;
    } else {
      bound_activity_status[i] = ACTIVE_ABOVE;
      jnlst_->Printf(J_MOREDETAILED, J_MAIN, "Upper active\n");
      num_variable_upper_bounds_active++;
    }
  }
  jnlst_->Printf(J_MOREDETAILED, J_MAIN, "\n");
  jnlst_->Printf(J_SUMMARY, J_MAIN,
                 "Number of active lower variable bounds....: %d\n",
                 num_variable_lower_bounds_active);
  jnlst_->Printf(J_SUMMARY, J_MAIN,
                 "Number of active upper variable bounds....: %d\n",
                 num_variable_upper_bounds_active);

  // Determine activity status for the constraints based on slacks and
  // multipliers
  double constraint_active_tol = 1e-4; // TODO, find more reliable tolerances
  jnlst_->Printf(J_DETAILED, J_MAIN,
                 "\nDetermine active constraints with tolderance %e:\n",
                 constraint_active_tol);

  int num_constraint_lower_bound_active = 0;
  int num_constraint_upper_bound_active = 0;
  for (int i = 0; i < num_constraints_; ++i) {
    // Check if this is a fixed constraint
    if (constraint_lower_bounds[i] == constraint_upper_bounds[i]) {
      jnlst_->Printf(J_MOREDETAILED, J_MAIN,
                     " Constraint %5d: Equality: lambda = %e ", i,
                     lambda_sol[i]);
      // constraint_activity_status[i] = ACTIVE_EQUALITY;
      if (lambda_sol[i] > active_mult_tol) {
        jnlst_->Printf(J_MOREDETAILED, J_MAIN, "Assume lower side is active");
        constraint_activity_status[i] = ACTIVE_BELOW;
        num_constraint_lower_bound_active++;
      } else if (lambda_sol[i] < -active_mult_tol) {
        jnlst_->Printf(J_MOREDETAILED, J_MAIN, "Assume upper side is active");
        constraint_activity_status[i] = ACTIVE_ABOVE;
        num_constraint_upper_bound_active++;
      } else {
        jnlst_->Printf(J_MOREDETAILED, J_MAIN, "Assume inactive");
        constraint_activity_status[i] = INACTIVE;
      }
      continue;
    }

    // Otherwise determine which side is active by looking at the ratio of the
    // slack and multipliers
    double slack_lower = g_sol[i] - constraint_lower_bounds[i];
    double slack_upper = constraint_upper_bounds[i] - g_sol[i];
    double mult_lower = max(1e-16, lambda_sol[i]);
    double mult_upper = max(1e-16, -lambda_sol[i]);
    double ratio_lower = slack_lower / mult_lower;
    double ratio_upper = slack_upper / mult_upper;
    jnlst_->Printf(J_MOREDETAILED, J_MAIN,
                   " constraint %5d: LOWER: slack/mult = %e/%e=%e", i,
                   slack_lower, mult_lower, ratio_lower);
    jnlst_->Printf(J_MOREDETAILED, J_MAIN, "  UPPER: slack/mult = %e/%e=%e: ",
                   i, slack_upper, mult_upper, ratio_upper);

    if (ratio_lower > constraint_active_tol &&
        ratio_upper > constraint_active_tol) {
      constraint_activity_status[i] = INACTIVE;
      jnlst_->Printf(J_MOREDETAILED, J_MAIN, "Inactive\n");
    } else if (ratio_lower < ratio_upper) {
      constraint_activity_status[i] = ACTIVE_BELOW;
      jnlst_->Printf(J_MOREDETAILED, J_MAIN, "Lower active\n");
      num_constraint_lower_bound_active++;
    } else {
      constraint_activity_status[i] = ACTIVE_ABOVE;
      jnlst_->Printf(J_MOREDETAILED, J_MAIN, "Upper active\n");
      num_constraint_upper_bound_active++;
    }
  }
  delete[] variable_lower_bounds;
  delete[] variable_upper_bounds;
  delete[] constraint_lower_bounds;
  delete[] constraint_upper_bounds;

  jnlst_->Printf(J_MOREDETAILED, J_MAIN, "\n");
  jnlst_->Printf(J_SUMMARY, J_MAIN,
                 "Number of active lower constraints........: %d\n",
                 num_constraint_lower_bound_active);
  jnlst_->Printf(J_SUMMARY, J_MAIN,
                 "Number of active upper constraints........: %d\n\n",
                 num_constraint_upper_bound_active);
}

void CrossoverSqpSolver::initial_solve(shared_ptr<SqpTNlp> sqp_tnlp,
                                     const string& options_file_name)
{
  // Get information about the dimensions of the current NLP.
  num_variables_;
  num_constraints_;
  int num_nonzeros_jacobian;
  int num_nonzeros_hessian;
  std::string nlp_name;
  bool retval = sqp_tnlp->get_nlp_info(num_variables_, num_constraints_,
                                       num_nonzeros_jacobian,
                                       num_nonzeros_hessian, nlp_name);
  assert(retval);

  // Create an Ipopt NLP for solve the SqpTNlp first
  SmartPtr<IpoptSqpNlp> ipopt_tnlp = new IpoptSqpNlp(sqp_tnlp);

  // Initialize the Ipopt application (process options?)
  ApplicationReturnStatus status =
      ipopt_app_->Initialize(options_file_name.c_str());
  if (status != Solve_Succeeded) {
    assert(false && "We have not handled failure in Ipopt initialization.");
  }

  // Now let's solve the problem with Ipopt
  status = ipopt_app_->OptimizeTNLP(ipopt_tnlp);
  jnlst_->Printf(J_SUMMARY, J_MAIN,
                 "\nInitial solve with Ipopt returned status = %d\n\n", status);

  // If there is an Ipopt error, we stop
  switch (status) {
    case Solve_Succeeded:
    case Solved_To_Acceptable_Level:
      // This is all good
      break;
    default:
      // There is some situation that we haven't covered yet
      assert(false && "Still need to implement what to do when Ipopt fails.");
  }

  // Try to find out what the active sets are
  ActivityStatus* bounds_activity_status = new ActivityStatus[num_variables_];
  ActivityStatus* constraints_activity_status =
      new ActivityStatus[num_constraints_];
  determine_activities_(bounds_activity_status, constraints_activity_status,
                        *ipopt_tnlp);

  // Create an SqpTNlp in which we can overwrite the working set and the starting point.
  shared_ptr<SqpInitSolveTNlp> init_solve_tnlp = make_shared<SqpInitSolveTNlp>(sqp_tnlp);

  // Set the initial working set
  init_solve_tnlp->set_initial_working_sets(num_variables_, bounds_activity_status,
                                     num_constraints_,
                                     constraints_activity_status);
  // Free memory
  delete[] bounds_activity_status;
  delete[] constraints_activity_status;

  // We also need to set the Ipopt solution as the starting point
  const double* initial_primal_variables = ipopt_tnlp->x_sol();
  const double* z_L = ipopt_tnlp->z_L_sol();
  const double* z_U = ipopt_tnlp->z_U_sol();
  const double* initial_constraint_multipliers = ipopt_tnlp->lambda_sol();

  // Compute the max norm of the multipliers so that we can set the penalty
  // parameter large enough
  double max_mult = 0.; // Make based on some option
  for (int i = 0; i < num_constraints_; ++i) {
    max_mult = max(max_mult, abs(initial_constraint_multipliers[i]));
  }

  // We need to combine the Ipopt bound multipliers into one
  double* initial_bound_multipliers = new double[num_variables_];
  for (int i = 0; i < num_variables_; ++i) {
    initial_bound_multipliers[i] = z_L[i] - z_U[i];
   // max_mult = max(max_mult, abs(initial_bound_multipliers[i]));
  }

  // Give that starting point to the SQP TNLP so that the SQP solver can get it.
  init_solve_tnlp->set_starting_point(
      num_variables_, initial_primal_variables, initial_bound_multipliers,
      num_constraints_, initial_constraint_multipliers);

  delete[] initial_bound_multipliers;

  // Determine initial penalty parameter
  jnlst_->Printf(J_SUMMARY, J_MAIN,
                 "\nMax norm of optimal multipliers after Ipopt = %e\n",
                 max_mult);

  max_mult = max(max_mult, 10.);
  const double init_penalty_parameter_factor = 2.;
  double initial_penalty_parameter = init_penalty_parameter_factor * max_mult;
  jnlst_->Printf(J_SUMMARY, J_MAIN,
                 "  Setting initial penalty parameter to %e\n\n",
                 initial_penalty_parameter);
  SmartPtr<OptionsList> options = ipopt_app_->Options();

  options->SetNumericValue("penalty_parameter_init_value",
                           initial_penalty_parameter);

  // Now call the SQP solver to get the active-set solution for this problem
  const string local_options_file_name = "";  //CK: this is a duplicate of the function argument, which clang++ rejects
  bool keep_output_file = true;
  sqp_solver_->optimize_nlp(init_solve_tnlp, local_options_file_name, keep_output_file);

  // Check if the optimization was conclude successfully
  exit_flag_ = sqp_solver_->get_exit_flag();
  if (exit_flag_ != OPTIMAL) {
    // make sure all allocated memory is deleted
    return;
  }

  // Compute the max-norm of the multipliers to find a good value for the
  // penalty parameter to be used later
  shared_ptr<const Vector> current_constraint_multipliers =
      sqp_solver_->get_current_constraint_multipliers();
  // shared_ptr<const Vector> current_bound_multipliers =
  //    sqp_solver_->get_current_bound_multipliers();
  //
  // max_mult = max(current_constraint_multipliers->calc_inf_norm(),
  //                current_bound_multipliers->calc_inf_norm());
  max_mult = current_constraint_multipliers->calc_inf_norm();
  jnlst_->Printf(J_SUMMARY, J_MAIN,
                 "\nMax norm of optimal multipliers after SQP solve = %e\n",
                 max_mult);

  max_mult = max(max_mult, 10.);
  initial_penalty_parameter = init_penalty_parameter_factor * max_mult;
  jnlst_->Printf(J_SUMMARY, J_MAIN,
                 "  Setting initial penalty parameter to %e\n",
                 initial_penalty_parameter);

  options->SetNumericValue("penalty_parameter_init_value",
                           initial_penalty_parameter);
}

void CrossoverSqpSolver::next_solve(shared_ptr<SqpTNlp> sqp_tnlp)
{
  // Call the sqp_solver object to solve the new problem
  sqp_solver_->reoptimize_nlp(sqp_tnlp);
}

void CrossoverSqpSolver::register_options_(
    SmartPtr<RegisteredOptions> reg_options)
{
}

} // END_NAMESPACE_SQPHOTSTART
