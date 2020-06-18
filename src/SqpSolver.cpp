/* Copyright (C) 2019
* All Rights Reserved.
*
* Authors: Xinyi Luo
* Date:2019-06
*/

#include "restartsqp/SqpSolver.hpp"
#include "restartsqp/MessageHandling.hpp"

#include <fstream>

/* TODOs:
 *
 * - check return values of sqp_tnlp calls
 * - decide whether to take out slack formulation (not maintained right now)
 * - do not introduce slack variables in QP that are not needed
 */

using namespace std;
using namespace Ipopt;

namespace RestartSqp {

DECLARE_STD_EXCEPTION(NEW_POINTS_WITH_INCREASE_OBJ_ACCEPTED);
DECLARE_STD_EXCEPTION(SMALL_TRUST_REGION);

enum WatchdogStatus
{
  WATCHDOG_INACTIVE = 0,
  WATCHDOG_READY = 1,
  WATCHDOG_IN_TRIAL_ITERATE = 2,
  WATCHDOG_SLEEPING = 3
};

enum StartingMode
{
  SM_PRIMAL_ONLY = 0,
  SM_PRIMAL_DUAL = 1,
  SM_WARM_START = 2
};

SqpSolver::SqpSolver(SmartPtr<RegisteredOptions> reg_options,
                     SmartPtr<OptionsList> options, SmartPtr<Journalist> jnlst)
 : current_output_file_name_("")
 , bound_type_(nullptr)
 , constraint_type_(nullptr)
 , reg_options_(reg_options)
 , options_(options)
 , jnlst_(jnlst)
 , force_warm_start_(false)
 , init_bound_activities_(nullptr)
 , init_constraint_activities_(nullptr)
{
  // Set the algorithm specific options
  bool skip_ipopt_options = true;
  register_options_(reg_options_, skip_ipopt_options);
}

SqpSolver::SqpSolver()
 : current_output_file_name_("")
 , bound_type_(nullptr)
 , constraint_type_(nullptr)
 , force_warm_start_(false)
 , init_bound_activities_(nullptr)
 , init_constraint_activities_(nullptr)
{
  // Some of the following code was taking from the Ipopt code.

  ////////////////////////////////////////////////////////////////////////////
  //              Set up output and options handling.                       //
  ////////////////////////////////////////////////////////////////////////////
  // Create journalist and set it up so that it prints to stdout.
  jnlst_ = new Journalist();
  SmartPtr<Journal> stdout_jrnl =
      jnlst_->AddFileJournal("console", "stdout", J_ITERSUMMARY);
  stdout_jrnl->SetPrintLevel(J_DBG, J_NONE);

  // Create a new options list
  options_ = new OptionsList();

  // Get the list of registered options to define the options for the
  // SQP algorithm .*/
  reg_options_ = new RegisteredOptions();

  // Finalize options list by giving it the list of registered options
  // and the journalist (for error message).
  options_->SetJournalist(jnlst_);
  options_->SetRegisteredOptions(reg_options_);

  // Set the algorithm specific options
  bool skip_ipopt_options = false;
  register_options_(reg_options_, skip_ipopt_options);
}

/**
 * Destructor
 */
SqpSolver::~SqpSolver()
{
  free_memory_();
}

/**
 *  Free any allocated memory.
 */
void SqpSolver::free_memory_()
{
  delete[] constraint_type_;
  constraint_type_ = nullptr;
  delete[] bound_type_;
  bound_type_ = nullptr;
  delete[] init_bound_activities_;
  init_bound_activities_ = nullptr;
  delete[] init_constraint_activities_;
  init_constraint_activities_ = nullptr;
}

void SqpSolver::initialize_options_(const string& options_file_name,
                                    bool keep_output_file)
{
  // Read options from file
  if (options_file_name != "") {
    try {
      ifstream is;
      is.open(options_file_name.c_str());
      bool allow_clobber = true; // The file can overwrite
      options_->ReadFromStream(*jnlst_, is, allow_clobber);
    } catch (...) {
      jnlst_->Printf(J_ERROR, J_MAIN, "Error readling options file.");
      assert(false && "Need to add proper exit for this");
    }
  }

  // If not file has been opened yet, reset keep_output_file name to
  // false since there is no output file to keep
  if (current_output_file_name_ == "") {
    keep_output_file = false;
  }

  if (!keep_output_file) {
    // Get output file name and file print level
    string new_output_file_name = "";
    int int_val;
    options_->GetIntegerValue("file_print_level", int_val, "");
    EJournalLevel file_print_level = (EJournalLevel)int_val;
    if (file_print_level != J_NONE) {
      options_->GetStringValue("output_file", new_output_file_name, "");
    }

    // If we already opened this file, there is nothing to do
    if (new_output_file_name != current_output_file_name_) {
      // There is currently no way to close a single Journal, so we close
      // everything now, including stdout.
      jnlst_->DeleteAllJournals();
      if (new_output_file_name != "") {
        // Open a new file journal
        jnlst_->AddFileJournal("output_file", new_output_file_name,
                               file_print_level);
      }
      current_output_file_name_ = new_output_file_name;
    }

    // Set the print level for stdout
    options_->GetIntegerValue("print_level", int_val, "");
    EJournalLevel print_level = (EJournalLevel)int_val;
    SmartPtr<Journal> stdout_jrnl = jnlst_->GetJournal("console");
    if (IsValid(stdout_jrnl)) {
      // Set print_level for stdout
      stdout_jrnl->SetAllPrintLevels(print_level);
    } else {
      jnlst_->AddFileJournal("console", "stdout", print_level);
    }
  }

  // Copy all options values into this object
  get_option_values_();
}

void SqpSolver::initialize_for_new_nlp_()
{
  // Free any existing objects or memory
  free_memory_();

  // Allocate all objects and memory needed in the loop.  This also creates the
  // QP and LP solver objects.
  allocate_memory_();

  // Determine the bound values
  sqp_nlp_->get_bounds_info(lower_variable_bounds_, upper_variable_bounds_,
                            lower_constraint_bounds_, upper_constraint_bounds_);

  // Determine inequality types
  classify_constraints_types_();
}

void SqpSolver::print_initial_output_()
{
  // Print problem statistics
  if (jnlst_->ProduceOutput(J_SUMMARY, J_MAIN)) {
    jnlst_->Printf(J_SUMMARY, J_MAIN, "Solving NLP %s\n\n",
                   sqp_nlp_->get_nlp_name().c_str());
    jnlst_->Printf(J_SUMMARY, J_MAIN,
                   "Number of variables.....................: %10d \n",
                   num_variables_);
    // Count the number of variables with bounds
    int num_vars_lower = 0;
    int num_vars_upper = 0;
    int num_vars_both = 0;
    int num_vars_fixed = 0;
    for (int i = 0; i < num_variables_; ++i) {
      switch (bound_type_[i]) {
        case BOUNDED_BELOW:
          num_vars_lower++;
          break;
        case BOUNDED_ABOVE:
          num_vars_upper++;
          break;
        case BOUNDED_BELOW_AND_ABOVE:
          num_vars_both++;
          break;
        case IS_EQUALITY:
          num_vars_fixed++;
          break;
        default:
          break;
      }
    }
    jnlst_->Printf(J_SUMMARY, J_MAIN,
                   "   with only lower bounds...............: %10d \n",
                   num_vars_lower);
    jnlst_->Printf(J_SUMMARY, J_MAIN,
                   "   with only upper bounds...............: %10d \n",
                   num_vars_upper);
    jnlst_->Printf(J_SUMMARY, J_MAIN,
                   "   with lower and upper bounds..........: %10d \n",
                   num_vars_both);
    jnlst_->Printf(J_SUMMARY, J_MAIN,
                   "   fixed................................: %10d \n",
                   num_vars_fixed);

    // Count the number of constraints of different types
    int num_cons_lower = 0;
    int num_cons_upper = 0;
    int num_cons_both = 0;
    int num_cons_equal = 0;
    for (int i = 0; i < num_constraints_; ++i) {
      switch (constraint_type_[i]) {
        case BOUNDED_BELOW:
          num_cons_lower++;
          break;
        case BOUNDED_ABOVE:
          num_cons_upper++;
          break;
        case BOUNDED_BELOW_AND_ABOVE:
          num_cons_both++;
          break;
        case IS_EQUALITY:
          num_cons_equal++;
          break;
        case UNBOUNDED:
          assert(false && "Invlid NLP: Constraint without bounds");
          break;
      }
    }
    jnlst_->Printf(J_SUMMARY, J_MAIN,
                   "Number of equality constraints..........: %10d \n",
                   num_cons_equal);
    jnlst_->Printf(J_SUMMARY, J_MAIN,
                   "Number of inequality constraints........: %10d \n",
                   num_constraints_ - num_cons_equal);
    jnlst_->Printf(J_SUMMARY, J_MAIN,
                   "   with only lower bounds...............: %10d \n",
                   num_cons_lower);
    jnlst_->Printf(J_SUMMARY, J_MAIN,
                   "   with only upper bounds...............: %10d \n",
                   num_cons_upper);
    jnlst_->Printf(J_SUMMARY, J_MAIN,
                   "   with lower and upper bounds..........: %10d \n\n",
                   num_cons_both);
  }

  static const char* s_qp_solver_name[2] = {"qpOASES", "QORE"};
  jnlst_->Printf(J_SUMMARY, J_MAIN,
                 "QP solver ..............................: %10s\n",
                 s_qp_solver_name[qp_solver_choice_]);
  jnlst_->Printf(J_SUMMARY, J_MAIN,
                 "Objective scaling factor................: %10.4e\n\n",
                 objective_scaling_factor_);

  // Print header of summary output table and the first line for the initial
  // iterate
  jnlst_->Printf(J_ITERSUMMARY, J_MAIN, "======================================"
                                        "======================================"
                                        "================================\n");
  jnlst_->Printf(J_ITERSUMMARY, J_MAIN,
                 "%6s %23s %9s %9s %8s %9s   %9s %9s  %5s %9s\n", "iter",
                 "objective", "||c_k||", "||p_k||", "Delta", "ratio", "pen par",
                 "QP_KKT", "QP it", "NLP_KKT");
  jnlst_->Printf(J_ITERSUMMARY, J_MAIN, "=================================="
                                        "=================================="
                                        "============="
                                        "===========================\n");

  double printable_obj_value =
      current_objective_value_ / objective_scaling_factor_;
  jnlst_->Printf(
      J_ITERSUMMARY, J_MAIN,
      "%6i %23.16e %10.3e %9.3e %9.3e %10.3e %9.3e %9.3e %5d %9.3e\n",
      solver_statistics_->num_sqp_iterations_, printable_obj_value,
      current_infeasibility_, 0., trust_region_radius_, 0.,
      current_penalty_parameter_, 0., 0, current_kkt_error_.worst_violation);
  // qp_solver_->get_QpOptimalStatus().KKT_error);

  double current_penalty_function_value =
      current_objective_value_ +
      current_penalty_parameter_ * current_infeasibility_;
  jnlst_->Printf(J_DETAILED, J_MAIN, "current_penalty_function_value: %23.16e\n", current_penalty_function_value);

  // Print some information about constant problem data
  if (jnlst_->ProduceOutput(J_VECTOR, J_MAIN)) {
    jnlst_->Printf(J_VECTOR, J_MAIN, "\n");
    jnlst_->Printf(J_VECTOR, J_MAIN, "Lower variable bounds:\n");
    lower_variable_bounds_->print("x_L", jnlst_, J_VECTOR, J_MAIN);
    jnlst_->Printf(J_VECTOR, J_MAIN, "\n");

    jnlst_->Printf(J_VECTOR, J_MAIN, "\n");
    jnlst_->Printf(J_VECTOR, J_MAIN, "Upper variable bounds:\n");
    upper_variable_bounds_->print("x_U", jnlst_, J_VECTOR, J_MAIN);
    jnlst_->Printf(J_VECTOR, J_MAIN, "\n");

    jnlst_->Printf(J_VECTOR, J_MAIN, "\n");
    jnlst_->Printf(J_VECTOR, J_MAIN, "Lower constraint bounds:\n");
    lower_constraint_bounds_->print("c_L", jnlst_, J_VECTOR, J_MAIN);
    jnlst_->Printf(J_VECTOR, J_MAIN, "\n");

    jnlst_->Printf(J_VECTOR, J_MAIN, "\n");
    jnlst_->Printf(J_VECTOR, J_MAIN, "Upper constraint bounds:\n");
    upper_constraint_bounds_->print("c_U", jnlst_, J_VECTOR, J_MAIN);
    jnlst_->Printf(J_VECTOR, J_MAIN, "\n\n");
  }

  // If desired, we print more information
  if (jnlst_->ProduceOutput(J_DETAILED, J_MAIN)) {
    jnlst_->Printf(J_DETAILED, J_MAIN, "\n");
    jnlst_->Printf(J_DETAILED, J_MAIN, "Values at starting point:\n\n");
    jnlst_->Printf(J_DETAILED, J_MAIN, "current_infeasibility.........: %e\n",
                   current_infeasibility_);
    jnlst_->Printf(J_DETAILED, J_MAIN, "trust region radius...........: %e\n",
                   trust_region_radius_);
    jnlst_->Printf(J_DETAILED, J_MAIN, "current objective value.......: %e\n",
                   current_objective_value_);
    jnlst_->Printf(J_DETAILED, J_MAIN, "penalty parameter.............: %e\n",
                   current_penalty_parameter_);
    jnlst_->Printf(J_DETAILED, J_MAIN, "\n");
  }

  // And if desired, we even print vector information
  if (jnlst_->ProduceOutput(J_VECTOR, J_MAIN)) {
    jnlst_->Printf(J_VECTOR, J_MAIN, "\n");
    jnlst_->Printf(J_VECTOR, J_MAIN, "Initial primal iterate:\n");
    current_iterate_->print("x_k", jnlst_, J_VECTOR, J_MAIN);
    jnlst_->Printf(J_VECTOR, J_MAIN, "\n");

    jnlst_->Printf(J_VECTOR, J_MAIN, "\n");
    jnlst_->Printf(J_VECTOR, J_MAIN,
                   "Initial values of the bound multipliers:\n");
    current_bound_multipliers_->print("bound multipliers", jnlst_, J_VECTOR,
                                      J_MAIN);
    jnlst_->Printf(J_VECTOR, J_MAIN, "\n");

    jnlst_->Printf(J_VECTOR, J_MAIN, "\n");
    jnlst_->Printf(J_VECTOR, J_MAIN,
                   "Initial values of the constraint multipliers:\n");
    current_constraint_multipliers_->print("constraint multipliers", jnlst_,
                                           J_VECTOR, J_MAIN);
    jnlst_->Printf(J_VECTOR, J_MAIN, "\n");

    jnlst_->Printf(J_VECTOR, J_MAIN, "\n");
    jnlst_->Printf(J_VECTOR, J_MAIN,
                   "Objective gradient at current iterate:\n");
    current_objective_gradient_->print("grad f(x_k)", jnlst_, J_VECTOR, J_MAIN);
    jnlst_->Printf(J_VECTOR, J_MAIN, "\n");

    jnlst_->Printf(J_VECTOR, J_MAIN, "\n");
    jnlst_->Printf(J_VECTOR, J_MAIN, "Constraint values at current iterate:\n");
    current_constraint_values_->print("c(x_k)", jnlst_, J_VECTOR, J_MAIN);
    jnlst_->Printf(J_VECTOR, J_MAIN, "\n");
  }

  // Finally, print matrix info if desired
  if (jnlst_->ProduceOutput(J_MATRIX, J_MAIN)) {
    jnlst_->Printf(J_MATRIX, J_MAIN, "\n");
    jnlst_->Printf(J_MATRIX, J_MAIN,
                   "Constraint Jacobian at current iterate:\n");
    current_constraint_jacobian_->print("Jac(x_k)", jnlst_, J_VECTOR, J_MAIN);
    jnlst_->Printf(J_MATRIX, J_MAIN, "\n");

    jnlst_->Printf(J_MATRIX, J_MAIN, "\n");
    jnlst_->Printf(J_MATRIX, J_MAIN,
                   "Lagrangian Hessian at current iterate:\n");
    current_lagrangian_hessian_->print("Hes(x_k)", jnlst_, J_VECTOR, J_MAIN);
    jnlst_->Printf(J_MATRIX, J_MAIN, "\n");
  }

  // Flush the buffer so that we see output even when the optimization is slow.
  jnlst_->FlushBuffer();

  // Clear the info string for the first iteration
  info_string_.clear();
}

void SqpSolver::calculate_search_direction_(
    shared_ptr<const Vector> qp_constraint_body)
{
  // Set the data for the QP corresponding to the current iterate
  setup_qp_(qp_constraint_body);

  // Solve the QP
  QpSolverExitStatus exit_status = qp_solver_->solve(solver_statistics_);
  if (exit_status != QPEXIT_OPTIMAL) {
    const string& nlp_name = sqp_nlp_->get_nlp_name();
    //qp_solver_->write_qp_data(nlp_name + "qpdata.log");

    assert(false && "Still need to decide how to handle QP solver error.");
    // exit_flag_ = qp_solver_->get_status();
    // break;
  }

  // Check if the QP found a bad local minimizer
  double qp_obj = qp_solver_->get_qp_objective();
  if (qp_obj > current_penalty_parameter_ * current_infeasibility_) {
    jnlst_->Printf(
        J_WARNING, J_MAIN, "WARNING: QP objective is %e, which is %e too large!\n",
        qp_obj, qp_obj - current_penalty_parameter_ * current_infeasibility_);
  }

  // Here, the trial_step_ vector has the length num_variables_, but the primal
  // solution of the QP solver includes the slack variables.
  trial_step_->copy_values(qp_solver_->get_primal_solution()->get_values());
  trial_step_->print("trial step", jnlst_, J_VECTOR, J_MAIN);

  // Compute the violation of the linearized constraints; this is used in
  // penalty update and computation of the predicted reduction.
  trial_model_infeasibility_ = calc_model_infeasibility_(trial_step_);

  // Get the multipliers from QP as trial multipliers
  trial_constraint_multipliers_->copy_vector(
      qp_solver_->get_constraint_multipliers());
  trial_bound_multipliers_->copy_vector(qp_solver_->get_bounds_multipliers());
}

void SqpSolver::print_iteration_output_()
{
  if (solver_statistics_->num_sqp_iterations_ % 10 == 0) {
    jnlst_->Printf(J_ITERSUMMARY, J_MAIN,
                   "%6s %23s %9s %9s %8s %9s   %9s %9s  %5s %9s\n", "iter",
                   "objective", "||c_k||", "||p_k||", "Delta", "ratio",
                   "pen par", "QP_KKT", "QP it", "NLP_KKT");
    jnlst_->Printf(J_ITERSUMMARY, J_MAIN, "=================================="
                                          "=================================="
                                          "============="
                                          "===========================\n");
  } else {
    jnlst_->Printf(J_DETAILED, J_MAIN,
                   "%6s %23s %9s %9s %8s %9s   %9s %9s  %5s %9s\n", "iter",
                   "objective", "||c_k||", "||p_k||", "Delta", "ratio",
                   "pen par", "QP_KKT", "QP it", "NLP_KKT");
    jnlst_->Printf(J_DETAILED, J_MAIN, "=================================="
                                       "=================================="
                                       "============="
                                       "===========================\n");
  }

  double qp_kkt_error = qp_solver_->get_qp_kkt_error();
  double model_ratio = actual_reduction_ / predicted_reduction_;
  int qp_iterations = qp_solver_->get_num_qp_iterations();
  double trial_step_norm = trial_step_->calc_inf_norm();
  double printable_obj_value =
      current_objective_value_ / objective_scaling_factor_;
  jnlst_->Printf(
      J_ITERSUMMARY, J_MAIN,
      "%6i %23.16e %10.3e %9.3e %9.3e %10.3e %9.3e %9.3e %5d %9.3e %s\n",
      solver_statistics_->num_sqp_iterations_, printable_obj_value,
      current_infeasibility_, trial_step_norm, trust_region_radius_,
      model_ratio, current_penalty_parameter_, qp_kkt_error, qp_iterations,
      current_kkt_error_.worst_violation, info_string_.c_str());

  double current_penalty_function_value =
      current_objective_value_ +
      current_penalty_parameter_ * current_infeasibility_;
  jnlst_->Printf(J_DETAILED, J_MAIN, "current_penalty_function_value: %23.16e\n", current_penalty_function_value);

  // If desired, we print more information
  if (jnlst_->ProduceOutput(J_DETAILED, J_MAIN)) {
    jnlst_->Printf(J_DETAILED, J_MAIN, "\n");
    jnlst_->Printf(J_DETAILED, J_MAIN,
                   "current objective value.......: %23.16e\n",
                   current_objective_value_);
    jnlst_->Printf(J_DETAILED, J_MAIN,
                   "trial objective value.........: %23.16e\n",
                   trial_objective_value_);
    jnlst_->Printf(J_DETAILED, J_MAIN,
                   "penalty parameter.............: %9.2e\n",
                   current_penalty_parameter_);
    jnlst_->Printf(J_DETAILED, J_MAIN,
                   "current_infeasibility.........: %9.2e\n",
                   current_infeasibility_);
    jnlst_->Printf(J_DETAILED, J_MAIN,
                   "trial_infeasibility...........: %9.2e\n",
                   trial_infeasibility_);
    jnlst_->Printf(J_DETAILED, J_MAIN,
                   "current model infeasibility...: %9.2e\n",
                   trial_model_infeasibility_);
    jnlst_->Printf(J_DETAILED, J_MAIN,
                   "trust region radius...........: %9.2e\n",
                   trust_region_radius_);
    jnlst_->Printf(J_DETAILED, J_MAIN,
                   "actual reduction..............: %9.2e\n",
                   actual_reduction_);
    jnlst_->Printf(J_DETAILED, J_MAIN,
                   "predicted reduction...........: %9.2e\n",
                   predicted_reduction_);
    jnlst_->Printf(J_DETAILED, J_MAIN, "\n");
  }

  // And if desired, we even print vector information
  if (jnlst_->ProduceOutput(J_VECTOR, J_MAIN)) {
    jnlst_->Printf(J_VECTOR, J_MAIN, "\n");
    jnlst_->Printf(J_VECTOR, J_MAIN, "Current iterate:\n");
    current_iterate_->print("x_k", jnlst_, J_VECTOR, J_MAIN);
    jnlst_->Printf(J_VECTOR, J_MAIN, "\n");

    jnlst_->Printf(J_VECTOR, J_MAIN, "\n");
    jnlst_->Printf(J_VECTOR, J_MAIN,
                   "Current values of the bound multipliers:\n");
    current_bound_multipliers_->print("bound multipliers", jnlst_, J_VECTOR,
                                      J_MAIN);
    jnlst_->Printf(J_VECTOR, J_MAIN, "\n");

    jnlst_->Printf(J_VECTOR, J_MAIN, "\n");
    jnlst_->Printf(J_VECTOR, J_MAIN,
                   "Current values of the constraint multipliers:\n");
    current_constraint_multipliers_->print("constraint multipliers", jnlst_,
                                           J_VECTOR, J_MAIN);
    jnlst_->Printf(J_VECTOR, J_MAIN, "\n");

    jnlst_->Printf(J_VECTOR, J_MAIN, "\n");
    jnlst_->Printf(J_VECTOR, J_MAIN,
                   "Objective gradient at current iterate:\n");
    current_objective_gradient_->print("grad f(x_k)", jnlst_, J_VECTOR, J_MAIN);
    jnlst_->Printf(J_VECTOR, J_MAIN, "\n");

    jnlst_->Printf(J_VECTOR, J_MAIN, "\n");
    jnlst_->Printf(J_VECTOR, J_MAIN, "Constraint values at current iterate:\n");
    current_constraint_values_->print("c(x_k)", jnlst_, J_VECTOR, J_MAIN);
    jnlst_->Printf(J_VECTOR, J_MAIN, "\n");
  }

  // Finally, print matrix info if desired
  if (jnlst_->ProduceOutput(J_MATRIX, J_MAIN)) {
    jnlst_->Printf(J_MATRIX, J_MAIN, "\n");
    jnlst_->Printf(J_MATRIX, J_MAIN,
                   "Constraint Jacobian at current iterate:\n");
    current_constraint_jacobian_->print("Jac(x_k)", jnlst_, J_VECTOR, J_MAIN);
    jnlst_->Printf(J_MATRIX, J_MAIN, "\n");

    jnlst_->Printf(J_MATRIX, J_MAIN, "\n");
    jnlst_->Printf(J_MATRIX, J_MAIN,
                   "Lagrangian Hessian at current iterate:\n");
    current_lagrangian_hessian_->print("Hes(x_k)", jnlst_, J_VECTOR, J_MAIN);
    jnlst_->Printf(J_MATRIX, J_MAIN, "\n");
  }

  // Flush the buffer so that we see output even when the optimization is slow.
  jnlst_->FlushBuffer();

  // Clear info string for next iteration.
  info_string_.clear();
}

void SqpSolver::optimize_nlp(shared_ptr<SqpTNlp> sqp_tnlp,
                             const string& options_file_name,
                             bool keep_output_file)
{
  // Initialize exit flag to UNKOWN to indicate that loop is not finished
  exit_flag_ = UNKNOWN_EXIT_STATUS;

  // Make NLP object directly accessible to all methods
  sqp_nlp_ = make_shared<SqpNlp>(sqp_tnlp);

  // Set the options based on the options file.
  initialize_options_(options_file_name, keep_output_file);

  // Set up everything needed to execute the loop: Read options, allocate
  // memory, create QP solver objects
  initialize_for_new_nlp_();

  // Initialize the iterates
  initialize_iterates_();

  // Now call the regular algorithm
  reoptimize_nlp(sqp_tnlp);
}

/**
 * @brief This is the main function to optimize the NLP given as the input
 *
 * @param nlp: the nlp reader that read data of the function to be minimized;
 */
void SqpSolver::reoptimize_nlp(shared_ptr<SqpTNlp> sqp_tnlp)
{
  // Initialize exit flag to UNKOWN to indicate that loop is not finished
  exit_flag_ = UNKNOWN_EXIT_STATUS;

  // Make NLP object directly accessible to all methods (does not anything new
  // if this is called from optimize_nlp()
  sqp_nlp_ = make_shared<SqpNlp>(sqp_tnlp);

  // Create objects that holds the solver statistics
  solver_statistics_ = make_shared<Statistics>();

  // Set the time at beginning of the algorithm
  cpu_time_at_start_ = get_cpu_time_since_start();
  wallclock_time_at_start_ = get_wallclock_time_since_start();

  // When we get here, the iterates have already been set, either in optimize_nlp
  // or from a previous solve.  All memory has been allocated, except for the
  // matrices (since the number of nonzeros might have changed) and the QP solver
  // objects.
  //
  // Initializes the QP and LP solvers
  initialize_qp_solvers_();

  // Initialize NLP function evaluations, bounds, and quantities computed from
  // NLP function values.
  compute_initial_values_();

  // Print initial output
  print_initial_output_();

  // Check if the starting point is already optimal.  This also computes the KKT
  // error for the current iterate

  // TODO: We need to solve at least one QP to get a working set to return
  // check_optimality_();

  // Main loop.  We catch exceptions to figure out final error code.
  try {
    while (solver_statistics_->num_sqp_iterations_ < max_num_iterations_ &&
           exit_flag_ == UNKNOWN_EXIT_STATUS) {

      // Wake up watchdog if it is time
      watchdog_sleep_iterations_++;
      if (watchdog_status_ == WATCHDOG_SLEEPING &&
          watchdog_sleep_iterations_ >= watchdog_min_wait_iterations_) {
        watchdog_status_ = WATCHDOG_READY;
      }

      // Solve the QP to get the trial step
      current_constraint_values_->print("cur c before main QP", jnlst_,
                                        J_VECTOR, J_MAIN);
      shared_ptr<const Vector> qp_constraint_body =
          current_constraint_values_;
      calculate_search_direction_(qp_constraint_body);

      // Update the penalty parameter if necessary.  In this process, a new
      // trial step might be computed.
      // We do not update the penalty parameter if we are in a watchdog trial step
      if (watchdog_status_ != WATCHDOG_IN_TRIAL_ITERATE) {
        update_penalty_parameter_();
      }

      // Calculate the trial point and evaluate the functions at that point.
      calc_trial_point_and_values_();

      // Check if the current iterate is acceptable
      perform_ratio_test_();

      // If the watchdog procedure is used, handle the different cases
      if (watchdog_status_ == WATCHDOG_READY ||
          watchdog_status_ == WATCHDOG_IN_TRIAL_ITERATE) {
        if (!trial_point_is_accepted_) {
          if (watchdog_status_ == WATCHDOG_READY) {
            // Activate the watchdog
            jnlst_->Printf(J_DETAILED, J_MAIN, "WATCHDOG: Activating watchdog.\n\n");
            watchdog_status_ = WATCHDOG_IN_TRIAL_ITERATE;
            info_string_ += "ws";
            // This stores the current iterates, search direction etc as backup
            store_watchdog_backups_();
            // We now temporarily accept the trial iterate
            trial_point_is_accepted_ = true;
          } else {
            assert(watchdog_status_ == WATCHDOG_IN_TRIAL_ITERATE);
            jnlst_->Printf(J_DETAILED, J_MAIN, "WATCHDOG: Watchdog iterate not accepted, restore original iterate.\n\n");
            info_string_ += "wr";
            // Restore values corresponding to iterate from which watchdog started.
            restore_watchdog_backups_();
            // Put watchdog to sleep
            watchdog_status_ = WATCHDOG_SLEEPING;
            // Reset the counter for consecutive itertions in which the watchdog has been sleeping.  We set it to -1 since it will be increased at the end of this loop by one
            watchdog_sleep_iterations_ = 0;
          }
        } else if (watchdog_status_ == WATCHDOG_IN_TRIAL_ITERATE) {
          // In this case the trial point was accepted and we reset the watchdog
          // flag
          jnlst_->Printf(J_DETAILED, J_MAIN, "WATCHDOG: Watchdog iterate accepted.\n\n");
          watchdog_status_ = WATCHDOG_READY;
          info_string_ += "wa";
          // and delete the backup
          delete_watchdog_backups_();
        }
      }

#if 0 // TODO : Need to integrate the second order correction with the watchdog
      // If requested by the options, do the second order correction step.
      if (!trial_point_is_accepted_ && perform_second_order_correction_step_) {
        second_order_correction_();
      }
#endif

      // Accept the new iterate
      if (trial_point_is_accepted_) {
        accept_trial_point_();
      }

      // Update the radius and the QP bounds if the radius has been changed
      solver_statistics_->increase_sqp_iteration_counter();

      // Compute the NLP KKT error and set exit_flag_ to indicate whether we
      // should stop..
      check_optimality_();

      // Print output
      print_iteration_output_();

      // exit the loop if required.
      if (exit_flag_ != UNKNOWN_EXIT_STATUS) {
        break;
      }

      // Update the trust region radius.  Even if we are in a trial watchdog
      // step, the trust region radius will be reduced since the predicted reduction
      // is negative.
      // TODO: Should we just leave TR radius in a trial step?  A short experiment
      //       seemed to indicate that this was inferior
      update_trust_region_radius_();
    }
  } catch (SQP_EXCEPTION_PENALTY_TOO_LARGE& exc) {
    exc.ReportException(*jnlst_, J_ERROR);
    exit_flag_ = PENALTY_TOO_LARGE;
  } catch (SQP_EXCEPTION_TRUST_REGION_TOO_SMALL& exc) {
    exc.ReportException(*jnlst_, J_ERROR);
    exit_flag_ = TRUST_REGION_TOO_SMALL;
  } catch (SQP_EXCEPTION_INFEASIBLE& exc) {
    exc.ReportException(*jnlst_, J_ERROR);
    exit_flag_ = QPERROR_INFEASIBLE;
  } catch (SQP_EXCEPTION_UNBOUNDED& exc) {
    exc.ReportException(*jnlst_, J_ERROR);
    exit_flag_ = QPERROR_UNBOUNDED;
  } catch (SQP_EXCEPTION_MAXITER& exc) {
    exc.ReportException(*jnlst_, J_ERROR);
    exit_flag_ = QPERROR_EXCEED_MAX_ITER;
  } catch (SQP_INVALID_WORKING_SET& exc) {
    exc.ReportException(*jnlst_, J_ERROR);
    exit_flag_ = INVALID_INITIAL_WORKING_SET;
  } catch (SQP_EXCEPTION_INTERNAL_ERROR& exc) {
    exc.ReportException(*jnlst_, J_ERROR);
    exit_flag_ = QPERROR_INTERNAL_ERROR;
  } catch (SQP_EXCEPTION_UNKNOWN& exc) {
    exc.ReportException(*jnlst_, J_ERROR);
    exit_flag_ = QPERROR_UNKNOWN;
  }

  // check if the current iterates get_status before exiting
  if (solver_statistics_->num_sqp_iterations_ == max_num_iterations_) {
    exit_flag_ = EXCEED_MAX_ITERATIONS;
  }

  // check if we are running out of CPU time
  double current_cpu_time = get_cpu_time_since_start();
  if (current_cpu_time - cpu_time_at_start_ > cpu_time_limit_) {
    exit_flag_ = EXCEED_MAX_CPU_TIME;
  }

  // check if we are running out of wallclock time
  double current_wallclock_time = get_wallclock_time_since_start();
  if (current_wallclock_time - wallclock_time_at_start_ >
      wallclock_time_limit_) {
    exit_flag_ = EXCEED_MAX_WALLCLOCK_TIME;
  }

  //    if (exitflag_ != OPTIMAL && exitflag_ != INVALID_NLP) {
  //        check_optimality();
  //    }

  // print the final summary message to the console
  print_final_stats_();
  jnlst_->FlushBuffer();

  // Return the information back to the user
  return_results_();
}

void SqpSolver::return_results_()
{
  // Get the activity status from the QP solver if one is available
  const ActivityStatus* bound_activity_status = NULL;
  const ActivityStatus* constraint_activity_status = NULL;
  if (qp_solver_->get_qp_solver_status() == QPEXIT_OPTIMAL) {
    bound_activity_status = qp_solver_->get_bounds_working_set();
    constraint_activity_status = qp_solver_->get_constraints_working_set();
  }

  // Undo any internal scaling for the final solution. */
  if (objective_scaling_factor_ != 1.) {
    current_bound_multipliers_->scale(1. / objective_scaling_factor_);
    current_constraint_multipliers_->scale(1. / objective_scaling_factor_);
    current_objective_value_ *= objective_scaling_factor_;
  }

  // Set the final value of the penalty parameter, unscaled
  solver_statistics_->set_final_penalty_parameter(current_penalty_parameter_/objective_scaling_factor_);

  sqp_nlp_->finalize_solution(
      exit_flag_, current_iterate_, current_bound_multipliers_,
      bound_activity_status, current_constraint_values_,
      current_constraint_multipliers_, constraint_activity_status,
      current_objective_value_, solver_statistics_);
}

/**
 * @brief This is the function that checks if the current point is optimal, and
 * decides if to exit the loop or not
 * *@return if it decides the function is optimal, the class member _exitflag =
 * OPTIMAL
 * if it decides that there is an error during the function run or the
 *  function cannot be solved, it will assign _exitflag the	corresponding
 *  code according to the error type.
 */
void SqpSolver::check_optimality_()
{
  // FIXME: not sure if it is better to use the new multiplier or the old one

  const ActivityStatus* bounds_working_set = NULL;
  const ActivityStatus* constraints_working_set = NULL;

  bool is_nlp = true;
  current_kkt_error_ =
      calc_kkt_error_(jnlst_, J_LAST_LEVEL, is_nlp, lower_variable_bounds_,
                      upper_variable_bounds_, lower_constraint_bounds_,
                      upper_constraint_bounds_, current_objective_gradient_,
                      current_constraint_values_, current_constraint_jacobian_,
                      nullptr, current_iterate_, current_bound_multipliers_,
                      current_constraint_multipliers_, bounds_working_set,
                      constraints_working_set);

  jnlst_->Printf(J_DETAILED, J_MAIN, "KKT Error in iteration %d:\n",
                 solver_statistics_->num_sqp_iterations_);
  jnlst_->Printf(J_DETAILED, J_MAIN, "    Primal infeasibility: %9.2e\n",
                 current_kkt_error_.primal_infeasibility);
  jnlst_->Printf(J_DETAILED, J_MAIN, "      Dual infeasibility: %9.2e\n",
                 current_kkt_error_.dual_infeasibility);
  jnlst_->Printf(J_DETAILED, J_MAIN, "   Complementarity error: %9.2e\n",
                 current_kkt_error_.complementarity_violation);
  jnlst_->Printf(J_DETAILED, J_MAIN, "           Overall error: %9.2e\n\n",
                 current_kkt_error_.worst_violation);

  bool is_optimal = true;
  if (current_kkt_error_.primal_infeasibility > opt_tol_primal_feasibility_) {
    is_optimal = false;
  } else if (current_kkt_error_.dual_infeasibility >
             opt_tol_dual_feasibility_) {
    is_optimal = false;
  } else if (current_kkt_error_.complementarity_violation >
             opt_tol_complementarity_) {
    is_optimal = false;
  }

  if (is_optimal) {
    exit_flag_ = OPTIMAL;
  }
}

void SqpSolver::calc_trial_point_and_values_()
{
  // Compute the trial iterate
  trial_iterate_->set_to_sum_of_vectors(1., current_iterate_, 1., trial_step_);

  // Evaluate the objective at the trial point
  bool retval = eval_f_(trial_iterate_, trial_objective_value_);
  if (!retval) {
    // We indicate that this is not an acceptable trial point by making the
    // objective function huge
    trial_objective_value_ = SqpInf;
  }

  // Evaluate the constraints at the trial point
  retval = eval_constraints_(trial_iterate_, trial_constraint_values_);
  assert(retval && "eval_constraints returned false");

  // Compute the constraint violation at the trial point
  trial_infeasibility_ =
      compute_constraint_violation_(trial_iterate_, trial_constraint_values_);
}

/**
 * @brief This function shifts the initial starting point to be feasible to the
 * bound constraints
 * @param x initial starting point
 * @param lower_variable_bounds_ lower bound constraints
 * @param upper_variable_bounds_ upper bound constraints
 */
// AW: This is only used here, so we can make it invisible to the outside
static void
shift_starting_point_(shared_ptr<Vector> x,
                      shared_ptr<const Vector> lower_variable_bounds_,
                      shared_ptr<const Vector> upper_variable_bounds_)
{
  for (int i = 0; i < x->get_dim(); i++) {
    assert(lower_variable_bounds_->get_value(i) <=
           upper_variable_bounds_->get_value(i));
    if (lower_variable_bounds_->get_value(i) > x->get_value(i)) {
      x->set_value(i, lower_variable_bounds_->get_value(i));
      //            x->set_value(i,
      //            lower_variable_bounds_->value(i)+0.5*(upper_variable_bounds_->value(i)-lower_variable_bounds_->value(i)));
    } else if (x->get_value(i) > upper_variable_bounds_->get_value(i)) {
      x->set_value(i, upper_variable_bounds_->get_value(i));
      //           x->set_value(i,
      //           upper_variable_bounds_->value(i)-0.5*(upper_variable_bounds_->value(i)-lower_variable_bounds_->value(i)));
    }
  }
}

void SqpSolver::initialize_iterates_()
{
  // Check if there are any degrees of freedom
  if (num_variables_ < num_equality_constraints_) {
    jnlst_->Printf(J_ERROR, J_MAIN, "ERROR: More equality constraints than variables.\n");
    THROW_EXCEPTION(SQP_NLP, "Too many equality constraints");
  }

  // Check if options was overwritten
  char mode = starting_mode_;
  if (force_warm_start_) {
    mode = SM_WARM_START;
    // Reset the overwrite flag
    force_warm_start_ = false;
  }

  // Get the primal dual starting point depending on the selected starting mode
  if (mode == SM_PRIMAL_ONLY) {
    // Get only the primal starting point
    sqp_nlp_->get_starting_point(current_iterate_, nullptr, nullptr);

    // Set the multipliers to zero
    current_bound_multipliers_->set_to_zero();
    current_constraint_multipliers_->set_to_zero();
  }
  else {
    // Get the primal-dual starting point
    sqp_nlp_->get_starting_point(current_iterate_, current_bound_multipliers_,
                                 current_constraint_multipliers_);
  }

  // Check if an initial working set is available.  If so, give it to the QP
  // solver
  bool use_initial_working_set;
  if (mode == SM_WARM_START) {
    use_initial_working_set = true;
  }
  else {
    use_initial_working_set = false;
  }

  if (use_initial_working_set) {
    // Make sure that the NLP can actually provide a starting point.
    if (!sqp_nlp_->use_initial_working_set()) {
      jnlst_->Printf(J_ERROR, J_MAIN, "ERROR: Warm start is chosen, but NLP does not provide initial working set.\n");
      THROW_EXCEPTION(SQP_INVALID_STARING_POINT, "Warm start is chosen, but NLP does not provide initial working set");
    }

    jnlst_->Printf(J_SUMMARY, J_MAIN,
                   "\nUser provided an initial working set:\n");
    // Store the working set in this object so that it can be given to the
    // QP solver when that one is initialized
    init_bound_activities_ = new ActivityStatus[num_variables_];
    init_constraint_activities_ =
        new ActivityStatus[num_constraints_];

    sqp_nlp_->get_initial_working_sets(num_variables_, init_bound_activities_,
                                       num_constraints_,
                                       init_constraint_activities_);

    // Let's count the number of active bounce to make sure there are not too many
    int num_active_lower_var = 0;
    int num_active_upper_var = 0;
    for (int i = 0; i < num_variables_; ++i) {
      if (init_bound_activities_[i] == ACTIVE_BELOW) {
        num_active_lower_var++;
      } else if (init_bound_activities_[i] == ACTIVE_ABOVE) {
        num_active_upper_var++;
      }
    }
    // Count the number of active inequality constraints
    int num_active_lower_con = 0;
    int num_active_upper_con = 0;
    int num_equality = 0;
    for (int i = 0; i < num_constraints_; ++i) {
      if (constraint_type_[i] == IS_EQUALITY) {
        if (init_constraint_activities_[i] != INACTIVE) {
          num_equality++;
        }
      } else if (init_constraint_activities_[i] == ACTIVE_BELOW) {
        num_active_lower_con++;
      } else if (init_constraint_activities_[i] == ACTIVE_ABOVE) {
        num_active_upper_con++;
      }
    }

    int num_total_active = num_active_lower_var + num_active_upper_var + num_active_lower_con + num_active_upper_con + num_equality;
    printf("num_total_active = %d num_var = %d\n", num_total_active, num_variables_);
    if (num_total_active > num_variables_) {
      THROW_EXCEPTION(SQP_INVALID_WORKING_SET,
                      "Too many constraints are active.");
    }

    // If desired, print the statistics of user-provided working set
    if (jnlst_->ProduceOutput(J_SUMMARY, J_MAIN)) {
      jnlst_->Printf(J_SUMMARY, J_MAIN,
                     "  Number of variables active at lower bound......: %d\n",
                     num_active_lower_var);
      jnlst_->Printf(J_SUMMARY, J_MAIN,
                     "  Number of variables active at upper bound......: %d\n",
                     num_active_upper_var);
      jnlst_->Printf(J_SUMMARY, J_MAIN,
                     "  Number of inequalities active at lower bound...: %d\n",
                     num_active_lower_con);
      jnlst_->Printf(J_SUMMARY, J_MAIN,
                     "  Number of inequalities active at upper bound...: %d\n",
                     num_active_upper_con);
      jnlst_->Printf(J_SUMMARY, J_MAIN,
                     "  Number of active equalities ...................: %d\n",
                     num_equality);
    }

  } else {
    jnlst_->Printf(J_SUMMARY, J_MAIN,
                   "\nUser did not provide an initial working set.\n");
  }
  jnlst_->Printf(J_SUMMARY, J_MAIN, "\n");

  // Initalize algorithmic quantities
  trust_region_radius_ = trust_region_init_size_;
  current_penalty_parameter_ = penalty_parameter_init_value_;

  if (disable_trust_region_) {
    trust_region_radius_ = 1e3;
  }
}

void SqpSolver::compute_initial_values_()
{
  // Determine the bound values
  sqp_nlp_->get_bounds_info(lower_variable_bounds_, upper_variable_bounds_,
                            lower_constraint_bounds_, upper_constraint_bounds_);

  // shift starting point to satisfy the bound constraint
  if (!slack_formulation_) {
    shift_starting_point_(current_iterate_, lower_variable_bounds_,
                          upper_variable_bounds_);
  }

  // Compute function and derivative values at the starting point
  bool retval = eval_f_(current_iterate_, current_objective_value_);
  assert(retval && "eval_f returned false at starting point");
  retval = eval_gradient_(current_iterate_, current_objective_gradient_);
  assert(retval && "eval_gradient returned false");
  retval = eval_constraints_(current_iterate_, current_constraint_values_);
  assert(retval && "eval_constraints returned false");
  retval =
      get_hessian_structure_(current_iterate_, current_constraint_multipliers_,
                             current_lagrangian_hessian_);
  assert(retval && "get_hessian_structure returned false");
  retval = eval_hessian_(current_iterate_, current_constraint_multipliers_,
                         current_lagrangian_hessian_);
  assert(retval && "eval_hessian returned false");
  retval =
      get_jacobian_structure_(current_iterate_, current_constraint_jacobian_);
  assert(retval && "get_jacobian_structure returned false");
  retval = eval_jacobian_(current_iterate_, current_constraint_jacobian_);
  assert(retval && "eval_jacobian returned false");

  current_infeasibility_ = compute_constraint_violation_(
      current_iterate_, current_constraint_values_);

  // Initialize the flag that tracks the status of the watchdog procedure
  if (watchdog_min_wait_iterations_ == 0) {
    watchdog_status_ = WATCHDOG_INACTIVE;
    watchdog_sleep_iterations_ = watchdog_min_wait_iterations_;
  }
  else {
    watchdog_status_ = WATCHDOG_READY;
    watchdog_sleep_iterations_ = 0;
  }

  // We need to make sure that complementarity holds for the multipliers
  for (int i = 0; i < num_variables_; ++i) {
    if (current_iterate_->get_value(i) > lower_variable_bounds_->get_value(i)) {
      current_bound_multipliers_->set_value(
          i, min(0., current_bound_multipliers_->get_value(i)));
    }
    if (current_iterate_->get_value(i) < upper_variable_bounds_->get_value(i)) {
      current_bound_multipliers_->set_value(
          i, max(0., current_bound_multipliers_->get_value(i)));
    }
  }
  for (int i = 0; i < num_constraints_; ++i) {
    if (current_constraint_values_->get_value(i) >
        lower_constraint_bounds_->get_value(i)) {
      current_constraint_multipliers_->set_value(
          i, min(0., current_constraint_multipliers_->get_value(i)));
    }
    if (current_constraint_values_->get_value(i) <
        upper_constraint_bounds_->get_value(i)) {
      current_constraint_multipliers_->set_value(
          i, max(0., current_constraint_multipliers_->get_value(i)));
    }
  }

}
void SqpSolver::initialize_qp_solvers_()
{
  // Determine the problem size
  shared_ptr<const SqpNlpSizeInfo> nlp_sizes = sqp_nlp_->get_problem_sizes();
  int num_variables = nlp_sizes->get_num_variables();
  int num_constraints = nlp_sizes->get_num_constraints();

  if (num_variables != num_variables_ ||
      num_constraints != num_constraints_) {
    jnlst_->Printf(J_ERROR, J_MAIN, "Error: Size of NLP has changed between restart calls.\n");
    THROW_EXCEPTION(SQP_INVALID_STARING_POINT, "Size of NLP has changed between restart calls");
  }

  // Create (new) matrix objects
  current_constraint_jacobian_ =
      make_shared<SparseTripletMatrix>(nlp_sizes->get_num_nonzeros_jacobian(),
                                       num_constraints_, num_variables_, false);
  current_lagrangian_hessian_ =
      make_shared<SparseTripletMatrix>(nlp_sizes->get_num_nonzeros_hessian(),
                                       num_variables_, num_variables_, true);

  // Create the (new) QP solver objects.  This also reads the option values
  qp_solver_ =
      make_shared<QpHandler>(nlp_sizes, QP_TYPE_QP, slack_formulation_,
                             sqp_nlp_->get_nlp_name(), jnlst_, options_);
  lp_solver_ =
      make_shared<QpHandler>(nlp_sizes, QP_TYPE_LP, slack_formulation_,
                             sqp_nlp_->get_nlp_name(), jnlst_, options_);

  // For the LP, the objective will never change.  We set the NLP gradient to
  // zero, and sjust sum up all slack variables by setting the penalty parameter
  // to 1.
  lp_solver_->set_linear_qp_objective_coefficients_to_zero();
  lp_solver_->update_penalty_parameter(1.);

  // Set flag that makes sure that the QP solver will be initialized at the
  // first call
  qp_solver_initialized_ = false;

  // If a initial working set was given, provide it now to the QP solver
  if (init_bound_activities_) {
    assert(init_constraint_activities_);

    qp_solver_->set_initial_working_sets(init_bound_activities_,
                                         init_constraint_activities_);

    // Delete the initial working set so that it won't be used again later
    delete[] init_bound_activities_;
    init_bound_activities_ = nullptr;
    delete[] init_constraint_activities_;
    init_constraint_activities_ = nullptr;
  }
}

/**
 * @brief alloocate memory for class members.
 * This function initializes all the shared pointer which will be used in the
 * Algorithm::Optimize, and it copies all parameters that might be changed
 * during
 * the run of the function Algorithm::Optimize.
 *
 * @param nlp: the nlp reader that read data of the function to be minimized;
 */
void SqpSolver::allocate_memory_()
{
  // Determine the problem size
  shared_ptr<const SqpNlpSizeInfo> nlp_sizes = sqp_nlp_->get_problem_sizes();
  num_variables_ = nlp_sizes->get_num_variables();
  num_constraints_ = nlp_sizes->get_num_constraints();

  // Allocate memory for the arrays that hold the constraint types
  constraint_type_ = new ConstraintType[num_constraints_];
  bound_type_ = new ConstraintType[num_variables_];

  // Create Vector objects that hold the bounds
  lower_variable_bounds_ = make_shared<Vector>(num_variables_);
  upper_variable_bounds_ = make_shared<Vector>(num_variables_);
  lower_constraint_bounds_ = make_shared<Vector>(num_constraints_);
  upper_constraint_bounds_ = make_shared<Vector>(num_constraints_);

  // Create Vector objects that hold the iterate values
  current_iterate_ = make_shared<Vector>(num_variables_);
  current_constraint_multipliers_ = make_shared<Vector>(num_constraints_);
  current_bound_multipliers_ = make_shared<Vector>(num_variables_);

  trial_iterate_ = make_shared<Vector>(num_variables_);
  trial_step_ = make_shared<Vector>(num_variables_);
  trial_constraint_multipliers_ = make_shared<Vector>(num_constraints_);
  trial_bound_multipliers_ = make_shared<Vector>(num_variables_);

  // Create objects that hold evaluated quantities
  current_constraint_values_ = make_shared<Vector>(num_constraints_);
  trial_constraint_values_ = make_shared<Vector>(num_constraints_);
  current_objective_gradient_ = make_shared<Vector>(num_variables_);
}

/**
* @brief This function calculates the infeasibility for given x_k and c_k with
respect
* to their corresponding bounds
* @return current_infeasibility_ = ||-max(c_k-c_u),0||_1 +||-min(c_k-c_l),0||_1+
                  ||-max(x_k-upper_variable_bounds_),0||_1
+||-min(x_k-lower_variable_bounds_),0||_1
*/
double SqpSolver::compute_constraint_violation_(
    shared_ptr<const Vector> variable_values,
    shared_ptr<const Vector> constraint_values)
{
  double current_infeasibility_ = 0.0;
  for (int i = 0; i < constraint_values->get_dim(); i++) {
    if (constraint_values->get_value(i) <
        lower_constraint_bounds_->get_value(i))
      current_infeasibility_ += (lower_constraint_bounds_->get_value(i) -
                                 constraint_values->get_value(i));
    else if (constraint_values->get_value(i) >
             upper_constraint_bounds_->get_value(i))
      current_infeasibility_ += (constraint_values->get_value(i) -
                                 upper_constraint_bounds_->get_value(i));
  }

  // TODO: I don't think we need to compute the following in the non-slack
  // formulation
  if (variable_values != nullptr) {
    for (int i = 0; i < current_iterate_->get_dim(); i++) {
      if (variable_values->get_value(i) < lower_variable_bounds_->get_value(i))
        current_infeasibility_ += (lower_variable_bounds_->get_value(i) -
                                   variable_values->get_value(i));
      else if (variable_values->get_value(i) >
               upper_variable_bounds_->get_value(i))
        current_infeasibility_ += (variable_values->get_value(i) -
                                   upper_variable_bounds_->get_value(i));
    }
  }

  return current_infeasibility_;
}

double SqpSolver::calc_model_infeasibility_(shared_ptr<const Vector> step)
{
  double infeas = 0.;
  // Compute the value of the linearized constraint body
  if (num_constraints_ > 0) {
    // Compute Jacobian vector product with step
    shared_ptr<Vector> constraint_body = make_shared<Vector>(num_constraints_);

    constraint_body->copy_vector(current_constraint_values_);
    current_constraint_jacobian_->multiply(step, constraint_body);

    const double* c = constraint_body->get_values();
    const double* c_L = lower_constraint_bounds_->get_values();
    const double* c_U = upper_constraint_bounds_->get_values();

    // Add up violation of the individual constraints
    for (int i = 0; i < num_constraints_; ++i) {
      infeas += max(0., c_L[i] - c[i]);
      infeas += max(0., c[i] - c_U[i]);
    }
  }

  // Add the violation of the bound constraints, in case we deal with the slack
  // formulation
  if (slack_formulation_) {
    const double* x_k = current_iterate_->get_values();
    const double* p = step->get_values();
    const double* x_L = lower_variable_bounds_->get_values();
    const double* x_U = upper_variable_bounds_->get_values();

    // Add up violation of the individual bounds
    for (int i = 0; i < num_variables_; ++i) {
      double x_p = x_k[i] + p[i];
      infeas += max(0., x_L[i] - x_p);
      infeas += max(0., x_p - x_U[i]);
    }
  }

  return infeas;
}

/**
 * @brief This function will set up the data for the QP subproblem
 *
 * It will initialize all the data at once at the beginning of the Algorithm.
 * After
 * that, the data in the QP problem will be updated according to the class
 * member QPinfoFlag_
 */

void SqpSolver::setup_qp_(shared_ptr<const Vector> qp_constraint_body)
{
  /** If this is the first QP that needs to be solved, set the initial values,
   * including matrix structures. */
  if (!qp_solver_initialized_) {
    // Set all QP quantities for the first time
    qp_solver_->set_jacobian(current_constraint_jacobian_);
    qp_solver_->set_hessian(current_lagrangian_hessian_);
    qp_solver_->set_bounds(trust_region_radius_, lower_variable_bounds_,
                           upper_variable_bounds_, current_iterate_,
                           lower_constraint_bounds_, upper_constraint_bounds_,
                           qp_constraint_body);
    qp_solver_->set_linear_qp_objective_coefficients(
        current_objective_gradient_, current_penalty_parameter_);

    // Initialize the update tracker
    qp_update_tracker_.reset();

    // Nothing more needs to be done.
    qp_solver_initialized_ = true;

    // Store the value of the penalty parameter that was just used so we know it
    // next time to decide if we need to update it.
    last_qp_penalty_parameter_ = current_penalty_parameter_;

    return;
  }

  // Check if the penalty parameter value has changed
  if (current_penalty_parameter_ != last_qp_penalty_parameter_) {
    qp_update_tracker_.trigger_penalty_parameter_update();
  }

  // Debug check: Make sure that an update is necessary.
  assert(qp_update_tracker_.need_update() && "QP is not changed.");

  // Update those quantities that need to be updated
  if (qp_update_tracker_.need_jacobian_update()) {
    qp_solver_->set_jacobian(current_constraint_jacobian_);
  }
  if (qp_update_tracker_.need_hessian_update()) {
    qp_solver_->set_hessian(current_lagrangian_hessian_);
  }

  if (qp_update_tracker_.need_bounds_update()) {
    qp_solver_->set_bounds(trust_region_radius_, lower_variable_bounds_,
                           upper_variable_bounds_, current_iterate_,
                           lower_constraint_bounds_, upper_constraint_bounds_,
                           qp_constraint_body);
  } else if (qp_update_tracker_.need_trust_region_radius_decrease()) {
    qp_solver_->decrease_trust_region(trust_region_radius_);
  }

  if (qp_update_tracker_.need_gradient_update()) {
    qp_solver_->set_linear_qp_objective_coefficients(
        current_objective_gradient_, current_penalty_parameter_);
  } else if (qp_update_tracker_.need_penalty_parameter_update()) {
    qp_solver_->update_penalty_parameter(current_penalty_parameter_);
  }

  // Reset the update tracker
  qp_update_tracker_.reset();

  // Store the value of the penalty parameter that was just used so we know it
  // next time to decide if we need to update it.
  last_qp_penalty_parameter_ = current_penalty_parameter_;
}

double SqpSolver::numerical_error_buffer_()
{
  // For round-off error, increase the predicted reduction by a tiny
  // amount
  // TODO: Make these options
  double increase_factor = 1e-10;
  double increase =
      increase_factor *
      max(1., max(fabs(current_objective_value_), current_infeasibility_));
  return increase;
}

void SqpSolver::setup_lp_()
{
  // The objectives does not change, it is always the sum of the slack
  // variables, and has been set in the initialization method.

  // All we have to do is update the constraints (bounds and Jacobian)
  lp_solver_->set_bounds(trust_region_radius_, lower_variable_bounds_,
                         upper_variable_bounds_, current_iterate_,
                         lower_constraint_bounds_, upper_constraint_bounds_,
                         current_constraint_values_);
  lp_solver_->set_jacobian(current_constraint_jacobian_);
}

void SqpSolver::increase_penalty_parameter_()
{
  if (current_penalty_parameter_ >= penalty_parameter_max_value_) {
    jnlst_->Printf(J_ERROR, J_MAIN, "Penalty parameter becomes too large: %e\n",
                   current_penalty_parameter_);
    THROW_EXCEPTION(SQP_EXCEPTION_PENALTY_TOO_LARGE,
                    "Penalty parameter becomes too large.")
  }

  // Increase the penalty parameter
  current_penalty_parameter_ =
      min(penalty_parameter_max_value_,
          current_penalty_parameter_ * penalty_parameter_increase_factor_);

  // Increase counter for penalty parameter trials
  solver_statistics_->try_new_penalty_parameter();

  jnlst_->Printf(J_ITERSUMMARY, J_MAIN, "  Solving QP for penalty %e\n",
                 current_penalty_parameter_);

  // Solve the QP for the new penalty parameter
  shared_ptr<const Vector> qp_constraint_body = current_constraint_values_;
  calculate_search_direction_(qp_constraint_body);
}

void SqpSolver::update_predicted_reduction_()
{
  // Start with the gradient part
  predicted_reduction_ =
      -current_objective_gradient_->calc_inner_product(trial_step_);

  // Add the quadratic term
  if (current_lagrangian_hessian_) {
    shared_ptr<Vector> Hx = make_shared<Vector>(num_variables_);
    Hx->set_to_zero();
    current_lagrangian_hessian_->multiply(trial_step_, Hx);
    predicted_reduction_ -= 0.5 * Hx->calc_inner_product(trial_step_);
  }

  // Finally add the infeasibility part
  predicted_reduction_ += current_penalty_parameter_ *
                          (current_infeasibility_ - trial_model_infeasibility_);
  jnlst_->Printf(J_DETAILED, J_MAIN, "Computing predicted reduction: %23.16e\n", predicted_reduction_);

  predicted_reduction_ += numerical_error_buffer_();
}

/**
 *
 * @brief This function performs the ratio test to determine if we should accept
 * the trial point
 *
 * The ratio is calculated by
 * (P_1(x_k;\rho)-P_1( x_trial;\rho))/(q_k(0;\rho)-q_k(p_k;rho), where
 * P_1(x,rho) = f(x) + rho* infeasibility_measure is the l_1 merit function and
 * q_k(p; rho) = f_k+ g_k^Tp +1/2 p^T H_k p+rho* infeasibility_measure_model is
 * the
 * quadratic model at x_k.
 * The trial point  will be accepted if the ratio >=
 * trust_region_ratio_accept_tol.
 * If it is accepted, the function will also updates the gradient, Jacobian
 * information by reading from nlp_ object. The corresponding flags of class
 * member
 * QPinfoFlag_ will set to be true.
 */
void SqpSolver::perform_ratio_test_()
{
  // Compute the actual reduction in the merit function
  double current_penalty_function_value =
      current_objective_value_ +
      current_penalty_parameter_ * current_infeasibility_;

  double trial_penalty_function_value =
      trial_objective_value_ +
      current_penalty_parameter_ * trial_infeasibility_;

  jnlst_->Printf(J_DETAILED, J_MAIN, "\nIn ratio test: current penalty function value..: %23.16e\n", current_penalty_function_value);
  jnlst_->Printf(J_DETAILED, J_MAIN, "               trial penalty function value....: %23.16e\n", trial_penalty_function_value);
  jnlst_->Printf(J_DETAILED, J_MAIN, "               current objective value.........: %23.16e\n", current_objective_value_);
  jnlst_->Printf(J_DETAILED, J_MAIN, "               trial objective value...........: %23.16e\n", trial_objective_value_);
  jnlst_->Printf(J_DETAILED, J_MAIN, "               current infeasibility...........: %23.16e\n", current_infeasibility_);
  jnlst_->Printf(J_DETAILED, J_MAIN, "               trial infeasibility.............: %23.16e\n", trial_infeasibility_);
  jnlst_->Printf(J_DETAILED, J_MAIN, "               current penalty parameter.......: %23.16e\n\n", current_penalty_parameter_);

  if (disable_trust_region_) {
    trial_point_is_accepted_ = true;
    predicted_reduction_ = 1.;
    actual_reduction_ = 1.;
    trial_model_infeasibility_ = 0.;
    return;
  }

  // The predicted reduction has already been computed for regular iterations.
  // If we are in the trial watchdog iteration, we need to use the one from the
  // original iteration in which the watchdog was triggered
  if (watchdog_status_ == WATCHDOG_IN_TRIAL_ITERATE) {
    predicted_reduction_ = backup_predicted_reduction_;
  }

  // If we are in the trial watchdog iteration, we need to use the original value
  // for the current penatly value
  if (watchdog_status_ == WATCHDOG_IN_TRIAL_ITERATE) {
    current_penalty_function_value =
          backup_current_objective_value_ +
          current_penalty_parameter_ * backup_current_infeasibility_;
    jnlst_->Printf(J_DETAILED, J_MAIN, "\nIn watchdog  : current penalty function value..: %23.16e\n", current_penalty_function_value);
    jnlst_->Printf(J_DETAILED, J_MAIN, "               trial penalty function value....: %23.16e\n", trial_penalty_function_value);
    jnlst_->Printf(J_DETAILED, J_MAIN, "               backup objective value..........: %23.16e\n", backup_current_objective_value_);
    jnlst_->Printf(J_DETAILED, J_MAIN, "               trial objective value...........: %23.16e\n", trial_objective_value_);
    jnlst_->Printf(J_DETAILED, J_MAIN, "               backup infeasibility............: %23.16e\n", backup_current_infeasibility_);
    jnlst_->Printf(J_DETAILED, J_MAIN, "               trial infeasibility.............: %23.16e\n", trial_infeasibility_);
    jnlst_->Printf(J_DETAILED, J_MAIN, "               current penalty parameter.......: %23.16e\n\n", current_penalty_parameter_);
  }

  actual_reduction_ =
      current_penalty_function_value - trial_penalty_function_value;

  // For round-off error, increase the actual reduction by a tiny
  // amount
  actual_reduction_ += numerical_error_buffer_();

  jnlst_->Printf(J_DETAILED, J_MAIN, "\nRatio test:\n");
  jnlst_->Printf(J_DETAILED, J_MAIN, "   predicted reduction: %9.2e\n",
                 predicted_reduction_);
  jnlst_->Printf(J_DETAILED, J_MAIN, "      actual reduction: %9.2e\n\n",
                 actual_reduction_);

  if (predicted_reduction_ <= 0.) {
    //    myQP_->WriteQPData(problem_name_+"qpdata.log");
    jnlst_->Printf(J_ERROR, J_MAIN, "Predicted reduction is non-positive: %e\n",
                   predicted_reduction_);
    exit_flag_ = PRED_REDUCTION_NEGATIVE;
    return;
  }

  if (actual_reduction_ >=
      (trust_region_ratio_accept_tol_ * predicted_reduction_)) {
    trial_point_is_accepted_ = true; // no need to calculate the SOC direction
  } else {
    trial_point_is_accepted_ = false;
  }
}

void SqpSolver::accept_trial_point_()
{
  assert(trial_point_is_accepted_);

  // succesfully update
  // copy information already calculated from the trial point
  current_infeasibility_ = trial_infeasibility_;
  current_objective_value_ = trial_objective_value_;
  current_iterate_->copy_vector(trial_iterate_);
  current_constraint_values_->copy_vector(trial_constraint_values_);

  // Get the new multiplier iterate from the QP solution
  current_constraint_multipliers_->copy_vector(trial_constraint_multipliers_);
  current_bound_multipliers_->copy_vector(trial_bound_multipliers_);

  // Evaluate the derivatives at the new iterate
  bool retval = eval_gradient_(current_iterate_, current_objective_gradient_);
  assert(retval && "eval_gradient returned false");
  retval = eval_jacobian_(current_iterate_, current_constraint_jacobian_);
  assert(retval && "eval_jacobian returned false");
  retval = eval_hessian_(current_iterate_, current_constraint_multipliers_,
                         current_lagrangian_hessian_);
  assert(retval && "eval_hessian returned false");

  // Tell QP solver that all data has changed
  qp_update_tracker_.trigger_all_updates();
}

/**
 * @brief Update the trust region radius.
 *
 * This function update the trust-region radius when the ratio calculated by the
 * ratio test is smaller than eta_c or bigger than
 * trust_region_ratio_increase_tol and the
 * search_direction
 * hits the trust-region bounds.
 * If ratio<eta_c, the trust region radius will decrease by the parameter
 * gamma_c, to be gamma_c* delta_
 * If ratio_test> trust_region_ratio_increase_tol and delta_= norm_p_k_ the
 * trust-region radius will be
 * increased by the parameter gamma_c.
 *
 * If trust region radius has changed, the corresponding flags will be set to be
 * true;
 */

void SqpSolver::update_trust_region_radius_()
{
  if (disable_trust_region_) {
    return;
  }
  if (actual_reduction_ <
      trust_region_ratio_decrease_tol_ * predicted_reduction_) {

    // The actual reduction is not sufficiently large, so let's decrease the
    // trust region radius
    double trial_step_size = trial_step_->calc_inf_norm();
    trust_region_radius_ = trust_region_decrease_factor_ *
                           min(trust_region_radius_, trial_step_size);

    // Tell the QP solver that the radius was decreased
    qp_update_tracker_.trigger_trust_region_radius_decrease();

  } else {
    // Check if actual reduction is sufficiently large.
    if (actual_reduction_ >
        trust_region_ratio_increase_tol_ * predicted_reduction_) {
      // If so, we increase the trust region if it was active
      // TODO: Find good replacement for opt_tol_ here
      if (fabs(trust_region_radius_ - trial_step_->calc_inf_norm()) <
          opt_tol_) {
        trust_region_radius_ =
            min(trust_region_increase_factor_ * trust_region_radius_,
                trust_region_max_value_);
        // Tell the QP solver that the variable bounds need to be updated
        qp_update_tracker_.trigger_bounds_update();
      }
    }
  }

  // if the trust-region becomes too small, throw the error message
  if (trust_region_radius_ < trust_region_min_value_) {
    jnlst_->Printf(J_ERROR, J_MAIN,
                   "Trust region radius is becoming too small: %e\n",
                   trust_region_radius_);
    exit_flag_ = TRUST_REGION_TOO_SMALL;
    THROW_EXCEPTION(SMALL_TRUST_REGION, SMALL_TRUST_REGION_MSG);
  }
}

/**
 *
 * @brief This function checks how each constraint specified by the nlp readers
 * are
 * bounded.
 * If there is only upper bounds for a constraint, c_i(x)<=c^i_u, then
 * cons_type_[i]= BOUNDED_ABOVE
 * If there is only lower bounds for a constraint, c_i(x)>=c^i_l, then
 * cons_type_[i]= BOUNDED_BELOW
 * If there are both upper bounds and lower bounds, c^i_l<=c_i(x)<=c^i_u, and
 * c^i_l<c^i_u then cons_type_[i]= BOUNDED,
 * If there is no constraints on all
 * of c_i(x), then cons_type_[i]= UNBOUNDED;
 *
 * The same rules are also applied to the bound-constraints.
 */

static ConstraintType classify_single_constraint(double lower_bound,
                                                 double upper_bound)
{
  assert(lower_bound <= upper_bound);
  if (lower_bound > -SqpInf && upper_bound < SqpInf) {
    if (upper_bound == lower_bound) {
      return IS_EQUALITY;
    } else {
      return BOUNDED_BELOW_AND_ABOVE;
    }
  } else if (lower_bound > -SqpInf && upper_bound >= SqpInf) {
    return BOUNDED_BELOW;
  } else if (upper_bound < SqpInf && lower_bound <= -SqpInf) {
    return BOUNDED_ABOVE;
  } else {
    return UNBOUNDED;
  }
}

void SqpSolver::classify_constraints_types_()
{
  // Determine which of the variables have lower and/or upper bounds
  for (int i = 0; i < num_variables_; i++) {
    bound_type_[i] =
        classify_single_constraint(lower_variable_bounds_->get_value(i),
                                   upper_variable_bounds_->get_value(i));
  }
  // Determine which of the constraints have lower and/or upper bounds
  // At the same time, count the equality constraints;
  num_equality_constraints_ = 0;
  for (int i = 0; i < num_constraints_; i++) {
    constraint_type_[i] =
        classify_single_constraint(lower_constraint_bounds_->get_value(i),
                                   upper_constraint_bounds_->get_value(i));
    if (constraint_type_[i] == IS_EQUALITY) {
      num_equality_constraints_++;
    }
  }
}

/**
 * @brief update the penalty parameter for the algorithm.
 *
 */
void SqpSolver::update_penalty_parameter_()
{
  if (disable_trust_region_) {
    return;
  }

  jnlst_->Printf(J_DETAILED, J_MAIN,
                 "\nCheck if the penalty parameter (%e) has to be increased.\n",
                 current_penalty_parameter_);

  // Consider the violation of the linearized cnostraints at the trial step
  jnlst_->Printf(J_DETAILED, J_MAIN, "  trial model infeasibility: %e\n",
                 trial_model_infeasibility_);

  // If there is no violation of the linear constraints, we do not need to
  // compare with feasibility progress from LP
  if (trial_model_infeasibility_ <= penalty_update_tol_) {
    jnlst_->Printf(J_DETAILED, J_MAIN, "    This is below the tolerance (%e), "
                                       "so no solution of the feasibility LP "
                                       "necessary.\n",
                   penalty_update_tol_);
  } else {
    jnlst_->Printf(J_DETAILED, J_MAIN, "     This is larger than the tolerance "
                                       "(%e), so we need to determine best "
                                       "possible improvement of feasibility.\n",
                   penalty_update_tol_);
    jnlst_->Printf(J_ITERSUMMARY, J_MAIN, "  Solving LP for penalty %e\n",
                   current_penalty_parameter_);

    // Solve the LP that minimizes the linearized constraint vilation within the
    // trust region
    setup_lp_();
    QpSolverExitStatus exit_status = lp_solver_->solve(solver_statistics_);
    if (exit_status != QPEXIT_OPTIMAL) {
      const string& nlp_name = sqp_nlp_->get_nlp_name();
      //lp_solver_->write_qp_data(nlp_name + "qpdata.log");
      assert(false && "Still need to decide how to handle QP solver error.");
      // exit_flag_ = qp_solver_->get_status();
      // break;
    }

    // Extract the step from the LP solution
    shared_ptr<Vector> lp_step = make_shared<Vector>(num_variables_);
    lp_step->copy_values(lp_solver_->get_primal_solution()->get_values());

    // Compute the violation of the linearized constraints
    double lp_model_infeasibility = calc_model_infeasibility_(lp_step);

    jnlst_->Printf(J_DETAILED, J_MAIN, "  lp_model_infeasibility: %e\n",
                   lp_model_infeasibility);

    // If the LP solution has no constraint violation, we want to increase the
    // penalty parameter so large that the QP solution also has no violation of
    // the linearlized constraints
    if (lp_model_infeasibility <= penalty_update_tol_) {

      jnlst_->Printf(J_DETAILED, J_MAIN,
                     "    This is below the tolerance (%e), "
                     "so we need to increase the penalty "
                     "parameter until the trial model "
                     "infeasibility is (almost) zero\n",
                     penalty_update_tol_);

      // try to increase the penalty parameter to a number such that the
      // infeasibility measure of QP model with such penalty parameter
      // becomes zero
      while (trial_model_infeasibility_ > penalty_update_tol_) {

        // increase the penalty parameter and recompute the trial step and
        // the corresponding violation of the linearized constraints
        increase_penalty_parameter_();

        jnlst_->Printf(
            J_DETAILED, J_MAIN,
            "    new penalty parameter: %e with trial model infeasibility %e\n",
            current_penalty_parameter_, trial_model_infeasibility_);
      }
    } else {

      // try to increase the penalty parameter to a number such that
      // the incurred reduction for the QP model is to a ratio to the
      // maximum possible reduction for current linear model.
      double lp_infeasibility_reduction =
          current_infeasibility_ - lp_model_infeasibility;
      double candidate_infeasibility_reduction =
          current_infeasibility_ - trial_model_infeasibility_;

      jnlst_->Printf(J_DETAILED, J_MAIN,
                     "    This is above the tolerance (%e), so we need to "
                     "increase the penalty parameter until the trial model "
                     "infeasibility is smaller than a fraction of the LP model "
                     "infeasibility (%e).\n",
                     penalty_update_tol_, eps1_ * lp_infeasibility_reduction);

      jnlst_->Printf(J_DETAILED, J_MAIN,
                     "  candidate infeasibility reduction: %e\n",
                     candidate_infeasibility_reduction);
      while (candidate_infeasibility_reduction <
             eps1_ * lp_infeasibility_reduction) {
        // increase the penalty parameter and recompute the trial step and
        // the corresponding violation of the linearized constraints
        increase_penalty_parameter_();

        candidate_infeasibility_reduction =
            current_infeasibility_ - trial_model_infeasibility_;

        jnlst_->Printf(J_DETAILED, J_MAIN, "  new penalty parameter %e with "
                                           "infeasibility reduction: %e\n",
                       current_penalty_parameter_,
                       candidate_infeasibility_reduction);
      }
    }
  }

  // New we need to make sure that the predicted reduction is at least a
  // fraction of the predicted infeasibility reduction

  // Compute the reduction in the model of the penalty function
  update_predicted_reduction_();

  // Compute the reduction in the violation of the linearized constraints
  double predicted_infeasibility_reduction =
      current_infeasibility_ - trial_model_infeasibility_;

  // Increase this a bit to take care of potential numerical error
  predicted_infeasibility_reduction += numerical_error_buffer_();

  jnlst_->Printf(J_DETAILED, J_MAIN, "  Make sure that predicted reduction >= "
                                     "(%e) * predicted infeasibility "
                                     "reduction:\n",
                 eps2_);
  jnlst_->Printf(J_DETAILED, J_MAIN,
                 "    pred reduction = %e, pred infeas reduction = %e\n",
                 predicted_reduction_, predicted_infeasibility_reduction);

  // Increase the penalty parameter until condition is satisfied
  while (predicted_reduction_ < eps2_ * current_penalty_parameter_ *
                                    predicted_infeasibility_reduction) {

    jnlst_->Printf(J_DETAILED, J_MAIN, "    Test not satisfied (%23.16e < "
                                       "%23.16e), increase penalty "
                                       "parameter.\n",
                   predicted_reduction_, eps2_ * current_penalty_parameter_ *
                                             predicted_infeasibility_reduction);
    increase_penalty_parameter_();

    // Update the predicted reduction to match the new penalty parameter
    update_penalty_parameter_();

    // Update reduction in the violation of the linearized constraints
    predicted_infeasibility_reduction =
        current_infeasibility_ - trial_model_infeasibility_;
  }
  jnlst_->Printf(J_DETAILED, J_MAIN,
                 "   Penalty parameter now large enough.\n");
}

/**
 * @brief Use the Ipopt Reference Options and set it to default values.
 */
void SqpSolver::register_options_(SmartPtr<RegisteredOptions> reg_options,
                                  bool skip_ipopt_options)
{
  // Options related to output
  //
  // Skip this if Ipopt has already registered these options
  if (!skip_ipopt_options) {
    reg_options->SetRegisteringCategory("Output");
    reg_options->AddBoundedIntegerOption(
        "print_level", "Output verbosity level.", 0, J_LAST_LEVEL - 1,
        J_ITERSUMMARY,
        "Sets the default verbosity level for console output. The "
        "larger this value the more detailed is the output.");

    reg_options->AddStringOption1(
        "output_file",
        "File name of desired output file (leave unset for no file output).",
        "", "*", "Any acceptable standard file name",
        "NOTE: This option only works when read from the sqp.opt options file! "
        "An output file with this name will be written (leave unset for no "
        "file output).  The verbosity level is by default set to "
        "\"print_level\", "
        "but can be overridden with \"file_print_level\".  The file name is "
        "changed to use only small letters.");
    reg_options->AddBoundedIntegerOption(
        "file_print_level", "Verbosity level for output file.", 0,
        J_LAST_LEVEL - 1, J_ITERSUMMARY,
        "NOTE: This option only works when read from the sqp.opt options file! "
        "Determines the verbosity level for the file specified by "
        "\"output_file\".  By default it is the same as \"print_level\".");
  }

  reg_options->SetRegisteringCategory("Starting Point");
  reg_options->AddStringOption3(
      "starting_mode", "Specifies how much starting information is available.",
      "primal-dual",
      "primal", "only primal point is provided",
      "primal-dual", "primal and dual variables are provided.",
      "warm-start", "primal-dual starting point and intial working set is provided.");

  reg_options->SetRegisteringCategory("Trust-region");
  reg_options->AddBoundedNumberOption(
      "trust_region_ratio_decrease_tol",
      "trust-region parameter for the ratio test triggering decrease.", 0.,
      true, 1., true, 1e-8,
      "If ratio <= trust_region_ratio_decrease_tol, then the trust-region "
      "radius for the next "
      "iteration will be decreased for the next iteration.  This must not be "
      "smaller than trust_region_ratio_accept_tol.");
  reg_options->AddBoundedNumberOption(
      "trust_region_ratio_accept_tol",
      "trust-region parameter for the ratio test.", 0., true, 1., true, 1e-8,
      "The trial point will be accepted if ratio >= "
      "trust_region_ratio_accept_tol. ");
  reg_options->AddBoundedNumberOption(
      "trust_region_ratio_increase_tol",
      "trust-region parameter for the ratio test.", 0., true, 1., true, 1e-8,
      "If ratio >= trust_region_ratio_increase_tol and the search direction "
      "hits the  "
      "trust-region boundary, the trust-region radius will "
      "be increased for the next iteration.  This must not be smaller than "
      "trust_region_ratio_decrease_tol.");
  reg_options->AddBoundedNumberOption(
      "trust_region_decrease_factor",
      "Factor used to reduce the trust-region size.", 0., true, 1., true, 0.5,
      "If the trust-region radius is going to be decreased, "
      "then it will be multiplied by the value of this options.");
  reg_options->AddLowerBoundedNumberOption(
      "trust_region_increase_factor",
      "Factor used to increase the trust-region size.", 1., true, 2.0,
      "If the trust-region radius is going to be "
      "increased, then it will be set as gamma_e*delta,"
      "where delta is current trust-region radius.");

  reg_options->AddLowerBoundedNumberOption("trust_region_init_size",
                                           "Initial trust-region radius",
                                           0., true, 10.);
  reg_options->AddLowerBoundedNumberOption(
      "trust_region_max_value", "Maximum value of trust-region radius "
                                "allowed for the radius update",
      0., true, 1e10);
  reg_options->AddLowerBoundedNumberOption(
      "trust_region_min_value", "Minimum value of trust-region radius "
                                "allowed for the radius update",
      0., true, 1e-16);
  reg_options->AddStringOption2(
      "disable_trust_region", "Indicates whether trust region should be used.",
      "no", "no", "do usual trust region method", "yes",
      "every trial step will be accepted");
  reg_options->AddLowerBoundedIntegerOption(
      "watchdog_min_wait_iterations", "Minimum number of watchdog wait iterations.",
      0, 10,
      "Specifies the number of iterations the watchdog should sleep after a rejected trial step.  0 switches the watchdog technique off.");

  reg_options->SetRegisteringCategory("Penalty Update");
  reg_options->AddLowerBoundedNumberOption(
      "penalty_parameter_init_value", "Initial value of the penalty parameter.",
      0., true, 10.0, "");
  reg_options->AddLowerBoundedNumberOption(
      "penalty_update_tol", "some tolerance.", 0., true, 1e-8, "");
  reg_options->AddLowerBoundedNumberOption(
      "penalty_parameter_increase_factor",
      "Factor by which penatly parameter is increased.", 1., true, 10., "");
  reg_options->AddNumberOption("eps1", "penalty update parameter something",
                               0.1, "");
  reg_options->AddNumberOption("eps1_change_parm",
                               "penalty update parameter something", 0.1, "");
  reg_options->AddNumberOption("eps2", "penalty update parameter something",
                               1.0e-6, "");
  reg_options->AddNumberOption("print_level_penalty_update",
                               "print level for penalty update", 0);
  reg_options->AddNumberOption("penalty_parameter_max_value",
                               "Maximum value of the penalty parameter", 1e12);
  reg_options->AddIntegerOption(
      "penalty_iter_max",
      "maximum number of penalty paramter update allowed in a "
      "single iteration in the main algorithm",
      200);
  reg_options->AddIntegerOption(
      "penalty_iter_max_total",
      "maximum number of penalty paramter update allowed "
      "in total",
      100);

  reg_options->SetRegisteringCategory("Optimality Test");
  reg_options->AddIntegerOption("testOption_NLP",
                                "Level of Optimality test for "
                                "NLP",
                                0);
  reg_options->AddStringOption2(
      "auto_gen_tol", "Tell the algorithm to automatically"
                      "generate the tolerance level for optimality test "
                      "based on information from NLP",
      "no", "no", "will use user-defined values of tolerance for"
                  " the optimality test",
      "yes", "will automatically generate the tolerance "
             "level for the optimality test");

  const double default_tol = 1e-6;
  reg_options->AddNumberOption("active_set_tol", "",
                               1.0e-5); // TODO: make lower bounded options
  reg_options->AddNumberOption("opt_tol", "", default_tol);
  reg_options->AddNumberOption("opt_tol_complementarity", "", default_tol);
  reg_options->AddNumberOption("opt_tol_dual_feasibility", " ", default_tol);
  reg_options->AddNumberOption("opt_tol_primal_feasibility", " ", default_tol);
  reg_options->AddNumberOption("opt_tol_stationarity_feasibility", "",
                               default_tol);
  reg_options->AddNumberOption("opt_second_tol", " ", default_tol);

  reg_options->AddLowerBoundedNumberOption(
      "cpu_time_limit", "CPU time limit", 0., true, 1e10,
      "Time limit measured in CPU time (in seconds)");
  reg_options->AddLowerBoundedNumberOption(
      "wallclock_time_limit", "Wallclock time limit", 0., true, 1e10,
      "Time limit measured in wallclock time (in seconds)");

  reg_options->SetRegisteringCategory("General");
  reg_options->AddNumberOption("objective_scaling_factor",
                               "scaling factor for the objective function", 1.,
                               "");
  reg_options->AddNumberOption("step_size_tol",
                               "the smallest stepsize can be accepted"
                               "before concluding convergence",
                               1.0e-15);
  reg_options->AddIntegerOption("max_num_iterations",
                                "Maximum number of iteration for the algorithm",
                                3000);
  reg_options->AddStringOption2(
      "perform_second_order_correction",
      "Tells the algorithm to calculate the second-order correction step "
      "during the main iteration",
      "no", "no", "not calculate the soc steps", "yes",
      "will calculate the soc steps");
  reg_options->AddStringOption2("slack_formulation",
                                "Indicates whether the SQP solver addresses "
                                "the formulation that introduces slack "
                                "variables for the variable bounds.",
                                "no", "no", "keep variables within bounds",
                                "yes", "permit bounds to be violated", "");

  reg_options->SetRegisteringCategory("QPsolver");
  reg_options->AddIntegerOption("testOption_QP",
                                "Level of Optimality test for QP", -99);
  reg_options->AddIntegerOption("qp_solver_max_num_iterations",
                                "maximum number of iteration for the "
                                "QP solver in solving each QP",
                                100000);
  reg_options->AddIntegerOption("lp_solver_max_num_iterations",
                                "maximum number of iteration for the "
                                "LP solver in solving each LP",
                                100000);
  reg_options->AddIntegerOption("qp_solver_print_level",
                                "print level for QP solver", 0);
  reg_options->AddStringOption2("qp_solver",
                                "QP solver used for step computation.", "qore",
                                "qpoases", "", "qore", "");

  //    reg_options->AddStringOption("QPsolverChoice",
  //		    "The choice of QP solver which will be used in the
  // Algorithm",
  //		    "qpOASES");

  reg_options->SetRegisteringCategory("LPsolver");
  //    reg_options->AddStringOption("LPsolverChoice",
  //		    "The choice of LP solver which will be used in the
  // Algorithm",
  //		    "qpOASES");

  reg_options->AddStringOption2(
      "qore_init_primal_variables", "Specifies whether QORE should initialize "
                                    "the primal soluiton (the step and the "
                                    "penatly variables) to zero.",
      "no", "no", "reuse the internal solution from the most recent solve.",
      "yes", "initialize the values to zero.");
  reg_options->AddLowerBoundedNumberOption(
      "qore_hessian_regularization", "Regularization parameter for the QP", 0.,
      false, 0., "Number that is added to the diagonal of the QP Hessian");
  reg_options->AddStringOption2(
      "qore_dump_file", "Indicates whether QORE should write a dump file.",
      "no", "no", "Do not write dump file.",
      "yes", "Write dump files (for QP and LP separately).");

  reg_options->AddIntegerOption("testOption_LP",
                                "Level of Optimality test for LP", -99);
  reg_options->AddNumberOption("iter_malower_variable_bounds_p",
                               "maximum number of iteration for the "
                               "LP solver in solving each LP",
                               100);
  reg_options->AddNumberOption("print_level_lp", "print level for LP solver",
                               0);
}

void SqpSolver::get_option_values_()
{
  string sval;
  options_->GetStringValue("starting_mode", sval, "");
  if (sval == "primal") {
    starting_mode_ = SM_PRIMAL_ONLY;
  }
  else if (sval == "primal-dual") {
    starting_mode_ = SM_PRIMAL_DUAL;
  }
  else if (sval == "warm-start") {
    starting_mode_ = SM_WARM_START;
  }
  else {
    assert(false && "Invalid string value for starting_mode option");
  }

  options_->GetIntegerValue("max_num_iterations", max_num_iterations_, "");
  options_->GetNumericValue("cpu_time_limit", cpu_time_limit_, "");
  options_->GetNumericValue("wallclock_time_limit", wallclock_time_limit_, "");

  options_->GetNumericValue("trust_region_init_size", trust_region_init_size_,
                            "");
  options_->GetNumericValue("trust_region_max_value", trust_region_max_value_,
                            "");
  options_->GetNumericValue("trust_region_min_value", trust_region_min_value_,
                            "");
  options_->GetNumericValue("trust_region_ratio_decrease_tol",
                            trust_region_ratio_decrease_tol_, "");
  options_->GetNumericValue("trust_region_ratio_accept_tol",
                            trust_region_ratio_accept_tol_, "");
  options_->GetNumericValue("trust_region_ratio_increase_tol",
                            trust_region_ratio_increase_tol_, "");
  options_->GetNumericValue("trust_region_decrease_factor",
                            trust_region_decrease_factor_, "");
  options_->GetNumericValue("trust_region_increase_factor",
                            trust_region_increase_factor_, "");
  options_->GetBoolValue("disable_trust_region", disable_trust_region_, "");
  options_->GetIntegerValue("watchdog_min_wait_iterations",
                             watchdog_min_wait_iterations_, "");


  options_->GetNumericValue("penalty_parameter_init_value",
                            penalty_parameter_init_value_, "");
  options_->GetNumericValue("penalty_update_tol", penalty_update_tol_, "");
  options_->GetNumericValue("penalty_parameter_increase_factor",
                            penalty_parameter_increase_factor_, "");
  options_->GetNumericValue("penalty_parameter_max_value",
                            penalty_parameter_max_value_, "");
  options_->GetNumericValue("objective_scaling_factor",
                            objective_scaling_factor_, "");
  options_->GetNumericValue("eps1", eps1_, "");
  options_->GetNumericValue("eps1_change_parm", eps1_change_parm_, "");
  options_->GetNumericValue("eps2", eps2_, "");
  options_->GetIntegerValue("penalty_iter_max", penalty_iter_max_, "");

  options_->GetBoolValue("perform_second_order_correction",
                         perform_second_order_correction_step_, "");
  options_->GetBoolValue("slack_formulation", slack_formulation_, "");

  options_->GetNumericValue("active_set_tol", active_set_tol_, "");
  options_->GetNumericValue("opt_tol", opt_tol_, "");
  options_->GetNumericValue("opt_tol_primal_feasibility",
                            opt_tol_primal_feasibility_, "");
  options_->GetNumericValue("opt_tol_dual_feasibility",
                            opt_tol_dual_feasibility_, "");
  options_->GetNumericValue("opt_tol_stationarity_feasibility",
                            opt_tol_stationarity_feasibility_, "");
  options_->GetNumericValue("opt_tol_complementarity", opt_tol_complementarity_,
                            "");

  /** QP solver usde for ??? */
  int enum_int;
  options_->GetEnumValue("qp_solver", enum_int, "");
  qp_solver_choice_ = QpSolver(enum_int);
}

void SqpSolver::second_order_correction_()
{
  jnlst_->Printf(J_DETAILED, J_MAIN, "\nTry second-order correction step.\n");

  // Store backup of original trial step for the case that the second-order
  // correction step is not accepted.
  Vector trial_step_backup(*trial_step_);
  double trial_model_infeasibility_backup = trial_model_infeasibility_;
  Vector trial_iterate_backup(*trial_iterate_);
  double trial_objective_value_backup = trial_objective_value_;
  Vector trial_constraint_values_backup(*trial_constraint_values_);
  double trial_infeasibility_backup = trial_infeasibility_;
  Vector trial_constraint_multipliers_backup(*trial_constraint_multipliers_);
  Vector trial_bound_multipliers_backup(*trial_bound_multipliers_);
  double predicted_reduction_backup = predicted_reduction_;
  double actual_reduction_backup = actual_reduction_;

  // Compute new constraint values by adding the current step
  shared_ptr<Vector> qp_constraint_body =
      make_shared<Vector>(*current_constraint_values_);
  current_constraint_jacobian_->multiply(trial_step_, qp_constraint_body);

  // Setup the new QP and solve it
  qp_update_tracker_.trigger_bounds_update();
  calculate_search_direction_(qp_constraint_body);

  // We need to trigger a bounds update because the next QP will be solved with
  // the original constraint body
  qp_update_tracker_.trigger_bounds_update();

  // Compute new trial point
  calc_trial_point_and_values_();

  // Check the ratio text again, but using the (not-updated) predicted reduction
  // from the original step
  perform_ratio_test_();

  // TAKE ME OUT!
  trial_point_is_accepted_ = false;
  // If second-order correction step is not accepted, restore the original
  // values
  if (!trial_point_is_accepted_) {
    trial_step_->copy_vector(trial_step_backup);
    trial_model_infeasibility_ = trial_model_infeasibility_backup;
    trial_iterate_->copy_vector(trial_iterate_backup);
    trial_objective_value_ = trial_objective_value_backup;
    trial_constraint_values_->copy_vector(trial_constraint_values_backup);
    trial_infeasibility_ = trial_infeasibility_backup;
    predicted_reduction_ = predicted_reduction_backup;
    actual_reduction_ = actual_reduction_backup;
    trial_constraint_multipliers_->copy_vector(
        trial_constraint_multipliers_backup);
    trial_bound_multipliers_->copy_vector(trial_bound_multipliers_backup);
  }
}

////////////////////////////////////////////////////////////////////////////

void SqpSolver::print_final_stats_()
{
  jnlst_->Printf(J_SUMMARY, J_MAIN, "======================================"
                                    "===================================="
                                    "==================================\n");

  // Determine string describing the exit status
  string exit_message;
  switch (exit_flag_) {
    case OPTIMAL:
      exit_message = "Optimal solution found.";
      break;
    case PRED_REDUCTION_NEGATIVE:
      exit_message = "Error: Predict reduction is negative.";
      break;
    case INVALID_NLP:
      exit_message = "Error: Invalid NLP.";
      break;
    case EXCEED_MAX_ITERATIONS:
      exit_message = "Maximum number of iterations exceeded.";
      break;
    case EXCEED_MAX_CPU_TIME:
      exit_message = "CPU time limit exceeded.";
      break;
    case EXCEED_MAX_WALLCLOCK_TIME: // TODO NEXT
      exit_message = "Wallclock time limit exceeded.";
      break;
    case TRUST_REGION_TOO_SMALL:
      exit_message = "Trust region becomes too small.";
      break;
    case PENALTY_TOO_LARGE:
      exit_message = "Penalty parameter becomes too large.";
      break;
    case QPERROR_INFEASIBLE:
      exit_message = "Error: QP solver claims that QP is infeasible.";
      break;
    case QPERROR_UNBOUNDED:
      exit_message = "Error: QP solver claims that QP is unbounded.";
      break;
    case QPERROR_EXCEED_MAX_ITER:
      exit_message = "Error: QP solver exceeded internal iteration limit.";
      break;
    case QPERROR_UNKNOWN:
      exit_message = "Error: Unknown QP solver error.";
      break;
    default:
      exit_message =
          "Error: exit_flag has uncaught value " + to_string(exit_flag_) + ".";
      break;
  }

  // Determine the number of active constraints (only if solve was successful)
  int num_active_lower_bounds = 0;
  int num_active_upper_bounds = 0;
  int num_active_lower_inequalities = 0;
  int num_active_upper_inequalities = 0;

  if (exit_flag_ == OPTIMAL) {
    const ActivityStatus* bound_activity_status;
    bound_activity_status = qp_solver_->get_bounds_working_set();
    for (int i = 0; i < num_variables_; ++i) {
      if (bound_activity_status[i] == ACTIVE_ABOVE) {
        num_active_upper_bounds++;
      } else if (bound_activity_status[i] == ACTIVE_BELOW) {
        num_active_lower_bounds++;
      }
    }

    const ActivityStatus* constraint_activity_status;
    constraint_activity_status = qp_solver_->get_constraints_working_set();
    for (int i = 0; i < num_constraints_; ++i) {
      if (constraint_type_[i] == IS_EQUALITY) {
        continue;
      }

      if (constraint_activity_status[i] == ACTIVE_ABOVE) {
        num_active_upper_inequalities++;
      } else if (constraint_activity_status[i] == ACTIVE_BELOW) {
        num_active_lower_inequalities++;
      }
    }
  }

  // Print the exit status
  jnlst_->Printf(J_SUMMARY, J_MAIN,
                 "\nExit status............................:  %s\n",
                 exit_message.c_str());
  jnlst_->Printf(J_SUMMARY, J_MAIN,
                 "Number of Variables....................:  %d\n",
                 num_variables_);
  jnlst_->Printf(J_SUMMARY, J_MAIN,
                 "Number of Equality Constraints.........:  %d\n",
                 num_equality_constraints_);
  jnlst_->Printf(J_SUMMARY, J_MAIN,
                 "Number of Inquality Constraints........:  %d\n\n",
                 num_constraints_ - num_equality_constraints_);

  jnlst_->Printf(J_SUMMARY, J_MAIN,
                 "Number of Major Iterations.............:  %d\n",
                 solver_statistics_->num_sqp_iterations_);
  jnlst_->Printf(J_SUMMARY, J_MAIN,
                 "Number of QP Solver Iterations.........:  %d\n\n",
                 solver_statistics_->num_qp_iterations_);

  jnlst_->Printf(J_SUMMARY, J_MAIN,
                 "Final Objectives.......................: %23.16e\n",
                 current_objective_value_);
  jnlst_->Printf(J_SUMMARY, J_MAIN,
                 "Constraint Violation...................: %23.16e\n",
                 current_kkt_error_.primal_infeasibility);
  jnlst_->Printf(J_SUMMARY, J_MAIN,
                 "Dual Infeasibility.....................: %23.16e\n",
                 current_kkt_error_.dual_infeasibility);
  jnlst_->Printf(J_SUMMARY, J_MAIN,
                 "Complmentarity Violation...............: %23.16e\n",
                 current_kkt_error_.complementarity_violation);
  jnlst_->Printf(J_SUMMARY, J_MAIN,
                 "||c_k||................................: %23.16e\n\n",
                 current_infeasibility_);

  if (exit_flag_ == OPTIMAL) {
    jnlst_->Printf(J_SUMMARY, J_MAIN,
                   "Number of active lower bounds..........:  %d\n",
                   num_active_lower_bounds);
    jnlst_->Printf(J_SUMMARY, J_MAIN,
                   "Number of active upper bounds..........:  %d\n",
                   num_active_upper_bounds);
    jnlst_->Printf(J_SUMMARY, J_MAIN,
                   "Number of active lower inequalities....:  %d\n",
                   num_active_lower_inequalities);
    jnlst_->Printf(J_SUMMARY, J_MAIN,
                   "Number of active upper inequalities....:  %d\n\n",
                   num_active_upper_inequalities);
  }

  double cpu_time_now = get_cpu_time_since_start();
  double wallclock_time_now = get_wallclock_time_since_start();

  jnlst_->Printf(J_SUMMARY, J_MAIN,
                 "CPU time used..........................: %12.4f secs\n",
                 cpu_time_now - cpu_time_at_start_);
  jnlst_->Printf(J_SUMMARY, J_MAIN,
                 "Wall clock time passed.................: %12.4f secs\n",
                 wallclock_time_now - wallclock_time_at_start_);

  jnlst_->Printf(J_SUMMARY, J_MAIN, DOUBLE_LONG_DIVIDER);
}

bool SqpSolver::eval_f_(shared_ptr<const Vector> x, double& obj_value)
{
  bool retval = sqp_nlp_->eval_f(x, obj_value);
  if (retval) {
    obj_value *= objective_scaling_factor_;
  }
  return retval;
}

bool SqpSolver::eval_constraints_(shared_ptr<const Vector> x,
                                  shared_ptr<Vector> constraints)
{
  bool retval = sqp_nlp_->eval_constraints(x, constraints);
  return retval;
}

bool SqpSolver::eval_gradient_(shared_ptr<const Vector> x,
                               shared_ptr<Vector> gradient)
{
  bool retval = sqp_nlp_->eval_gradient(x, gradient);
  if (retval && objective_scaling_factor_ != 1.) {
    gradient->scale(objective_scaling_factor_);
  }
  return retval;
}

bool SqpSolver::get_jacobian_structure_(
    shared_ptr<const Vector> x, shared_ptr<SparseTripletMatrix> Jacobian)
{
  bool retval = sqp_nlp_->get_jacobian_structure(x, Jacobian);
  return retval;
}

bool SqpSolver::eval_jacobian_(shared_ptr<const Vector> x,
                               shared_ptr<SparseTripletMatrix> Jacobian)
{
  bool retval = sqp_nlp_->eval_jacobian(x, Jacobian);
  return retval;
}

bool SqpSolver::get_hessian_structure_(shared_ptr<const Vector> x,
                                       shared_ptr<const Vector> lambda,
                                       shared_ptr<SparseTripletMatrix> Hessian)
{
  // We assume here that the values of lambda do not matter and do not scale
  // lambda
  bool retval = sqp_nlp_->get_hessian_structure(x, lambda, Hessian);
  return retval;
}

bool SqpSolver::eval_hessian_(shared_ptr<const Vector> x,
                              shared_ptr<const Vector> lambda,
                              shared_ptr<SparseTripletMatrix> Hessian)
{
  bool retval =
      sqp_nlp_->eval_hessian(x, lambda, objective_scaling_factor_, Hessian);
  return retval;
}

void SqpSolver::store_watchdog_backups_()
{
  backup_current_iterate_ = make_shared<Vector>(*current_iterate_);
  backup_current_constraint_multipliers_ = make_shared<Vector>(*current_constraint_multipliers_);
  backup_current_bound_multipliers_ = make_shared<Vector>(*current_bound_multipliers_);

  backup_current_objective_value_ = current_objective_value_;
  backup_current_objective_gradient_ = make_shared<Vector>(*current_objective_gradient_);
  backup_current_constraint_values_ = make_shared<Vector>(*current_constraint_values_);
  backup_current_constraint_jacobian_ = make_shared<SparseTripletMatrix>(*current_constraint_jacobian_);
  backup_current_lagrangian_hessian_ = make_shared<SparseTripletMatrix>(*current_lagrangian_hessian_);

  backup_current_infeasibility_ = current_infeasibility_;
  backup_predicted_reduction_ = predicted_reduction_;
  backup_current_penalty_parameter_ = current_penalty_parameter_;
  backup_trust_region_radius_ = trust_region_radius_;

  backup_trial_step_ = make_shared<Vector>(*trial_step_);
  backup_trial_constraint_multipliers_ = make_shared<Vector>(*trial_constraint_multipliers_);
  backup_trial_bound_multipliers_ = make_shared<Vector>(*trial_bound_multipliers_);
  backup_trial_model_infeasibility_ = trial_model_infeasibility_;
}

void SqpSolver::delete_watchdog_backups_()
{
  backup_current_iterate_ .reset();
  backup_current_constraint_multipliers_.reset();
  backup_current_bound_multipliers_.reset();

  backup_current_objective_gradient_.reset();
  backup_current_constraint_values_.reset();
  backup_current_constraint_jacobian_.reset();
  backup_current_lagrangian_hessian_.reset();

  backup_trial_step_.reset();
  backup_trial_constraint_multipliers_.reset();
  backup_trial_bound_multipliers_.reset();
}

void SqpSolver::restore_watchdog_backups_()
{
  current_iterate_ = backup_current_iterate_;
  current_constraint_multipliers_ = backup_current_constraint_multipliers_;
  current_bound_multipliers_ = backup_current_bound_multipliers_;

  current_objective_value_ = backup_current_objective_value_;
  current_objective_gradient_ = backup_current_objective_gradient_;
  current_constraint_values_ = backup_current_constraint_values_;
  current_constraint_jacobian_ = backup_current_constraint_jacobian_;
  current_lagrangian_hessian_ = backup_current_lagrangian_hessian_;

  current_infeasibility_ = backup_current_infeasibility_;
  predicted_reduction_ = backup_predicted_reduction_;
  current_penalty_parameter_ = backup_current_penalty_parameter_;
  trust_region_radius_ = backup_trust_region_radius_;

  trial_step_ = backup_trial_step_;
  trial_constraint_multipliers_ = backup_trial_constraint_multipliers_;
  trial_bound_multipliers_ = backup_trial_bound_multipliers_;
  trial_model_infeasibility_ = backup_trial_model_infeasibility_;

  // Tell QP solver that all data has changed
  qp_update_tracker_.trigger_all_updates();
}
} // END_NAMESPACE_SQPHOTSTART
