/* Copyright (C) 2019
* All Rights Reserved.
*
* Authors: Xinyi Luo
* Date:2019-06
*/

#include "sqphot/SqpAlgorithm.hpp"
#include "sqphot/MessageHandling.hpp"
#include "sqphot/SQPDebug.hpp"

#include <fstream>

using namespace std;
using namespace Ipopt;

namespace SQPhotstart {

DECLARE_STD_EXCEPTION(NEW_POINTS_WITH_INCREASE_OBJ_ACCEPTED);
DECLARE_STD_EXCEPTION(SMALL_TRUST_REGION);

/**
 * Default Constructor
 */
SqpAlgorithm::SqpAlgorithm()
 : constraint_activity_status_(nullptr)
 , current_output_file_name_("")
 , bound_activity_status_(nullptr)
 , constraint_type_(nullptr)
 , bound_type_(nullptr)
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
  SmartPtr<RegisteredOptions> reg_options = new RegisteredOptions();

  // Set the algorithm specific options
  register_options_(reg_options);

  // Finalize options list by giving it the list of registered options
  // and the journalist (for error message).
  options_->SetJournalist(jnlst_);
  options_->SetRegisteredOptions(reg_options);
}

/**
 * Destructor
 */
SqpAlgorithm::~SqpAlgorithm()
{
  free_memory_();
}

/**
 *  Free any allocated memory.
 */
void SqpAlgorithm::free_memory_()
{
  delete[] constraint_type_;
  constraint_type_ = nullptr;
  delete[] bound_type_;
  bound_type_ = nullptr;
  delete[] bound_activity_status_;
  bound_activity_status_ = nullptr;
  delete[] constraint_activity_status_;
  constraint_activity_status_ = nullptr;
}

void SqpAlgorithm::initialize_for_new_nlp_(const string& options_file_name)
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

  // Copy all options values into this object
  get_option_values_();

  // Free any existing objects or memory
  free_memory_();

  // Allocate all objects and memory needed in the loop.  This also creates the
  // QP and LP solver objects.
  allocate_memory_();
}

void SqpAlgorithm::print_initial_output_()
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
          assert(false && "should not get here");
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

  // Print header of summary output table and the first line for the initial
  // iterate
  jnlst_->Printf(J_ITERSUMMARY, J_MAIN, "======================================"
                                        "======================================"
                                        "==========================\n");
  jnlst_->Printf(J_ITERSUMMARY, J_MAIN,
                 "%6s  %23s    %9s    %9s    %9s    %9s    %9s  \n", "iter",
                 "f", "||p_k||", "||c_k||", "Delta", "rho", "QP_KKT_Error");
  jnlst_->Printf(J_ITERSUMMARY, J_MAIN, "======================================"
                                        "===================================="
                                        "============================\n");
  jnlst_->Printf(J_ITERSUMMARY, J_MAIN,
                 "%6i  %23.16e  %9.6e  %9.6e  %9.6e  %9.6e  %9.6e\n",
                 solver_statistics_->num_sqp_iterations_, obj_value_,
                 trial_step_norm_, current_infeasibility_, trust_region_radius_,
                 penalty_parameter_,
                 qp_solver_->get_QpOptimalStatus().KKT_error);
}

void SqpAlgorithm::calculate_search_direction_()
{
  // Set the data for the QP corresponding to the current iterate
  setup_qp_();
  // for debugging
  //@{
  //    hessian_->print_full("hessian");
  //    jacobian_->print_full("jacobian");
  //@}

  // Solve the QP
  QpSolverExitStatus exit_status = qp_solver_->solve_qp(solver_statistics_);
  if (exit_status != QPEXIT_OPTIMAL) {
    const string& nlp_name = sqp_nlp_->get_nlp_name();
    qp_solver_->write_qp_data(nlp_name + "qpdata.log");
    assert(false && "Still need to decide how to handle QP solver error.");
    // exit_flag_ = qp_solver_->get_status();
    // break;
  }

  // Extract the trial step from the QP solution.
  extract_trial_step_();
}

/**
 * @brief This is the main function to optimize the NLP given as the input
 *
 * @param nlp: the nlp reader that read data of the function to be minimized;
 */
void SqpAlgorithm::optimize_nlp(shared_ptr<SqpNlpBase> sqp_nlp,
                                std::string options_file_name)
{
  // Make NLP object directly accessible to all methods
  sqp_nlp_ = sqp_nlp;

  // Store the

  // Set up everything needed to execute the loop: Read options, allocate
  // memory, create QP solver objects
  initialize_for_new_nlp_(options_file_name);

  // Set the time at beginning of the algorithm
  cpu_time_at_start_ = get_cpu_time_since_start();
  wallclock_time_at_start_ = get_wallclock_time_since_start();

  // Initializes the iterate and all related quantities
  initialize_iterates_();

  // Print initial output
  print_initial_output_();

  // Initialize exit flag to UNKOWN to indicate that loop is not finished
  exit_flag_ = UNKNOWN;
  while (solver_statistics_->num_sqp_iterations_ < max_num_iterations_ &&
         exit_flag_ == UNKNOWN) {

    // Solve the QP to get the search direction
    calculate_search_direction_();

    //        p_k_->print("p_k_");
    qp_obj_ = get_obj_QP();

    // Update the penalty parameter if necessary
    update_penalty_parameter_();

    // calculate the infinity norm of the search direction
    trial_step_norm_ = trial_step_->calc_inf_norm();

    get_trial_point_info();

    ratio_test();

    // Calculate the second-order-correction steps
    second_order_correction();

    // Update the radius and the QP bounds if the radius has been changed
    solver_statistics_->increase_sqp_iteration_counter();
    /* output some information to the console*/

    if (solver_statistics_->num_sqp_iterations_ % 10 == 0) {
      jnlst_->Printf(J_ITERSUMMARY, J_MAIN,
                     "%6s  %23s    %9s    %9s    %9s    %9s    %9s  \n", "iter",
                     "f", "||p_k||", "||c_k||", "Delta", "rho", "QP_KKT_Error");
      jnlst_->Printf(J_ITERSUMMARY, J_MAIN, "=================================="
                                            "=================================="
                                            "======"
                                            "============================\n");
    }
    jnlst_->Printf(J_ITERSUMMARY, J_MAIN,
                   "%6i  %23.16e  %9.6e  %9.6e  %9.6e  %9.6e  %9.6e\n",
                   solver_statistics_->num_sqp_iterations_, obj_value_,
                   trial_step_norm_, current_infeasibility_,
                   trust_region_radius_, penalty_parameter_,
                   qp_solver_->get_QpOptimalStatus().KKT_error);

    // check if the current iterates is optimal and decide to
    // exit the loop or not
    check_optimality();
    if (exit_flag_ != UNKNOWN) {
      break;
    }

    try {
      update_radius();
    } catch (SMALL_TRUST_REGION) {
      check_optimality();
      break;
    }
  }

  // check if the current iterates get_status before exiting
  if (solver_statistics_->num_sqp_iterations_ == max_num_iterations_)
    exit_flag_ = EXCEED_MAX_ITERATIONS;

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
  print_final_stats();
  jnlst_->FlushBuffer();
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
void SqpAlgorithm::check_optimality()
{
  // FIXME: not sure if it is better to use the new multiplier or the old one

  double primal_violation = 0;
  double dual_violation = 0;
  double compl_violation = 0;
  double stationarity_violation = 0;

  int i;
  get_multipliers();

  if (bound_multipliers_ == nullptr)
    /**-------------------------------------------------------**/
    /**               Check the KKT conditions                **/
    /**-------------------------------------------------------**/

    /**-------------------------------------------------------**/
    /**                    Identify Active Set                **/
    /**-------------------------------------------------------**/
    for (i = 0; i < num_constraints_; i++) {
      if (constraint_type_[i] == BOUNDED_ABOVE) {
        if (abs(upper_constraint_bounds_->get_value(i) -
                current_constraint_values_->get_value(i)) < active_set_tol_)
          constraint_activity_status_[i] = ACTIVE_ABOVE;
      } else if (constraint_type_[i] == BOUNDED_BELOW) {
        if (abs(current_constraint_values_->get_value(i) -
                lower_constraint_bounds_->get_value(i)) < active_set_tol_) {
          constraint_activity_status_[i] = ACTIVE_BELOW;
        }
      } else if (constraint_type_[i] == IS_EQUALITY) {
        if ((abs(upper_constraint_bounds_->get_value(i) -
                 current_constraint_values_->get_value(i)) < active_set_tol_) &&
            (abs(current_constraint_values_->get_value(i) -
                 lower_constraint_bounds_->get_value(i)) < active_set_tol_))
          constraint_activity_status_[i] = ACTIVE_BOTH_SIDES;
      } else {
        constraint_activity_status_[i] = INACTIVE;
      }
    }

  for (i = 0; i < num_variables_; i++) {
    if (bound_type_[i] == BOUNDED_ABOVE) {
      if (abs(upper_variable_bounds_->get_value(i) -
              current_iterate_->get_value(i)) < active_set_tol_)
        bound_activity_status_[i] = ACTIVE_ABOVE;
    } else if (bound_type_[i] == BOUNDED_BELOW) {
      if (abs(current_iterate_->get_value(i) -
              lower_variable_bounds_->get_value(i)) < active_set_tol_)
        bound_activity_status_[i] = ACTIVE_BELOW;
    } else if (bound_type_[i] == IS_EQUALITY) {
      if ((abs(upper_variable_bounds_->get_value(i) -
               current_iterate_->get_value(i)) < active_set_tol_) &&
          (abs(current_iterate_->get_value(i) -
               lower_variable_bounds_->get_value(i)) < active_set_tol_))
        bound_activity_status_[i] = ACTIVE_BOTH_SIDES;
    } else {
      bound_activity_status_[i] = INACTIVE;
    }
  }
  /**-------------------------------------------------------**/
  /**                    Primal Feasibility                 **/
  /**-------------------------------------------------------**/

  primal_violation += current_infeasibility_;

  /**-------------------------------------------------------**/
  /**                    Dual Feasibility                   **/
  /**-------------------------------------------------------**/

  //    multiplier_vars_->print("multiplier_vars");
  //    multiplier_cons_->print("multiplier_cons");
  //    x_k_->print("x_k");
  //    lower_variable_bounds__->print("lower_variable_bounds_");
  //    upper_variable_bounds__->print("upper_variable_bounds_");
  //
  //    c_k_->print("c_k");
  //    lower_constraint_bounds__->print("lower_constraint_bounds_");
  //    upper_constraint_bounds__->print("upper_constraint_bounds_");
  i = 0;
  opt_status_.dual_feasibility = true;
  while (i < num_variables_) {
    if (bound_type_[i] == BOUNDED_ABOVE) {
      dual_violation += max(bound_multipliers_->get_value(i), 0.0);
    } else if (bound_type_[i] == BOUNDED_BELOW) {
      dual_violation += -min(bound_multipliers_->get_value(i), 0.0);
    }
    i++;
  }

  i = 0;
  while (i < num_constraints_) {
    if (constraint_type_[i] == BOUNDED_ABOVE) {
      dual_violation += max(constraint_multipliers_->get_value(i), 0.0);
    } else if (constraint_type_[i] == BOUNDED_BELOW) {
      dual_violation += -min(constraint_multipliers_->get_value(i), 0.0);
    }
    i++;
  }

  /**-------------------------------------------------------**/
  /**                    Complemtarity                      **/
  /**-------------------------------------------------------**/
  //@{

  i = 0;
  while (i < num_constraints_) {
    if (constraint_type_[i] == BOUNDED_ABOVE) {
      compl_violation += abs(constraint_multipliers_->get_value(i) *
                             (upper_constraint_bounds_->get_value(i) -
                              current_constraint_values_->get_value(i)));
    } else if (constraint_type_[i] == BOUNDED_BELOW) {
      compl_violation += abs(constraint_multipliers_->get_value(i) *
                             (current_constraint_values_->get_value(i) -
                              lower_constraint_bounds_->get_value(i)));
    } else if (constraint_type_[i] == UNBOUNDED) {
      compl_violation += abs(constraint_multipliers_->get_value(i));
    }
    i++;
  }

  i = 0;
  while (i < num_variables_) {
    if (bound_type_[i] == BOUNDED_ABOVE) {
      compl_violation += abs(bound_multipliers_->get_value(i) *
                             (upper_variable_bounds_->get_value(i) -
                              current_iterate_->get_value(i)));
    } else if (bound_type_[i] == BOUNDED_BELOW) {
      compl_violation += abs(bound_multipliers_->get_value(i) *
                             (current_iterate_->get_value(i) -
                              lower_variable_bounds_->get_value(i)));
    } else if (bound_type_[i] == UNBOUNDED) {
      compl_violation += abs(bound_multipliers_->get_value(i));
    }
    i++;
  }
  //@{
  //    lower_constraint_bounds__->print("lower_constraint_bounds_");
  //    upper_constraint_bounds__->print("upper_constraint_bounds_");
  //    cout<<compl_violation<<endl;
  //    multiplier_vars_->print("multiplier_vars");
  //    multiplier_cons_->print("multiplier_constr");
  //@}

  //@}
  /**-------------------------------------------------------**/
  /**                    Stationarity                       **/
  /**-------------------------------------------------------**/
  //@{
  shared_ptr<Vector> difference = make_shared<Vector>(num_variables_);
  // the difference of g-J^T y -\lambda
  //
  //    grad_f_->print("grad_f_");
  //
  //    jacobian_->print_full("jacobian");
  //
  //    multiplier_cons_->print("multiplier_cons_");
  //    multiplier_vars_->print("multiplier_vars_");
  current_constraint_jacobian_->multiply_transpose(constraint_multipliers_,
                                                   difference);
  difference->add_vector(1., bound_multipliers_);
  difference->add_vector(-1., current_objective_gradient_);

  stationarity_violation = difference->calc_one_norm();
  //@}

  /**-------------------------------------------------------**/
  /**             Decide if x_k is optimal                  **/
  /**-------------------------------------------------------**/

  opt_status_.dual_violation = dual_violation;
  opt_status_.primal_violation = primal_violation;
  opt_status_.compl_violation = compl_violation;
  opt_status_.stationarity_violation = stationarity_violation;
  opt_status_.KKT_error = dual_violation + primal_violation + compl_violation +
                          stationarity_violation;
  //    printf("primal_violation = %23.16e\n",primal_violation);
  //    printf("dual_violation = %23.16e\n",dual_violation);
  //    printf("compl_violation = %23.16e\n",compl_violation);
  //    printf("statioanrity_violation = %23.16e\n",statioanrity_violation);
  //    printf("KKT error = %23.16e\n", opt_status_.KKT_error);

  opt_status_.primal_feasibility =
      primal_violation < opt_tol_primal_feasibility_;
  opt_status_.dual_feasibility = dual_violation < opt_tol_dual_feasibility_;
  opt_status_.complementarity = compl_violation < opt_tol_complementarity_;
  opt_status_.stationarity =
      stationarity_violation < opt_tol_stationarity_feasibility_;

  if (opt_status_.primal_feasibility && opt_status_.dual_feasibility &&
      opt_status_.complementarity && opt_status_.stationarity) {
    opt_status_.first_order_opt = true;
    exit_flag_ = OPTIMAL;
  } else {
// if it is not optimal
#ifdef DEBUG
#ifdef CHECK_TERMINATION

    EJournalLevel debug_print_level = old_options_->debug_print_level;
    SmartPtr<Journal> debug_jrnl = jnlst_->GetJournal("Debug");
    if (IsNull(debug_jrnl)) {
      debug_jrnl = jnlst_->AddFileJournal("Debug", "debug.out", J_ITERSUMMARY);
    }
    debug_jrnl->SetAllPrintLevels(debug_print_level);
    debug_jrnl->SetPrintLevel(J_DBG, J_ALL);
    jnlst_->Printf(J_ALL, J_DBG, DOUBLE_DIVIDER);
    jnlst_->Printf(J_ALL, J_DBG, "           Iteration  %i\n", stats_->iter);
    jnlst_->Printf(J_ALL, J_DBG, DOUBLE_DIVIDER);
    grad_f_->print("grad_f", jnlst_, J_MOREDETAILED, J_DBG);
    jnlst_->Printf(J_ALL, J_DBG, SINGLE_DIVIDER);
    upper_constraint_bounds__->print("upper_constraint_bounds_", jnlst_,
                                     J_MOREDETAILED, J_DBG);
    lower_constraint_bounds__->print("lower_constraint_bounds_", jnlst_,
                                     J_MOREDETAILED, J_DBG);
    c_k_->print("c_k", jnlst_, J_MOREDETAILED, J_DBG);
    multiplier_cons_->print("multiplier_cons", jnlst_, J_MOREDETAILED, J_DBG);
    jnlst_->Printf(J_ALL, J_DBG, SINGLE_DIVIDER);
    upper_variable_bounds__->print("upper_variable_bounds_", jnlst_,
                                   J_MOREDETAILED, J_DBG);
    lower_variable_bounds__->print("lower_variable_bounds_", jnlst_,
                                   J_MOREDETAILED, J_DBG);
    x_k_->print("x_k", jnlst_, J_MOREDETAILED, J_DBG);
    multiplier_vars_->print("multiplier_vars", jnlst_, J_MOREDETAILED, J_DBG);
    jnlst_->Printf(J_ALL, J_DBG, SINGLE_DIVIDER);
    jacobian_->print_full("jacobian", jnlst_, J_MOREDETAILED, J_DBG);
    hessian_->print_full("hessian", jnlst_, J_MOREDETAILED, J_DBG);
    difference->print("stationarity gap", jnlst_, J_MOREDETAILED, J_DBG);
    jnlst_->Printf(J_ALL, J_DBG, SINGLE_DIVIDER);
    jnlst_->Printf(J_ALL, J_DBG, "Feasibility      ");
    jnlst_->Printf(J_ALL, J_DBG, "%i\n", opt_status_.primal_feasibility);
    jnlst_->Printf(J_ALL, J_DBG, "Dual Feasibility ");
    jnlst_->Printf(J_ALL, J_DBG, "%i\n", opt_status_.dual_feasibility);
    jnlst_->Printf(J_ALL, J_DBG, "Stationarity     ");
    jnlst_->Printf(J_ALL, J_DBG, "%i\n", opt_status_.stationarity);
    jnlst_->Printf(J_ALL, J_DBG, "Complementarity  ");
    jnlst_->Printf(J_ALL, J_DBG, "%i\n", opt_status_.complementarity);
    jnlst_->Printf(J_ALL, J_DBG, SINGLE_DIVIDER);

#endif
#endif
  }
}

void SqpAlgorithm::get_trial_point_info()
{

  trial_iterate_->set_to_sum_of_vectors(1., current_iterate_, 1., trial_step_);

  // Calculate f_trial, c_trial and current_infeasibility__trial for the trial
  // points
  // x_trial
  sqp_nlp_->eval_f(trial_iterate_, obj_value_trial_);
  sqp_nlp_->eval_constraints(trial_iterate_, trial_constraint_values_);

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

void SqpAlgorithm::initialize_iterates_()
{
  // Determine the bound values
  sqp_nlp_->get_bounds_info(lower_variable_bounds_, upper_variable_bounds_,
                            lower_constraint_bounds_, upper_constraint_bounds_);

  // Determine inequality types
  classify_constraints_types_();

  // Determine the starting point
  sqp_nlp_->get_starting_point(current_iterate_, constraint_multipliers_);

// shift starting point to satisfy the bound constraint
#ifndef NEW_FORMULATION
  shift_starting_point_(current_iterate_, lower_variable_bounds_,
                        upper_variable_bounds_);
#endif

  // Compute function and derivative values at the starting point
  sqp_nlp_->eval_f(current_iterate_, obj_value_);
  sqp_nlp_->eval_gradient(current_iterate_, current_objective_gradient_);
  sqp_nlp_->eval_constraints(current_iterate_, current_constraint_values_);
  sqp_nlp_->get_hessian_structure(current_iterate_, constraint_multipliers_,
                                  current_lagrangian_hessian_);
  sqp_nlp_->eval_hessian(current_iterate_, constraint_multipliers_,
                         current_lagrangian_hessian_);
  sqp_nlp_->get_jacobian_structure(current_iterate_,
                                   current_constraint_jacobian_);
  sqp_nlp_->eval_jacobian(current_iterate_, current_constraint_jacobian_);

  // Initalize algorithmic quantities
  trust_region_radius_ = trust_region_init_value_;
  penalty_parameter_ = penalty_parameter_init_value_;
  trial_step_norm_ = 0.0;

#ifdef NEW_FORMULATION
  current_infeasibility_ = compute_constraint_violation_(
      current_iterate_, current_constraint_values_);
#else
  current_infeasibility_ = compute_constraint_violation_(
      current_iterate_, current_constraint_values_);
#endif

  // Set flag that makes sure that the QP solver will be initialized at the
  // first call
  qp_solver_initialized_ = false;
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
void SqpAlgorithm::allocate_memory_()
{
  // Determine the problem size
  shared_ptr<const SqpNlpSizeInfo> nlp_sizes = sqp_nlp_->get_problem_sizes();
  num_variables_ = nlp_sizes->get_num_variables();
  num_constraints_ = nlp_sizes->get_num_constraints();

  // Allocate memory for the arrays that hold the constraint types
  constraint_type_ = new ConstraintType[num_constraints_];
  bound_type_ = new ConstraintType[num_variables_];

  // Allocate memory for the arrays that hold the activity status
  bound_activity_status_ = new ActivityStatus[num_variables_];
  constraint_activity_status_ = new ActivityStatus[num_constraints_];

  // Create Vector objects that hold the bounds
  lower_variable_bounds_ = make_shared<Vector>(num_variables_);
  upper_variable_bounds_ = make_shared<Vector>(num_variables_);
  lower_constraint_bounds_ = make_shared<Vector>(num_constraints_);
  upper_constraint_bounds_ = make_shared<Vector>(num_constraints_);

  // Create Vector objects that hold the iterate values
  current_iterate_ = make_shared<Vector>(num_variables_);
  trial_iterate_ = make_shared<Vector>(num_variables_);
  constraint_multipliers_ = make_shared<Vector>(num_constraints_);
  bound_multipliers_ = make_shared<Vector>(num_variables_);
  trial_step_ = make_shared<Vector>(num_variables_);

  // Create objects that hold evaluated quantities
  current_constraint_values_ = make_shared<Vector>(num_constraints_);
  trial_constraint_values_ = make_shared<Vector>(num_constraints_);
  current_objective_gradient_ = make_shared<Vector>(num_variables_);
  current_constraint_jacobian_ =
      make_shared<SpTripletMat>(nlp_sizes->get_num_nonzeros_jacobian(),
                                num_constraints_, num_variables_, false);
  current_lagrangian_hessian_ =
      make_shared<SpTripletMat>(nlp_sizes->get_num_nonzeros_hessian(),
                                num_variables_, num_variables_, true);

  // Create objects that holds the solver statistics
  solver_statistics_ = make_shared<Statistics>();

  // Create the QP solver objects.  This also reads the option values
  qp_solver_ = make_shared<QPhandler>(nlp_sizes, QP, jnlst_, options_);
  lp_solver_ = make_shared<QPhandler>(nlp_sizes, LP, jnlst_, options_);
}

/**
* @brief This function calculates the infeasibility for given x_k and c_k with
respect
* to their corresponding bounds
* @return current_infeasibility_ = ||-max(c_k-c_u),0||_1 +||-min(c_k-c_l),0||_1+
                  ||-max(x_k-upper_variable_bounds_),0||_1
+||-min(x_k-lower_variable_bounds_),0||_1
*/
double SqpAlgorithm::compute_constraint_violation_(
    shared_ptr<const Vector> variable_values,
    shared_ptr<const Vector> constraint_values)
{
  double current_infeasibility_ = 0.0;
  for (int i = 0; i < current_constraint_values_->get_dim(); i++) {
    if (constraint_values->get_value(i) <
        lower_constraint_bounds_->get_value(i))
      current_infeasibility_ += (lower_constraint_bounds_->get_value(i) -
                                 constraint_values->get_value(i));
    else if (constraint_values->get_value(i) >
             upper_constraint_bounds_->get_value(i))
      current_infeasibility_ += (constraint_values->get_value(i) -
                                 upper_constraint_bounds_->get_value(i));
  }

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

/**
 * @brief This function extracts the search direction for NLP from the QP
 * subproblem
 * solved, and copies it to the class member _p_k
 */
void SqpAlgorithm::extract_trial_step_()
{
  // Here, the trial_step_ vector has the length num_variables_, but the primal
  // solution of the QP solver includes the slack variables.
  trial_step_->copy_values(qp_solver_->get_primal_solution()->get_values());
}

/**
 * @brief This function gets the Lagragian multipliers for constraints and for
 * for the
 * bounds from QPsolver
 */

void SqpAlgorithm::get_multipliers()
{
  if (qp_solver_choice_ == QORE || qp_solver_choice_ == QPOASES) {
    constraint_multipliers_->copy_vector(
        qp_solver_->get_constraints_multipliers());
    bound_multipliers_->copy_values(qp_solver_->get_bounds_multipliers()
                                        ->get_values()); // AW would be good to
    // separate different parts
    // of the QP multipliers?
  } else if (qp_solver_choice_ == GUROBI || qp_solver_choice_ == CPLEX) {
    constraint_multipliers_->copy_vector(
        qp_solver_->get_constraints_multipliers());
    shared_ptr<Vector> tmp_vec_nVar = make_shared<Vector>(num_variables_);
    current_constraint_jacobian_->multiply_transpose(constraint_multipliers_,
                                                     tmp_vec_nVar);
    current_lagrangian_hessian_->multiply(trial_step_, bound_multipliers_);
    bound_multipliers_->add_vector(1., current_objective_gradient_);
    bound_multipliers_->add_vector(-1., tmp_vec_nVar);
  }
}

/**
 * @brief This function will set up the data for the QP subproblem
 *
 * It will initialize all the data at once at the beginning of the Algorithm.
 * After
 * that, the data in the QP problem will be updated according to the class
 * member QPinfoFlag_
 */

void SqpAlgorithm::setup_qp_()
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
                           current_constraint_values_);
    qp_solver_->set_gradient(current_objective_gradient_, penalty_parameter_);

    // Initialize the update tracker
    qp_update_tracker_.reset();

    // Nothing more needs to be done.
    qp_solver_initialized_ = true;
    return;
  }

  // Debug check: Make sure that an update is necessary.
  assert(qp_update_tracker_.need_update() && "QP is not changed.");

  // Update those quantities that need to be updated
  if (qp_update_tracker_.need_jacobian_update()) {
    qp_solver_->update_A(current_constraint_jacobian_);
  }
  if (qp_update_tracker_.need_hessian_update()) {
    qp_solver_->update_H(current_lagrangian_hessian_);
  }
  if (qp_update_tracker_.need_bounds_update()) {
    qp_solver_->update_bounds(
        trust_region_radius_, lower_variable_bounds_, upper_variable_bounds_,
        current_iterate_, lower_constraint_bounds_, upper_constraint_bounds_,
        current_constraint_values_);
  } else if (qp_update_tracker_.need_trust_region_radius_update()) {
    qp_solver_->update_delta(trust_region_radius_, lower_variable_bounds_,
                             upper_variable_bounds_, current_iterate_);
  }
  if (qp_update_tracker_.need_penalty_parameter_update()) {
    qp_solver_->update_penalty(penalty_parameter_);
  }
  if (qp_update_tracker_.need_gradient_update()) {
    qp_solver_->update_grad(current_objective_gradient_);
  }

  // Reset the update tracker
  qp_update_tracker_.reset();
}

void SqpAlgorithm::setupLP()
{
  lp_solver_->set_bounds(trust_region_radius_, lower_variable_bounds_,
                         upper_variable_bounds_, current_iterate_,
                         lower_constraint_bounds_, upper_constraint_bounds_,
                         current_constraint_values_);
  lp_solver_->set_gradient(penalty_parameter_);
  lp_solver_->set_jacobian(current_constraint_jacobian_);
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
void SqpAlgorithm::ratio_test()
{

  double P1_x = obj_value_ + penalty_parameter_ * current_infeasibility_;
  double P1_x_trial =
      obj_value_trial_ + penalty_parameter_ * trial_infeasibility_;

  actual_reduction_ = P1_x - P1_x_trial;
  pred_reduction_ = penalty_parameter_ * current_infeasibility_ - get_obj_QP();

#ifdef DEBUG
#ifdef CHECK_TR_ALG
  SmartPtr<Journal> debug_jrnl = jnlst_->GetJournal("Debug");
  if (IsNull(debug_jrnl)) {
    debug_jrnl = jnlst_->AddFileJournal("Debug", "debug.out", J_ITERSUMMARY);
  }
  int debug_print_level = 0;
  debug_jrnl->SetAllPrintLevels((EJournalLevel)debug_print_level);
  debug_jrnl->SetPrintLevel(J_DBG, J_ALL);

  jnlst_->Printf(J_ALL, J_DBG, "\n");
  jnlst_->Printf(J_ALL, J_DBG, SINGLE_DIVIDER);
  jnlst_->Printf(J_ALL, J_DBG, "       The actual reduction is %23.16e\n",
                 actual_reduction_);
  jnlst_->Printf(J_ALL, J_DBG, "       The pred reduction is   %23.16e\n",
                 pred_reduction_);
  double ratio = actual_reduction_ / pred_reduction_;
  jnlst_->Printf(J_ALL, J_DBG, "       The calculated ratio is %23.16e\n",
                 ratio);
  jnlst_->Printf(J_ALL, J_DBG, "       The correct decision is ");
  if (ratio >= trust_region_ratio_accept_tol_)
    jnlst_->Printf(J_ALL, J_DBG, "to ACCEPT the trial point\n");
  else
    jnlst_->Printf(J_ALL, J_DBG, "to REJECT the trial point and change the "
                                 "trust-region radius\n");
  jnlst_->Printf(J_ALL, J_DBG, "       The TRUE decision is ");

  if (actual_reduction_ >= (trust_region_ratio_accept_tol_ * pred_reduction_)) {
    jnlst_->Printf(J_ALL, J_DBG, "to ACCEPT the trial point\n");
  } else
    jnlst_->Printf(J_ALL, J_DBG, "to REJECT the trial point and change the "
                                 "trust-region radius\n");
  jnlst_->Printf(J_ALL, J_DBG, SINGLE_DIVIDER);
  jnlst_->Printf(J_ALL, J_DBG, "\n");

#endif
#endif

#if 1
  if (actual_reduction_ >= (trust_region_ratio_accept_tol_ * pred_reduction_) &&
      actual_reduction_ >= -opt_tol_)
#else
  if (pred_reduction_ < -1.0e-8) {
    //    myQP_->WriteQPData(problem_name_+"qpdata.log");
    exitflag_ = PRED_REDUCTION_NEGATIVE;
    return;
  }
  if (actual_reduction_ >= (trust_region_ratio_accept_tol_ * pred_reduction_))
#endif
  {
    // succesfully update
    // copy information already calculated from the trial point
    current_infeasibility_ = trial_infeasibility_;

    obj_value_ = obj_value_trial_;
    current_iterate_->copy_vector(trial_iterate_);
    current_constraint_values_->copy_vector(trial_constraint_values_);
    // update function information by reading from sqp_nlp_ object
    get_multipliers();
    sqp_nlp_->eval_gradient(current_iterate_, current_objective_gradient_);
    sqp_nlp_->eval_jacobian(current_iterate_, current_constraint_jacobian_);
    sqp_nlp_->eval_hessian(current_iterate_, constraint_multipliers_,
                           current_lagrangian_hessian_);

    qp_update_tracker_.trigger_all_updates();

    isaccept_ = true; // no need to calculate the SOC direction
  } else {
    isaccept_ = false;
  }
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

void SqpAlgorithm::update_radius()
{
  if (actual_reduction_ < trust_region_ratio_decrease_tol_ * pred_reduction_) {
    trust_region_radius_ = trust_region_decrease_factor_ * trust_region_radius_;
    qp_update_tracker_.trigger_trust_region_radius_update();
    // decrease the trust region radius. gamma_c is the parameter in options_
    // object
  } else {
    // printf("delta_ = %23.16e, ||p_k|| = %23.16e\n",delta_,p_k_->inf_norm());
    if (actual_reduction_ >
            trust_region_ratio_increase_tol_ * pred_reduction_ &&
        (opt_tol_ >
         fabs(trust_region_radius_ - trial_step_->calc_inf_norm()))) {
      trust_region_radius_ =
          min(trust_region_increase_factor_ * trust_region_radius_,
              trust_region_max_value_);
      qp_update_tracker_.trigger_trust_region_radius_update();
    }
  }

  // if the trust-region becomes too small, throw the error message

  if (trust_region_radius_ < trust_region_min_value_) {

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
  if (lower_bound > -INF && upper_bound < INF) {
    if ((upper_bound - lower_bound) < 1.0e-8) {
      return IS_EQUALITY;
    } else {
      return BOUNDED_BELOW_AND_ABOVE;
    }
  } else if (lower_bound > -INF && upper_bound > INF) {
    return BOUNDED_BELOW;
  } else if (upper_bound < INF && lower_bound < -INF) {
    return BOUNDED_ABOVE;
  } else {
    return UNBOUNDED;
  }
}

void SqpAlgorithm::classify_constraints_types_()
{

  for (int i = 0; i < num_constraints_; i++) {
    constraint_type_[i] =
        classify_single_constraint(lower_constraint_bounds_->get_value(i),
                                   upper_constraint_bounds_->get_value(i));
  }
  for (int i = 0; i < num_variables_; i++) {
    bound_type_[i] =
        classify_single_constraint(lower_variable_bounds_->get_value(i),
                                   upper_variable_bounds_->get_value(i));
  }
}

/**
 * @brief update the penalty parameter for the algorithm.
 *
 */
void SqpAlgorithm::update_penalty_parameter_()
{

  current_model_infeasibility_ = qp_solver_->get_model_infeasibility();

  // prin/tf("current_infeasibility__model =
  // %23.16e\n",current_infeasibility__model_);
  // printf("current_infeasibility__ = %23.16e\n",current_infeasibility__);
  if (current_model_infeasibility_ > penalty_update_tol_) {
    double current_model_infeasibility_tmp =
        current_model_infeasibility_;      // temporarily store the value
    double rho_trial = penalty_parameter_; // the temporary trial value for rho
    setupLP();

    QpSolverExitStatus exit_status = lp_solver_->solve_lp(solver_statistics_);
    if (exit_status != QPEXIT_OPTIMAL) {
      const string& nlp_name = sqp_nlp_->get_nlp_name();
      qp_solver_->write_qp_data(nlp_name + "qpdata.log");
      assert(false && "Still need to decide how to handle QP solver error.");
      // exit_flag_ = qp_solver_->get_status();
      // break;
    }

    // calculate the current_infeasibility_ of the LP
    double lp_model_infeasibility = lp_solver_->get_model_infeasibility();

    //     printf("lp_model_infeasibility = %23.16e\n",lp_model_infeasibility);
    if (lp_model_infeasibility <= penalty_update_tol_) {
      // try to increase the penalty parameter to a number such that the
      // infeasibility measure of QP model with such penalty parameter
      // becomes zero

      while (current_model_infeasibility_ > penalty_update_tol_) {
        if (rho_trial >= penalty_parameter_max_value_) {
          break;
        } // TODO:safeguarded procedure...put here for now

        rho_trial =
            min(penalty_parameter_max_value_,
                rho_trial * penalty_parameter_increase_factor_); // increase rho

        solver_statistics_->try_new_penalty_parameter();

        qp_solver_->update_penalty(rho_trial);

        QpSolverExitStatus exit_status =
            qp_solver_->solve_qp(solver_statistics_);

        if (exit_status != QPEXIT_OPTIMAL) {
          assert(false &&
                 "Still need to decide how to handle QP solver error.");
          // exit_flag_ = qp_solver_->get_status();
          // break;
        }

        // recalculate the infeasibility measure of the model by
        // calculating the one norm of the slack variables

        current_model_infeasibility_ = qp_solver_->get_model_infeasibility();
      }
    } else {
      while ((current_infeasibility_ - current_model_infeasibility_ <
                  eps1_ * (current_infeasibility_ - lp_model_infeasibility) &&
              (solver_statistics_->num_trial_penalty_parameters_ <
               penalty_iter_max_))) {

        if (rho_trial >= penalty_parameter_max_value_) {
          break;
        }

        // try to increase the penalty parameter to a number such that
        // the incurred reduction for the QP model is to a ratio to the
        // maximum possible reduction for current linear model.
        rho_trial = min(penalty_parameter_max_value_,
                        rho_trial * penalty_parameter_increase_factor_);

        solver_statistics_->try_new_penalty_parameter();

        qp_solver_->update_penalty(rho_trial);

        QpSolverExitStatus exit_status =
            qp_solver_->solve_qp(solver_statistics_);
        if (exit_status != QPEXIT_OPTIMAL) {
          assert(false &&
                 "Still need to decide how to handle QP solver error.");
          // exit_flag_ = qp_solver_->get_status();
          // break;
        }
        // recalculate the infeasibility measure of the model by
        // calculating the one norm of the slack variables

        current_model_infeasibility_ = qp_solver_->get_model_infeasibility();
      }
    }
    // if any change occurs
    if (rho_trial > penalty_parameter_) {
      if (rho_trial * current_infeasibility_ - get_obj_QP() >=
          eps2_ * rho_trial *
              (current_infeasibility_ - current_model_infeasibility_)) {
        //                    printf("rho_trial = %23.16e\n",rho_trial);
        //                    printf("current_infeasibility__ =
        //                    %23.16e\n",current_infeasibility__);
        solver_statistics_->penalty_parameter_increased();

        eps1_ += (1 - eps1_) * eps1_change_parm_;

        // use the new solution as the search direction
        extract_trial_step_();
        penalty_parameter_ = rho_trial; // update to the class variable
        get_trial_point_info();
        qp_obj_ = get_obj_QP();
        double P1_x = obj_value_ + penalty_parameter_ * current_infeasibility_;
        double P1_x_trial =
            obj_value_trial_ + penalty_parameter_ * trial_infeasibility_;

        actual_reduction_ = P1_x - P1_x_trial;
        pred_reduction_ =
            penalty_parameter_ * current_infeasibility_ - get_obj_QP();

      } else {
        //                    printf("rho_trial = %23.16e\n",rho_trial);
        //                    printf("current_infeasibility__ =
        //                    %23.16e\n",current_infeasibility__);
        solver_statistics_->trial_penalty_parameter_not_accepted();
        current_model_infeasibility_ = current_model_infeasibility_tmp;
        qp_update_tracker_.trigger_penalty_parameter_update();
      }
    } else if (current_infeasibility_ < opt_tol_primal_feasibility_ * 1.0e-3) {
      //            rho_ = rho_*0.5;
      //            shared_ptr<Vector> sol_tmp = make_shared<Vector>(nVar_ + 2
      //            * nCon_);
      //            myQP_->update_penalty(rho_);
      //
      //            try {
      //                myQP_->solveQP(stats_, old_options_);
      //            }
      //            catch (QP_NOT_OPTIMAL) {
      //                handle_error("QP NOT OPTIMAL");
      //            }
      //
      //            //recalculate the infeasibility measure of the model by
      //            // calculating the one norm of the slack variables
      //
      //            current_infeasibility__model_ =
      //            myQP_->get_current_infeasibility__model();

      //            p_k_->copy_vector(sol_tmp);
      //            get_trial_point_info();
      //            qp_obj = get_obj_QP();
    }
  }
}

/**
 * @brief Use the Ipopt Reference Options and set it to default values.
 */
void SqpAlgorithm::register_options_(SmartPtr<RegisteredOptions> reg_options)
{
  // Options related to output
  reg_options->SetRegisteringCategory("Output");
  reg_options->AddBoundedIntegerOption(
      "print_level", "Output verbosity level.", 0, J_LAST_LEVEL - 1,
      J_ITERSUMMARY, "Sets the default verbosity level for console output. The "
                     "larger this value the more detailed is the output.");

  reg_options->AddStringOption1(
      "output_file",
      "File name of desired output file (leave unset for no file output).",
      "sqp.opt", "*", "Any acceptable standard file name",
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

  reg_options->SetRegisteringCategory("Trust-region");
  reg_options->AddBoundedNumberOption(
      "trust_region_ratio_decrease_tol",
      "trust-region parameter for the ratio test triggering decrease.", 0.,
      true, 1., true, 0.25,
      "If ratio <= trust_region_ratio_decrease_tol, then the trust-region "
      "radius for the next "
      "iteration will be decreased for the next iteration.");
  reg_options->AddBoundedNumberOption(
      "trust_region_ratio_accept_tol",
      "trust-region parameter for the ratio test.", 0., true, 1., true, 1.0e-8,
      "The trial point will be accepted if ratio >= "
      "trust_region_ratio_accept_tol. ");
  reg_options->AddBoundedNumberOption(
      "trust_region_ratio_increase_tol",
      "trust-region parameter for the ratio test.", 0., true, 1., true, 0.75,
      "If ratio >= trust_region_ratio_increase_tol and the search direction "
      "hits the  "
      "trust-region boundary, the trust-region radius will "
      "be increased for the next iteration.");
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

  reg_options->AddLowerBoundedNumberOption("trust_region_init_value",
                                           "Initial trust-region radius value",
                                           0., true, 1.0);
  reg_options->AddLowerBoundedNumberOption(
      "trust_region_max_value", "Maximum value of trust-region radius "
                                "allowed for the radius update",
      0., true, 1e10);
  reg_options->AddLowerBoundedNumberOption(
      "trust_region_min_value", "Minimum value of trust-region radius "
                                "allowed for the radius update",
      0., true, 1e-16);

  reg_options->SetRegisteringCategory("Penalty Update");
  reg_options->AddLowerBoundedNumberOption(
      "penalty_parameter_init_value", "Initial value of the penalty parameter.",
      0., true, 1.0, "");
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
                               "Maximum value of the penalty parameter", 1.0e8);
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
  reg_options->AddNumberOption("active_set_tol", "",
                               1.0e-5); // TODO: make lower bounded options
  reg_options->AddNumberOption("opt_tol", "", 1.0e-8);
  reg_options->AddNumberOption("opt_tol_complementarity", "", 1.0e-4);
  reg_options->AddNumberOption("opt_tol_dual_feasibility", " ", 1.0e-4);
  reg_options->AddNumberOption("opt_tol_primal_feasibility", " ", 1.0e-4);
  reg_options->AddNumberOption("opt_tol_stationarity_feasibility", "", 1e-4);
  reg_options->AddNumberOption("opt_second_tol", " ", 1.0e-8);

  reg_options->AddLowerBoundedNumberOption(
      "cpu_time_limit", "CPU time limit", 0., true, 1e10,
      "Time limit measured in CPU time (in seconds)");
  reg_options->AddLowerBoundedNumberOption(
      "wallclock_time_limit", "Wallclock time limit", 0., true, 1e10,
      "Time limit measured in wallclock time (in seconds)");

  reg_options->SetRegisteringCategory("General");
  reg_options->AddNumberOption("step_size_tol",
                               "the smallest stepsize can be accepted"
                               "before concluding convergence",
                               1.0e-15);
  reg_options->AddIntegerOption("max_num_iterations",
                                "Maximum number of iteration for the algorithm",
                                3000);
  reg_options->AddStringOption2(
      "perform_second_order_correction_step",
      "Tells the algorithm to calculate the second-order correction step "
      "during the main iteration"
      "yes",
      "no", "not calculate the soc steps", "yes",
      "will calculate the soc steps", "");

  reg_options->SetRegisteringCategory("QPsolver");
  reg_options->AddIntegerOption("testOption_QP",
                                "Level of Optimality test for QP", -99);
  reg_options->AddIntegerOption("qp_solver_max_num_iterations",
                                "maximum number of iteration for the "
                                "QP solver in solving each QP",
                                1000);
  reg_options->AddIntegerOption("lp_solver_max_num_iterations",
                                "maximum number of iteration for the "
                                "LP solver in solving each LP",
                                1000);
  reg_options->AddIntegerOption("qp_solver_print_level",
                                "print level for QP solver", 0);
  reg_options->AddStringOption4(
      "qp_solver_choice", "QP solver used for step computation.", "qore",
      "qpoases", "", "qore", "", "gurobi", "", "cplex", "");

  //    reg_options->AddStringOption("QPsolverChoice",
  //		    "The choice of QP solver which will be used in the
  // Algorithm",
  //		    "qpOASES");

  reg_options->SetRegisteringCategory("LPsolver");
  //    reg_options->AddStringOption("LPsolverChoice",
  //		    "The choice of LP solver which will be used in the
  // Algorithm",
  //		    "qpOASES");

  reg_options->AddIntegerOption("testOption_LP",
                                "Level of Optimality test for LP", -99);
  reg_options->AddNumberOption("iter_malower_variable_bounds_p",
                               "maximum number of iteration for the "
                               "LP solver in solving each LP",
                               100);
  reg_options->AddNumberOption("print_level_lp", "print level for LP solver",
                               0);
}

void SqpAlgorithm::get_option_values_()
{
  options_->GetIntegerValue("max_num_iterations", max_num_iterations_, "");
  options_->GetNumericValue("cpu_time_limit", cpu_time_limit_, "");
  options_->GetNumericValue("wallclock_time_limit", wallclock_time_limit_, "");

  options_->GetNumericValue("trust_region_init_value", trust_region_init_value_,
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

  options_->GetNumericValue("penalty_parameter_init_value",
                            penalty_parameter_init_value_, "");
  options_->GetNumericValue("penalty_update_tol", penalty_update_tol_, "");
  options_->GetNumericValue("penalty_parameter_increase_factor",
                            penalty_parameter_increase_factor_, "");
  options_->GetNumericValue("penalty_parameter_max_value",
                            penalty_parameter_max_value_, "");
  options_->GetNumericValue("eps1", eps1_, "");
  options_->GetNumericValue("eps1_change_parm", eps1_change_parm_, "");
  options_->GetNumericValue("eps2", eps2_, "");
  options_->GetIntegerValue("penalty_iter_max", penalty_iter_max_, "");

  options_->GetBoolValue("perform_second_order_correction_step",
                         perform_second_order_correction_step_, "");

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
  options_->GetEnumValue("qp_solver_choice", enum_int, "");
  qp_solver_choice_ = QpSolver(enum_int);
}

void SqpAlgorithm::second_order_correction()
{
  // FIXME: check correctness
  if ((!isaccept_) && perform_second_order_correction_step_) {
    isaccept_ = false;

#ifdef DEBUG
//#ifdef CHECK_SOC
//        EJournalLevel debug_print_level = old_options_->debug_print_level;
//        SmartPtr<Journal> debug_jrnl = jnlst_->GetJournal("Debug");
//        if (IsNull(debug_jrnl)) {
//            debug_jrnl = jnlst_->AddFileJournal("Debug", "debug.out",
//            J_ITERSUMMARY);
//        }
//        debug_jrnl->SetAllPrintLevels(debug_print_level);
//        debug_jrnl->SetPrintLevel(J_DBG, J_ALL);
//
//
//        cout << "       Entering the SOC STEP calculation\n";
//        cout << "       class member isaccept is ";
//        if (isaccept_)
//            cout << "TRUE" << endl;
//        else
//            cout << "FALSE" << endl;
//        cout << "---------------------------------------------------------\n";
//
//        cout << endl;
//#endif
#endif

    shared_ptr<Vector> p_k_tmp = make_shared<Vector>(
        num_variables_); // for temporarily storing data for p_k
    p_k_tmp->copy_vector(trial_step_);

    shared_ptr<Vector> s_k =
        make_shared<Vector>(num_variables_); // for storing solution for SOC
    shared_ptr<Vector> tmp_sol =
        make_shared<Vector>(num_variables_ + 2 * num_constraints_);

    double norm_p_k_tmp = trial_step_norm_;
    double qp_obj_tmp = qp_obj_;
    shared_ptr<Vector> Hp = make_shared<Vector>(num_variables_);
    current_lagrangian_hessian_->multiply(trial_step_, Hp); // Hp = H_k*p_k
    Hp->add_vector(1., current_objective_gradient_);        //(H_k*p_k+g_k)
    qp_solver_->update_grad(Hp);
    qp_solver_->update_bounds(
        trust_region_radius_, lower_variable_bounds_, upper_variable_bounds_,
        trial_iterate_, lower_constraint_bounds_, upper_constraint_bounds_,
        trial_constraint_values_);
    trial_step_norm_ = trial_step_->calc_inf_norm();

    QpSolverExitStatus exit_status = qp_solver_->solve_qp(solver_statistics_);

    if (exit_status != QPEXIT_OPTIMAL) {
      const string& nlp_name = sqp_nlp_->get_nlp_name();
      qp_solver_->write_qp_data(nlp_name + "qpdata.log");
      assert(false && "Still need to decide how to handle QP solver error.");
      // exit_flag_ = qp_solver_->get_status();
      // break;
    }

    tmp_sol->copy_vector(qp_solver_->get_primal_solution());
    s_k->copy_vector(tmp_sol);

    qp_obj_ = get_obj_QP() +
              (qp_obj_tmp - penalty_parameter_ * current_model_infeasibility_);
    trial_step_->add_vector(1., s_k);
    get_trial_point_info();
    ratio_test();
    if (!isaccept_) {
      trial_step_->copy_vector(p_k_tmp);
      qp_obj_ = qp_obj_tmp;
      trial_step_norm_ = norm_p_k_tmp;
      qp_solver_->update_grad(current_objective_gradient_);
      qp_solver_->update_bounds(
          trust_region_radius_, lower_variable_bounds_, upper_variable_bounds_,
          current_iterate_, lower_constraint_bounds_, upper_constraint_bounds_,
          current_constraint_values_);
    }
  }
}

/**
*@brief get the objective value of QP from myQP object
*@relates QPhandler.hpp
*@return QP obejctive
*/

double SqpAlgorithm::get_obj_QP()
{
  return (qp_solver_->get_objective());
}

////////////////////////////////////////////////////////////////////////////

void SqpAlgorithm::print_final_stats()
{
  jnlst_->Printf(J_ITERSUMMARY, J_MAIN, "======================================"
                                        "===================================="
                                        "============================\n");

  // Determine string describing the exit status
  string exit_status;
  switch (exit_flag_) {
    case OPTIMAL:
      exit_status = "Optimal solution found.";
      break;
    case PRED_REDUCTION_NEGATIVE:
      exit_status = "Error: Predict reduction is negative.";
      break;
    case INVALID_NLP:
      exit_status = "Error: Invalid NLP.";
      break;
    case EXCEED_MAX_ITERATIONS:
      exit_status = "Maximum number of iterations exceeded.";
      break;
    case EXCEED_MAX_CPU_TIME:
      exit_status = "CPU time limit exceeded.";
      break;
    case EXCEED_MAX_WALLCLOCK_TIME: // TODO NEXT
      exit_status = "Wallclock time limit exceeded.";
      break;
    case TRUST_REGION_TOO_SMALL:
      exit_status = "Trust region becomes too small.";
      break;
    case QPERROR_INFEASIBLE:
      exit_status = "Error: QP solver claims that QP is infeasible.";
      break;
    case QPERROR_UNBOUNDED:
      exit_status = "Error: QP solver claims that QP is unbounded.";
      break;
    case QPERROR_EXCEED_MAX_ITER:
      exit_status = "Error: QP solver exceeded internal iteration limit.";
      break;
    case QPERROR_UNKNOWN:
      exit_status = "Error: Unknown QP solver error.";
      break;
#if 0
    case QP_OPTIMAL:
    case QPERROR_NOTINITIALISED:
    case CONVERGE_TO_NONOPTIMAL:
    case QPERROR_PREPARINGAUXILIARYQP:
    case QPERROR_AUXILIARYQPSOLVED:
    case QPERROR_PERFORMINGHOMOTOPY:
    case QPERROR_HOMOTOPYQPSOLVED:
      exit_status = "Should not appear as exit code.";
      break;
#endif
    default:
      exit_status =
          "Error: exit_flag has uncaught value " + to_string(exit_flag_) + ".";
      break;
  }

  // Print the exit status
  jnlst_->Printf(
      J_ITERSUMMARY, J_MAIN,
      "Exit status:                                                %23s\n",
      exit_status.c_str());

  jnlst_->Printf(
      J_ITERSUMMARY, J_MAIN,
      "Number of Variables                                         %23i\n",
      num_variables_);
  jnlst_->Printf(
      J_ITERSUMMARY, J_MAIN,
      "Number of Constraints                                       %23i\n",
      num_constraints_);
  jnlst_->Printf(
      J_ITERSUMMARY, J_MAIN,
      "Iterations:                                                 %23i\n",
      solver_statistics_->num_sqp_iterations_);
  jnlst_->Printf(
      J_ITERSUMMARY, J_MAIN,
      "QP Solver Iterations:                                       %23i\n",
      solver_statistics_->num_qp_iterations_);
  jnlst_->Printf(
      J_ITERSUMMARY, J_MAIN,
      "Final Objectives:                                           %23.16e\n",
      obj_value_);
  jnlst_->Printf(
      J_ITERSUMMARY, J_MAIN,
      "Primal Feasibility Violation                                %23.16e\n",
      opt_status_.primal_violation);

  jnlst_->Printf(
      J_ITERSUMMARY, J_MAIN,
      "Dual Feasibility Violation                                  %23.16e\n",
      opt_status_.dual_violation);
  jnlst_->Printf(
      J_ITERSUMMARY, J_MAIN,
      "Complmentarity Violation                                    %23.16e\n",
      opt_status_.compl_violation);
  jnlst_->Printf(
      J_ITERSUMMARY, J_MAIN,
      "Stationarity Violation                                      %23.16e\n",
      opt_status_.stationarity_violation);
  jnlst_->Printf(
      J_ITERSUMMARY, J_MAIN,
      "||p_k||                                                     %23.16e\n",
      trial_step_norm_);
  jnlst_->Printf(
      J_ITERSUMMARY, J_MAIN,
      "||c_k||                                                     %23.16e\n",
      current_infeasibility_);

  jnlst_->Printf(J_ITERSUMMARY, J_MAIN, DOUBLE_LONG_DIVIDER);
}

} // END_NAMESPACE_SQPHOTSTART
