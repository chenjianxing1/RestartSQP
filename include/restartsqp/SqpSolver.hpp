/* Copyright (C) 2019
 * All Rights Reserved.
 *
 * Authors: Xinyi Luo
 * Date:    2019-07
 */
#ifndef SQPSOLVER_HPP_
#define SQPSOLVER_HPP_

#include "restartsqp/Types.hpp"
#include "IpOptionsList.hpp"
#include "restartsqp/QpHandler.hpp"
#include "restartsqp/Statistics.hpp"

namespace RestartSqp {
/**
 *
 * @brief This is the class with method solve a NLP problem by using SQP(SL1QP)
 *
 *It can solve a problem in the following format
 *
 *	minimize 	f(x)
 *	subject   c_l<=c(x)<=c_u,
 *	      	  x_l<= x  <=x_u
 *
 *
 *To use this method, call @Optimize and input NLP class object.
 *
 */
class SqpSolver
{
  ///////////////////////////////////////////////////////////
  //                      PUBLIC METHODS                   //
  ///////////////////////////////////////////////////////////
public:
  /** @name constructor/destructor*/
  //@{
  /** Default Constructor*/
  SqpSolver();

  /** Constructor that can take in preinitialized options and journalist from
   * Ipopt. */
  SqpSolver(Ipopt::SmartPtr<Ipopt::RegisteredOptions> reg_options,
            Ipopt::SmartPtr<Ipopt::OptionsList> options,
            Ipopt::SmartPtr<Ipopt::Journalist> jnlst);

  /** Destructor*/
  ~SqpSolver();
  //@}

  /** Return journalist.  This can be used by caller. */
  Ipopt::SmartPtr<Ipopt::Journalist> get_jnlst()
  {
    return jnlst_;
  }

  /** Return options list.  This can be used by caller to set algorithmic
   * options. */
  Ipopt::SmartPtr<Ipopt::OptionsList> get_options_list()
  {
    return options_;
  }

  /** Return the list of registered options.  This can be used to add more
   * options that can be used by the called. */
  Ipopt::SmartPtr<Ipopt::RegisteredOptions> get_registered_options()
  {
    return reg_options_;
  }

  //@}

  /** @name Solution methods. */
  //@{
  /**
   * @brief This is the main method to optimize the NLP given as the input
   * for the first time.
   *
   * @param sqp_nlp: nlp object defining the NLP to be solved
   * options_file_name: Name of the file with the options
   * keep_output_file: If true, no new output file will be opened
   */
  void optimize_nlp(std::shared_ptr<SqpTNlp> sqp_tnlp,
                    const std::string& options_file_name = "sqp.opt",
                    bool keep_output_file = false);

  /**
   * @brief This is the method to reoptimize an NLP.  This method
   * remembers some algorithmic quantities from the previous optimization,
   * such as penalty parameter value.  It also retains the options from the
   * previous solve.  But it does not assume that the structure or data of
   * the NLP is the same as from the previous call.
   *
   * @param sqp_nlp: nlp object defining the NLP to be solved
   */
  void reoptimize_nlp(std::shared_ptr<SqpTNlp> sqp_tnlp);

  //@}

  /** @name Getters*/
  //@{
  inline SqpSolverExitStatus get_exit_flag() const
  {
    return exit_flag_;
  }

  /** This is the internal KKT error, using the internal scaling. */
  inline KktError get_kkt_error() const
  {
    return current_kkt_error_;
  }

  /** After successful termination, this is the optimal objective function
   * value.
   *  This is the value consistent with the original problem formulation,
   * independent
   *  of any internal scaling. */
  inline double get_final_objective() const
  {
    return current_objective_value_;
  }

  inline std::shared_ptr<Statistics> get_stats() const
  {
    return solver_statistics_;
  }

  inline int get_num_constr() const
  {
    return num_constraints_;
  }

  inline int get_num_var() const
  {
    return num_variables_;
  }

  /** After successful termination, these are the optimal primal variables. */
  inline std::shared_ptr<const Vector> get_current_primal_variables() const
  {
    return current_iterate_;
  }

  /** After successful termination, these are the optimal constraint
   *  multipliers.
   *  They are scaled consistent with the original problem formulation,
   *  independent of any internal scaling. */
  inline std::shared_ptr<const Vector>
  get_current_constraint_multipliers() const
  {
    return current_constraint_multipliers_;
  }

  /** After successful termination, these are the optimal bound multipliers.
   *  They are scaled consistent with the original problem formulation,
   *  independent of any internal scaling. */
  inline std::shared_ptr<const Vector> get_current_bound_multipliers() const
  {
    return current_bound_multipliers_;
  }

  inline const ActivityStatus* get_bounds_working_set() const
  {
    return qp_solver_->get_bounds_working_set();
  }

  inline const ActivityStatus* get_constraints_working_set() const
  {
    return qp_solver_->get_constraints_working_set();
  }

  //@}
  ///////////////////////////////////////////////////////////
  //                      PRIVATE METHODS                  //
  ///////////////////////////////////////////////////////////
private:
  /**
   * @brief set the default option values
   */
  void register_options_(Ipopt::SmartPtr<Ipopt::RegisteredOptions> reg_options,
                         bool skip_ipopt_options);

  /**
   * @brief Prepare data structures before we can start the algorithm loop.
   *
   *  Create QP objects
   */
  void initialize_for_new_nlp_();

  /**
   * @brief Initialize the options
   *
   *  Get the options values.
   */
  void initialize_options_(const std::string& options_file_name,
                           bool keep_output_file);

  /** Determine starting point and compute all corresponding quanties. */
  void initialize_iterates_();

  /** Print initial output
   *
   *  This prints the problem dimensions and initial header line.
   */
  void print_initial_output_();

  /** Print iteration output. */
  void print_iteration_output_();

  /** Calculate the search direction.
   *
   *  For this, the QP is solved.
   */
  void
  calculate_search_direction_(std::shared_ptr<const Vector> qp_constraint_body);

  /**
   * @brief This is the function that checks if the current point is optimal,
   * and
   * decides if to exit the loop or not
   *
   * If it decides the function is optimal, the class member _exitflag =
   * OPTIMAL
   * if it decides that there is an error during the function run or the
   *  function cannot be solved, it will assign _exitflag the	corresponding
   *  code according to the error type.
   */
  void check_optimality_();

  /**
  * @brief This function calculates the infeasibility for given x_k and c_k with
  respect
  * to their corresponding bounds
  * @return infea_measure = ||-max(c_k-c_u),0||_1 +||-min(c_k-c_l),0||_1+
                    ||-max(x_k-x_u),0||_1 +||-min(x_k-x_l),0||_1
  */
  double compute_constraint_violation_(
      std::shared_ptr<const Vector> variable_values,
      std::shared_ptr<const Vector> constraint_values_);

  /**
   * @brief This function calculates the infeasibility measure for  current
   * iterate x_k
   *
   *   infea_measure_trial =
   * norm(-max(c_trial-cu,0),1)+norm(-min(c_trial-cl,0),1)
   *
   *
   */

  /**
   * @brief Calculate the trial point based on current search direction,
   * x_trial = x_k+p_k,
   *       and get the funcion value, constraints value, and infeasibility
   * measure
   *       at the trial point
   */
  void calc_trial_point_and_values_();
  ;

  /**
   * @brief calculate the second order correction step and decide if accept the
   * calculated step or not.
   * It will be calculated only if the second_order_correction in options is set
   * to be true.
   *
   */

  void second_order_correction_();

  /**
   *
   * @brief This function performs the ratio test to determine if we should
   * accept
   * the trial point
   *
   * The ratio is calculated by
   * (P_1(x_k;\rho)-P_1( x_trial;\rho))/(q_k(0;\rho)-q_k(p_k;rho), where
   * P_1(x,rho) = f(x) + rho* infeasibility_measure is the l_1 merit function
   * and
   * q_k(p; rho) = f_k+ g_k^Tp +1/2 p^T H_k p+rho* infeasibility_measure_model
   * is
   * the quadratic model at x_k.
   * The trial point  will be accepted if the ratio >= eta_s.
   * If it is accepted, the function will also updates the gradient, Jacobian
   * information by reading from nlp_ object. The corresponding flags of class
   * member QPinfoFlag_ will set to be true.
   */
  void perform_ratio_test_();

  /**
   * @brief Update the trust region radius.
   *
   * This function update the trust-region radius when the ratio calculated by
   * the
   * ratio test is smaller than eta_c or bigger than eta_e and the
   * search_direction
   * hits the trust-region bounds.
   * If ratio<eta_c, the trust region radius will decrease by the parameter
   * gamma_c, to be gamma_c* delta_
   * If ratio > eta_e and delta_= norm_p_k_ the trust-region radius will be
   * increased by the parameter gamma_c.
   *
   * If trust region radius has changed, the corresponding flags will be set to
   * be
   * true;
   */
  void update_trust_region_radius_();

  /** Copy trial quantities to current quantities and evaluate gradients at new
   * iterate
   */
  void accept_trial_point_();

  /**
   * @brief Update the penalty parameter
   */
  void update_penalty_parameter_();

  /**
   * @brief This function will set up the data for the QP subproblem
   *
   * It will initialize all the data at once at the beginning of the Algorithm.
   * After
   * that, the data in the QP problem will be updated according to the class
   * member qp_update_tracker_
   */
  void setup_qp_(std::shared_ptr<const Vector> qp_constraint_body);

  /**
   * @brief This function will setup the data for the LP subproblem
   */
  void setup_lp_();

  /** Calculate constraint violated of the model in l1-norm, for a given step.
   */
  double calc_model_infeasibility_(std::shared_ptr<const Vector> step);

  /** Update the predicted reduction of the merit function for the current trial
   * step and penalty parameter */
  void update_predicted_reduction_();

  /** Return a small number close to machine precision relative to current data
   * that is added the quantities to prevent wrong decisions due to numerical
   * error. */
  double numerical_error_buffer_();

  /** Increase the penalty parameter and compute the new candidate step and the
   * corresponding violation of the linearized constraints. */
  void increase_penalty_parameter_();

  /** Store a backup of the relevant values from the current iteration (after QP solve and computation of predicted reduction).  This is used for the watchdog technique. */
  void store_watchdog_backups_();

  /** Restore a backup of the relevant values from the current iteration (after QP solve and computation of predicted reduction).  This is used for the watchdog technique. */
  void restore_watchdog_backups_();

  /** Delete backup vectors and matrices stored by store_watchdog_backups_(); */
  void delete_watchdog_backups_();

  /**
   * @brief alloocate memory for class members.
   * This function initializes all the shared pointer which will be used in the
   * Algorithm::Optimize, and it copies all parameters that might be changed
   * during
   * the run of the function Algorithm::Optimize.
   */
  void allocate_memory_();

  /** Free memory of any allocated arrays. */
  void free_memory_();

  /**
   *
   * @brief This function checks how each constraint specified by the nlp
   * readers are
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
  void classify_constraints_types_();

  void print_final_stats_();

  /** Call the finalize solution method to return results to the user. */
  void return_results_();
  //@}

  /** @name Wrappers for function evaluations that take care of scaling. */
  //@{
  /** Objective function. */
  bool eval_f_(std::shared_ptr<const Vector> x, double& obj_value);

  /** Constraints values. */
  bool eval_constraints_(std::shared_ptr<const Vector> x,
                         std::shared_ptr<Vector> constraints);

  /** Objective gradient. */
  bool eval_gradient_(std::shared_ptr<const Vector> x,
                      std::shared_ptr<Vector> gradient);

  /** Get the matrix structure of the Jacobian */
  bool get_jacobian_structure_(std::shared_ptr<const Vector> x,
                               std::shared_ptr<SparseTripletMatrix> Jacobian);

  /** Evaluate Jacobian at point x */
  bool eval_jacobian_(std::shared_ptr<const Vector> x,
                      std::shared_ptr<SparseTripletMatrix> Jacobian);

  /** Get the structure of the Hessian */
  bool get_hessian_structure_(std::shared_ptr<const Vector> x,
                              std::shared_ptr<const Vector> lambda,
                              std::shared_ptr<SparseTripletMatrix> Hessian);

  /** Evaluate Hessian of Lagragian function at  (x, lambda) */
  bool eval_hessian_(std::shared_ptr<const Vector> x,
                     std::shared_ptr<const Vector> lambda,
                     std::shared_ptr<SparseTripletMatrix> Hessian);
  //@}

  ///////////////////////////////////////////////////////////
  //              PRIVATE CLASS MEMBERS                    //
  ///////////////////////////////////////////////////////////
private:
  /** @name Hide unused default methods. */
  //@{
  /** Copy Constructor */
  SqpSolver(const SqpSolver&);

  /** Overloaded Equals Operator */
  void operator=(const SqpSolver&);
  //@}

  /** Get the option values from the options object. */
  void get_option_values_();

  /** SqpNlp object that evaluates all NLP information. */
  std::shared_ptr<SqpNlp> sqp_nlp_;

  /** Options and output infrastructure.
   *
   *  We use the tools from the Ipopt package for now.  The annoying
   *  thing is that we need to use Ipopt's SmartPtr's. */
  //@{
  /** Journal through which all output is directed. */
  Ipopt::SmartPtr<Ipopt::Journalist> jnlst_;
  /** Name of the currently opened output file. */
  std::string current_output_file_name_;

  /** List of options.
   *
   *  It is initialized in the constructor but then be retrieved by
   *  user of this object and changed.  With a few exceptions, the
   *  option values are retrieved when the Initialize method is
   *  called. */
  Ipopt::SmartPtr<Ipopt::OptionsList> options_;

  /** List of registered options.
   *
   *  This contains all options with default values, etc.
   */
  Ipopt::SmartPtr<Ipopt::RegisteredOptions> reg_options_;
  //@}

  /** Lower bounds for the variables. */
  std::shared_ptr<Vector> lower_variable_bounds_;
  /** Upper bounds for the variables. */
  std::shared_ptr<Vector> upper_variable_bounds_;
  /** Lower bounds for the constraints. */
  std::shared_ptr<Vector> lower_constraint_bounds_;
  /** Upper bounds for the constraints. */
  std::shared_ptr<Vector> upper_constraint_bounds_;

  /** Indicates for each bound whether it has lower and/or upper bounds. */
  ConstraintType* bound_type_;
  /** Indicates for each constraint whether it has lower and/or upper bounds. */
  ConstraintType* constraint_type_;

  /** Problem name. */
  std::string problem_name_;

  /** Number of variables. */
  int num_variables_;
  /** Number of constraints. */
  int num_constraints_;
  /** Number of equality contraints (subset of the constraints). */
  int num_equality_constraints_;

  /** Flag indicating whether this is the formulation that permits the bounds to
   *  be violated. */
  bool slack_formulation_;

  /** Exit flag.  If set to UNKNOWN the algorithm loop should still progress. */
  SqpSolverExitStatus exit_flag_;

  /** Primal current iterate. */
  std::shared_ptr<Vector> current_iterate_;
  /** Current estimate of the constraint multipliers.
   *
   *  During the SQP algorithm, these are the internal scaled values.  Once
   *  the algorithm terminated, these are the unscaled values, consistent with
   *  the original problem formulation.
   */
  std::shared_ptr<Vector> current_constraint_multipliers_;
  /** Current estimate of the bound multipliers.
  *
  *  During the SQP algorithm, these are the internal scaled values.  Once
  *  the algorithm terminated, these are the unscaled values, consistent with
  *  the original problem formulation.
  */
  std::shared_ptr<Vector> current_bound_multipliers_;

  /** Gradient of the objective function at the current iterate. */
  std::shared_ptr<Vector> current_objective_gradient_;
  /** Value of the constraint bodies at the current iterate. */
  std::shared_ptr<Vector> current_constraint_values_;
  /** Constraint Jacobian at the current iterate. */
  std::shared_ptr<SparseTripletMatrix> current_constraint_jacobian_;
  /** Hessian of the Lagrangian function at the current primal-dual iterate. */
  std::shared_ptr<SparseTripletMatrix> current_lagrangian_hessian_;

  /** Constraint violation at current iterate in l1-norm. */
  double current_infeasibility_;
  /** Constraint violation at trial iterate in l1-norm. */
  double trial_infeasibility_;
  /** Constraint violation of trial step in linear constraint model. */
  double trial_model_infeasibility_;

  /** KKT error for the current iterate. */
  KktError current_kkt_error_;

  /** Current value of the trust region radius. */
  double trust_region_radius_;
  /** Actual reduction in the penalty function with trial step. */
  double actual_reduction_;
  /** Predicted reduction in penalty function with trial step. */
  double predicted_reduction_;

  /** Value of the objective function at the current iterate x_k */
  double current_objective_value_;

  /** Current value of the penalty parameter. */
  double current_penalty_parameter_;

  /** Object to solve LPs. */
  std::shared_ptr<QpHandler> lp_solver_;
  /** Object to solve step computation QPs. */
  std::shared_ptr<QpHandler> qp_solver_;
  /** Object to track which QP data needs to be updated. */
  QpUpdateTracker qp_update_tracker_;
  /** Flag indicating whether the step computation QP has been initialized. */
  bool qp_solver_initialized_;
  /** Store the penalty parameter value used in the most recent QP solve. */
  double last_qp_penalty_parameter_;

  /** Flag indicating whether the trial point is accepted by ratio test. */
  bool trial_point_is_accepted_;

  /** Flag indicating in which phase the watchdog strategy is. */
  char watchdog_status_;

  /** Object that collects the solver statistics (iteration counts, etc) */
  std::shared_ptr<Statistics> solver_statistics_;

  /** Trial step computed by QP. */
  std::shared_ptr<Vector> trial_step_;
  /** Trial iterate */
  std::shared_ptr<Vector> trial_iterate_;
  /** Trial constraint multipliers. */
  std::shared_ptr<Vector> trial_constraint_multipliers_;
  /** Trial bound multipliers. */
  std::shared_ptr<Vector> trial_bound_multipliers_;
  /** Value of the objective function at the trial point x_k */
  double trial_objective_value_;
  /** Constraint values at trial point. */
  std::shared_ptr<Vector> trial_constraint_values_;

  /** CPU time at the beginning of the optimization run. */
  double cpu_time_at_start_;
  /** wallclock_ time at the beginning of the optimization run. */
  double wallclock_time_at_start_;

  /** @name Watchbog backup. */
  //@{
  std::shared_ptr<Vector> backup_current_iterate_;
  std::shared_ptr<Vector> backup_current_constraint_multipliers_;
  std::shared_ptr<Vector> backup_current_bound_multipliers_;

  double backup_current_objective_value_;
  std::shared_ptr<Vector> backup_current_objective_gradient_;
  std::shared_ptr<Vector> backup_current_constraint_values_;
  std::shared_ptr<SparseTripletMatrix> backup_current_constraint_jacobian_;
  std::shared_ptr<SparseTripletMatrix> backup_current_lagrangian_hessian_;

  double backup_current_infeasibility_;
  double backup_predicted_reduction_;

  std::shared_ptr<Vector> backup_trial_step_;
  std::shared_ptr<Vector> backup_trial_constraint_multipliers_;
  std::shared_ptr<Vector> backup_trial_bound_multipliers_;
  double backup_trial_model_infeasibility_;

  //@}

  /** @namd Algorithmic option values */
  //@{
  /** Maximum number of iterations. */
  int max_num_iterations_;
  /** Maximum cpu time for one solve (in seconds). */
  double cpu_time_limit_;
  /** Maximum wallclock time for one solve (in seconds). */
  double wallclock_time_limit_;

  /** Initial trust region radius. */
  double trust_region_init_value_;
  /** Maximal trust region radius. */
  double trust_region_max_value_;
  /** Minimum trust region radius. */
  double trust_region_min_value_;
  /** Trust region ratio decrease tolerance. */
  double trust_region_ratio_decrease_tol_;
  /** Trust region ratio acceptance tolerance. */
  double trust_region_ratio_accept_tol_;
  /** Trust region ratio increase tolerance. */
  double trust_region_ratio_increase_tol_;
  /** Factor by which the trust region is reduced. */
  double trust_region_decrease_factor_;
  /** Factor by which the trust region is increased. */
  double trust_region_increase_factor_;
  /** Flag indicating whether every trial step should be accepted. */
  bool disable_trust_region_;

  /** Initial penalty parameter value. */
  double penalty_parameter_init_value_;
  /** Some penatly parameter update tolernace. */ // TODO find out what this is
  double penalty_update_tol_;
  /** Factor by which penalty parameter is increased. */
  double penalty_parameter_increase_factor_;
  /** Maximum value of penalty parameter. */
  double penalty_parameter_max_value_;
  /** Some constant in penalty parameter update rule. */
  double eps1_;
  /** Some constant in penalty parameter update rule. */
  double eps1_change_parm_;
  /** Some constant in penalty parameter update rule. */
  double eps2_;
  /** Maximum number of penalty parameter updates. ??? */ // TODO
  int penalty_iter_max_;

  /** Flag indicating whether we do a second-order correction. */
  bool perform_second_order_correction_step_;

  /** Objective scaling factor. */
  double objective_scaling_factor_;

  /** Tolerance for active set identification. */ // TODO: This should be
                                                  // separate for constraints,
                                                  // relative, etc?
  double active_set_tol_;
  /** Overall optimality tolerance. */
  double opt_tol_;
  /** Optimality tolerance for primal feasibility. */
  double opt_tol_primal_feasibility_;
  /** Optimality tolerance for dual feasibility. */
  double opt_tol_dual_feasibility_;
  /** Optimality tolerance for stationarity feasibility. */ // TODO what is
                                                            // this?
  double opt_tol_stationarity_feasibility_;
  /** Optimality tolerance for complementarity. */
  double opt_tol_complementarity_;

  /** CPU time limit for one optimization run. */
  double max_cputime_;

  /** QP solver usde for ??? */
  QpSolver qp_solver_choice_;
  //@}

}; // END_OF_ALG_CLASS

} // END_NAMESPACE_SQPHOTSTART

#endif
