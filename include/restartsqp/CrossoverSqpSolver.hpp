/* Copyright (C) 2019
 * All Rights Reserved.
 *
 * Authors: Xinyi Luo
 * Date:    2019-07
 */
#ifndef CROSSOVERSQPSOLVER_HPP_
#define CROSSOVERSQPSOLVER_HPP_

#include "IpIpoptApplication.hpp"
#include "restartsqp/IpoptSqpNlp.hpp"
#include "restartsqp/SqpSolver.hpp"

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
class CrossoverSqpSolver
{
  ///////////////////////////////////////////////////////////
  //                      PUBLIC METHODS                   //
  ///////////////////////////////////////////////////////////
public:
  /** @name constructor/destructor*/
  //@{
  /** Default Constructor*/
  CrossoverSqpSolver();

  /** Destructor*/
  ~CrossoverSqpSolver();
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
    return sqp_solver_->get_options_list();
  }

  /** Return the options list for Ipopt. */
  Ipopt::SmartPtr<Ipopt::OptionsList> get_ipopt_options_list()
  {
    return ipopt_app_->Options();
  }

  /**
   * Solve a new NLP, first with interior point method, then after cross-over
   * with active-set method.
   */
  void initial_solve(std::shared_ptr<SqpTNlp> sqp_tnlp,
                     const std::__cxx11::string& options_file_name = "sqp.opt");

  /**
   * ReOptimize a problem by warm-starting from the old optimal solution and
   * given active set.
   *
   */
  void next_solve(std::shared_ptr<SqpTNlp> sqp_tnlp);
  //@}

  /** @name Getters*/
  //@{
  inline SqpSolverExitStatus get_exit_flag() const
  {
    return exit_flag_;
  }

  inline KktError get_kkt_error() const
  {
    return final_kkt_error_;
  }

  inline double get_final_objective() const
  {
    return final_objective_value_;
  }

  inline std::shared_ptr<Statistics> get_stats() const
  {
    return solver_statistics_;
  }

  inline int get_num_constr() const
  {
    return sqp_solver_->get_num_constr();
  }

  inline int get_num_var() const
  {
    return sqp_solver_->get_num_var();
  }

  inline std::shared_ptr<const Vector>
  get_current_primal_variables() const
  {
    return sqp_solver_->get_current_primal_variables();
  }

  inline std::shared_ptr<const Vector>
  get_current_constraint_multipliers() const
  {
    return sqp_solver_->get_current_constraint_multipliers();
  }

  inline std::shared_ptr<const Vector> get_current_bound_multipliers() const
  {
    return sqp_solver_->get_current_bound_multipliers();
  }

  inline const ActivityStatus* get_bounds_working_set() const
  {
    return sqp_solver_->get_bounds_working_set();
  }

  inline const ActivityStatus* get_constraints_working_set() const
  {
    return sqp_solver_->get_constraints_working_set();
  }

  //@}
  ///////////////////////////////////////////////////////////
  //                      PRIVATE METHODS                  //
  ///////////////////////////////////////////////////////////
private:
  /** Register additional options for the restart solver. */
  void register_options_(Ipopt::SmartPtr<Ipopt::RegisteredOptions> reg_options);

  /** Determine the set of active variable bounds and constraints.
   *
   *  This is the cross-over method.  This version is based on the ratio of the
   * Ipopt slacks and multipliers. */
  void determine_activities_(ActivityStatus* bound_activity_status,
                             ActivityStatus* constraint_activity_status,
                             IpoptSqpNlp& ipopt_tnlp);

  ///////////////////////////////////////////////////////////
  //              PRIVATE CLASS MEMBERS                    //
  ///////////////////////////////////////////////////////////
private:
  /** @name Hide unused default methods. */
  //@{
  /** Copy Constructor */
  CrossoverSqpSolver(const CrossoverSqpSolver&);

  /** Overloaded Equals Operator */
  void operator=(const CrossoverSqpSolver&);
  //@}

  /** SQP solver for the active-set solves. */
  std::shared_ptr<SqpSolver> sqp_solver_;

  /** Ipopt solver for the interior-point parts. */
  Ipopt::SmartPtr<Ipopt::IpoptApplication> ipopt_app_;

  /** Journal through which all output is directed. */
  Ipopt::SmartPtr<Ipopt::Journalist> jnlst_;

  /** Exit flag.  If set to UNKNOWN the algorithm loop should still progress. */
  SqpSolverExitStatus exit_flag_;

  /** Number of variables in the current NLP. */
  int num_variables_;

  /** Number of constraints in the current NLP. */
  int num_constraints_;

  /** KKT error from most recent solve. */
  KktError final_kkt_error_;

  /** Final objective function value from most recent solve. */
  double final_objective_value_;

  /** Solver statistics from most recent solve. */
  std::shared_ptr<Statistics> solver_statistics_;

}; // END_OF_ALG_CLASS

} // END_NAMESPACE_SQPHOTSTART

#endif
