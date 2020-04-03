/* Copyright (C) 2019
 * All Rights Reserved.
 *
 * Authors: Xinyi Luo
 * Date:    2019-07
 */
#ifndef LAZYSQPSOLVER_HPP_
#define LAZYSQPSOLVER_HPP_

#include "sqphot/CrossoverSqpSolver.hpp"
#include "sqphot/LazySqpTNlp.hpp"

namespace RestartSqp {
/**
 *
 * @brief This is an implementation of a lazy NLP solver that solve a problem by
 *a sequence of Sqp solves
 * by adding an increasing number of violated constraints.
 *
 * It can solve a problem in the following format
 *
 *	minimize 	f(x)
 *	subject   c_l<=c(x)<=c_u,
 *	      	  x_l<= x  <=x_u
 *
 */
class LazySqpSolver
{
  ///////////////////////////////////////////////////////////
  //                      PUBLIC METHODS                   //
  ///////////////////////////////////////////////////////////
public:
  /** @name constructor/destructor*/
  //@{
  /** Default Constructor*/
  LazySqpSolver();

#if 0
  /** Constructor that can take in preinitialized options and journalist from
   * Ipopt. */
  LazySqpSolver(Ipopt::SmartPtr<Ipopt::RegisteredOptions> reg_options,
                Ipopt::SmartPtr<Ipopt::OptionsList> options,
                Ipopt::SmartPtr<Ipopt::Journalist> jnlst);
#endif
  /** Destructor*/
  ~LazySqpSolver();
  //@}

  /** Return journalist.  This can be used by caller. */
  Ipopt::SmartPtr<Ipopt::Journalist> get_jnlst()
  {
    return crossover_sqp_solver_->get_jnlst();
  }

  /** Return options list.  This can be used by caller to set algorithmic
   * options. */
  Ipopt::SmartPtr<Ipopt::OptionsList> get_options_list()
  {
    return crossover_sqp_solver_->get_options_list();
  }
  //@}

  /** @name Solution methods. */
  //@{
  /**
   * @brief Solve the NLP in a lazy fashion.
   *
   * We called provides an initial set of constraints that should be considered.
   * Counting of constraints starts at 1.
   */
  void optimize_nlp(std::shared_ptr<SqpTNlp> sqp_tnlp,
                    int num_initial_constraints,
                    int* initial_constraint_indices,
                    const std::string& options_file_name = "sqp.opt");
  //@}

  /** @name Getters*/
  //@{
  inline SqpSolverExitStatus get_exit_flag() const
  {
    return crossover_sqp_solver_->get_exit_flag();
  }

  inline KktError get_kkt_error() const
  {
    return crossover_sqp_solver_->get_kkt_error();
  }

  inline double get_final_objective() const
  {
    return crossover_sqp_solver_->get_final_objective();
  }

  inline std::shared_ptr<Statistics> get_stats() const
  {
    return crossover_sqp_solver_->get_stats();
  }

  inline std::shared_ptr<const Vector>
  get_current_constraint_multipliers() const;

  inline std::shared_ptr<const Vector> get_current_bound_multipliers() const;
  //@}
  ///////////////////////////////////////////////////////////
  //                      PRIVATE METHODS                  //
  ///////////////////////////////////////////////////////////
private:
  /** Copy Constructor */
  LazySqpSolver(const LazySqpSolver&);

  /** Overloaded Equals Operator */
  void operator=(const LazySqpSolver&);

  /** Update the list of considered constraints. It returns the number of added
   *  constraints.  If zero, nothing was added.  A negative number indicates and error. */
  int update_considered_constraints_();

  ///////////////////////////////////////////////////////////
  //              PRIVATE CLASS MEMBERS                    //
  ///////////////////////////////////////////////////////////
private:
  /** SQP solver that will be used to solve the initial problem with Ipopt and
   * iterate with new sets of constraints.  */
  std::shared_ptr<CrossoverSqpSolver> crossover_sqp_solver_;

  /** Journalist for output */
  Ipopt::SmartPtr<Ipopt::Journalist> jnlst_;

  /** Original NLP. */
  std::shared_ptr<SqpTNlp> orig_sqp_tnlp_;

  /** Lazy constraints NLP. */
  std::shared_ptr<LazySqpTNlp> lazy_sqp_tnlp_;

  /** Number of variables in the original NLP. */
  int num_orig_variables_;

  /** Number of constraints in the original NLP. */
  int num_orig_constraints_;

  /** Number of constraints considered in the lazy NLP. */
  int num_considered_constraints_;

  /** Indices of the constraints considered in the lazy NLP.  Counting starts at
   * 1. */
  int* considered_constraints_indices_;

  /** Lower bounds for the original problem. */
  double* orig_constraint_lower_bounds_;

  /** Upper bounds for the original problem. */
  double* orig_constraint_upper_bounds_;

  /** Exit flag from the most recent NLP solve. */
  SqpSolverExitStatus exit_flag_;

  /** @name Algorithmic options. */
  //@{
  // Maximum number of Lazy NLP solves after the initial solve.
  int num_max_lazy_nlp_solves_;

  // Maximum number of constraints to be added per outer iteration
  int max_add_lazy_constraints_;
  //@}

}; // END_OF_ALG_CLASS

} // END_NAMESPACE_SQPHOTSTART

#endif
