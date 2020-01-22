/* Copyright (C) 2019
* All Rights Reserved.
*
* Authors: Xinyi Luo
* Date:2019-05
*/
#ifndef SQPHOTSTART_STATISTICS_HPP
#define SQPHOTSTART_STATISTICS_HPP

namespace RestartSqp {

class Statistics
{
public:
  /** Constructor*/
  Statistics()
  {
    reset();
  }

  /** Destructor*/
  ~Statistics()
  {
  }

  /** Reset all counters to 0. */
  void reset()
  {
    num_qp_iterations_ = 0;
    num_sqp_iterations_ = 0;
    num_trial_penalty_parameters_ = 0;
    num_penalty_parameter_failed_ = 0;
    num_penalty_parameter_increase_ = 0;
    num_second_order_corrections_ = 0;
  }

  /* Increase the SQP iteration counter by one. */
  inline void increase_sqp_iteration_counter()
  {
    num_sqp_iterations_++;
  }

  /* Increase the QP iteration counter by n. */
  inline void increase_qp_iteration_counter(int n = 1)
  {
    num_qp_iterations_ += n;
  }

  /* Increase the counter that keeps track of how many times a trial penalty
   * parameter was not accepted. */
  inline void trial_penalty_parameter_not_accepted()
  {
    num_penalty_parameter_failed_++;
  }

  /* Increase the counter that keeps track of how many times a new penalty
   * parameter values was tried. */
  void try_new_penalty_parameter()
  {
    num_trial_penalty_parameters_++;
  }

  /* Increase the counter that keeps track of how many times a new penalty
   * parameter values was accepted. */
  inline void penalty_parameter_increased()
  {
    num_penalty_parameter_increase_++;
  }

  /* Increase the counter for second order corrections. */
  inline void second_order_correction()
  {
    num_second_order_corrections_++;
  }

  /* Member Variables */
public:
  /** Number of QP iterations. */
  int num_qp_iterations_;

  /** Number of SQP iterations. */
  int num_sqp_iterations_;

  /** Number of times a new penalty parameter was tried. */
  int num_trial_penalty_parameters_;

  /** Number of times the penalty parameter was increased. */
  int num_penalty_parameter_failed_;

  /** Number of times a new penalty parameter was accepted. */
  int num_penalty_parameter_increase_;

  /** Number of second order correction steps. */
  int num_second_order_corrections_;
};

} // END_NAMESPACE_SQPHOTSTART

#endif /* SQPHOTSTART_STATS_HPP*/
