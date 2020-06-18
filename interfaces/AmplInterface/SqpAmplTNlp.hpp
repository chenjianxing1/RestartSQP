/* Copyright (C) 2019
 * All Rights Reserved.
 *
 * Authors: Xinyi Luo
 * Date:2019-06
 */
#ifndef SQPHOTSTART_SQPAMPLTNLP_HPP
#define SQPHOTSTART_SQPAMPLTNLP_HPP

#include "IpoptAmplTNlp.hpp"
#include "restartsqp/SqpIpoptTNlp.hpp"

namespace RestartSqp {
/**
 * This is part of SQPhotstart
 *
 * This makes an Ipopt::TNLP look like an SqpTNlp.
 *
 * IMPORTANT: The Lagrangian function here is defined as L(x,l) = f(x) - sum_i
 * l_ic_j(c)
 *
 *            THIS IS DIFFERENT FROM IPOPT's DEFINITION!!
 */
class SqpAmplTNlp : public SqpIpoptTNlp
{

public:
  /** @brief Constructor that an instance of Ipopt's TNLP */
  SqpAmplTNlp(Ipopt::SmartPtr<IpoptAmplTNlp> ipopt_ampl_tnlp,
              const std::string& nlp_name = "AMPL NLP")
   : SqpIpoptTNlp(ipopt_ampl_tnlp, nlp_name)
   , ipopt_ampl_tnlp_(ipopt_ampl_tnlp)
  {}

  /** Destructor*/
  ~SqpAmplTNlp() {}

  /*
   * @brief Get the starting point from the NLP object.
   */
  bool get_starting_point(int num_variables, bool init_primal_variables,
                          double* primal_variables, bool init_bound_multipliers,
                          double* bound_multipliers, int num_constraints,
                          bool init_constraint_multipliers,
                          double* constraint_multipliers) override
  {
    return ipopt_ampl_tnlp_->get_starting_point(
        num_variables, init_primal_variables, primal_variables,
        init_bound_multipliers, bound_multipliers, num_constraints,
        init_constraint_multipliers, constraint_multipliers);
  }

  /**
   * @brief Return the results of the optimization run to the user.
   */
  void finalize_solution(SqpSolverExitStatus status, int num_variables,
                         const double* primal_solution,
                         const double* bound_multipliers,
                         const ActivityStatus* bound_activity_status,
                         int num_constraints, const double* constraint_values,
                         const double* constraint_multipliers,
                         const ActivityStatus* constraint_activity_status,
                         double objective_value,
                         std::shared_ptr<const Statistics> stats) override
  {
    return ipopt_ampl_tnlp_->finalize_solution(
        status, num_variables, primal_solution, bound_multipliers,
        bound_activity_status, num_constraints, constraint_values,
        constraint_multipliers, constraint_activity_status, objective_value,
        stats);
  }

  /** No initial working set is available from an Ipopt TNLP */
  bool use_initial_working_set() override
  {
    return ipopt_ampl_tnlp_->use_initial_working_set();
  }

  /** Get the initial working set. */
  bool get_initial_working_sets(
      int num_variables, ActivityStatus* bounds_working_set,
      int num_constraints, ActivityStatus* constraints_working_set) override
  {
    return ipopt_ampl_tnlp_->get_initial_working_sets(
        num_variables, bounds_working_set, num_constraints,
        constraints_working_set);
  }

  /** Get the initial value of the penalty parameter. */
  double get_initial_penalty_parameter()
  {
    return ipopt_ampl_tnlp_->get_initial_penalty_parameter();
  }

private:
  /** @name Hide unused default methods. */
  //@{
  /** Default constructor*/
  SqpAmplTNlp();

  /** Copy Constructor */
  SqpAmplTNlp(const SqpAmplTNlp&);

  /** Overloaded Equals Operator */
  void operator=(const SqpAmplTNlp&);
  //@}

  /** Pointer to the object that can get the SQP specific information from the
   * AMPL interface. */
  Ipopt::SmartPtr<IpoptAmplTNlp> ipopt_ampl_tnlp_;
};

} // namespace RestartSqp

#endif // SQPHOTSTART_SqpIpoptNlp_HPP
