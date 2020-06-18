// Copyright (C)
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// This code is derived from AmplTNLP.hpp from the Ipopt project.
//
// Authors:

#ifndef SQPHOTSTART_IPOPTAMPLTNLP_HPP
#define SQPHOTSTART_IPOPTAMPLTNLP_HPP

#include "AmplTNLP.hpp"
#include "restartsqp/Statistics.hpp"
#include "restartsqp/Types.hpp"
#include <memory>

using namespace Ipopt;

namespace RestartSqp {

/** Special version of the Ipopt AMPL interface, derived from AmplTNLP.
 *  It uses most of the methods from the Ipopt AmplTNLP, except that it
 *  uses different suffixes to communicate bound multipliers and activities.
 */
class IpoptAmplTNlp : public Ipopt::AmplTNLP
{
public:
  /**@name Constructors/Destructors */
  //@{
  /** Constructor. */
  IpoptAmplTNlp(const SmartPtr<const Journalist>& jnlst,
                const SmartPtr<OptionsList> options, char**& argv,
                SmartPtr<AmplSuffixHandler> suffix_handler);

  /** Destructor. */
  ~IpoptAmplTNlp();
  //@}

  /** This is the RestartSqp version for getting the starting point.  We will
   * ignore the one for Ipopt. */
  bool get_starting_point(int num_variables, bool init_primal_variables,
                          double* primal_variables, bool init_bound_multipliers,
                          double* bound_multipliers, int num_constraints,
                          bool init_constraint_multipliers,
                          double* constraint_multipliers);

  /** This is the RestartSqp version for giving the solution to AMPL. */
  void finalize_solution(SqpSolverExitStatus status, int num_variables,
                         const double* primal_solution,
                         const double* bound_multipliers,
                         const ActivityStatus* bound_activity_status,
                         int num_constraints, const double* constraint_values,
                         const double* constraint_multipliers,
                         const ActivityStatus* constraint_activity_status,
                         double objective_value,
                         std::shared_ptr<const Statistics> stats);

  /** Return true if an active set has been given with suffixes. */
  bool use_initial_working_set();

  /** Retrieve the initial working set. */
  bool get_initial_working_sets(int num_variables,
                                ActivityStatus* bounds_working_set,
                                int num_constraints,
                                ActivityStatus* constraints_working_set);

  double get_initial_penalty_parameter();

private:
  /**@name Default Compiler Generated Methods
   * (Hidden to avoid implicit creation/calling).
   *
   * These methods are not implemented and
   * we do not want the compiler to implement
   * them for us, so we declare them private
   * and do not define them. This ensures that
   * they will not be implicitly created/called.
   */
  //@{
  /** Default Constructor */
  IpoptAmplTNlp();

  /** Copy Constructor */
  IpoptAmplTNlp(const IpoptAmplTNlp&);

  /** Default Assignment Operator */
  void operator=(const IpoptAmplTNlp&);
  //@}

  /** @name arrays that store the values of the suffixes.  ASL does not make a
   * copy so we need to keep them around. */
  //@{
  double* asl_primal_solution_;
  double* asl_bound_multipliers_;
  double* asl_constraint_multipliers_;
  int* asl_bounds_activities_;
  int* asl_constraints_activities_;
  double final_penalty_parameter_;
  //@}
};

} // namespace RestartSqp

#endif
