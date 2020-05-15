/* Copyright (C) 2019
* All Rights Reserved.
*
* This code is derived from the file AmplTNLP.cpp in the Ipopt project.
*
* Authors: Xinyi Luo
* Date:2019-07
*/
#include "IpoptAmplTNlp.hpp"

/* AMPL includes */
#include "asl.h"
#include "asl_pfgh.h"
#include "getstub.h"
#include <cassert>

using namespace std;
using namespace Ipopt;

namespace RestartSqp {

/** Default constructor*/
IpoptAmplTNlp::IpoptAmplTNlp(const SmartPtr<const Journalist>& jnlst,
                             const SmartPtr<OptionsList>       options,
                             char**&                           argv,
                             SmartPtr<AmplSuffixHandler>       suffix_handler)
 : AmplTNLP(jnlst, options, argv, suffix_handler)
 , asl_primal_solution_(nullptr)
 , asl_bound_multipliers_(nullptr)
 , asl_constraint_multipliers_(nullptr)
 , asl_bounds_activities_(nullptr)
 , asl_constraints_activities_(nullptr)
{
}

IpoptAmplTNlp::~IpoptAmplTNlp()
{
  delete[] asl_primal_solution_;
  delete[] asl_bound_multipliers_;
  delete[] asl_constraint_multipliers_;
  delete[] asl_bounds_activities_;
  delete[] asl_constraints_activities_;
}

bool IpoptAmplTNlp::get_starting_point(int num_variables, bool init_primal_variables, double* primal_variables,
                                       bool init_bound_multipliers, double* bound_multipliers, int num_constraints,
                                       bool init_constraint_multipliers, double* constraint_multipliers)
{
   ASL_pfgh* asl = asl_;

   // Set the initial point if given, otherwise use the projection of zero into
   // the bounds
   if (init_primal_variables) {
     for (int i = 0; i < num_variables; ++i) {
       primal_variables[i] = havex0[i] ? X0[i] : max(LUv[2 * i], min(LUv[2 * i + 1], 0.0));
      }
   }

  if (init_bound_multipliers) {
    // Modified for warm-start from AMPL
    const double* bound_mult = suffix_handler_->GetNumberSuffixValues("bound_mult", AmplSuffixHandler::Variable_Source);
    if (bound_mult) {
      for (int i = 0; i < num_variables; ++i) {
        bound_multipliers[i] = bound_mult[i];
      }
    }
    else {
      for (int i=0; i<num_variables; ++i) {
        bound_multipliers[i] = 0.;
      }
    }
  }

  if (init_constraint_multipliers) {
    for(int i = 0; i < num_constraints; ++i) {
      // TODO: Is the objective sign correct here?
      constraint_multipliers[i] = havepi0[i] ? obj_sign_ * pi0[i] : 0.0;
    }
  }

  return true;
}

void IpoptAmplTNlp::finalize_solution(SqpSolverExitStatus status, int num_variables,
                 const double* primal_solution,
                         const double* bound_multipliers,
                         const ActivityStatus* bound_activity_status,
                         int num_constraints, const double* constraint_values,
                         const double* constraint_multipliers,
                         const ActivityStatus* constraint_activity_status,
                         double objective_value,
                         std::shared_ptr<const Statistics> stats)
{
  ASL_pfgh* asl = asl_;

  // Determine the objective sign (this is private in the base class)
  int obj_sign = 1;
  if (n_obj > 0 && objtype[obj_no] != 0) {
    obj_sign = -1;
  }

  // Reserve memory for the ASL solution
  if (asl_primal_solution_ == nullptr) {
    assert(!asl_bound_multipliers_);
    assert(!asl_constraint_multipliers_);
    asl_primal_solution_ = new double[num_variables];
    asl_bound_multipliers_ = new double[num_variables];
    asl_constraint_multipliers_ = new double[num_constraints];
  }

  // Copy the primal variables
  for (int i=0; i<num_variables; ++i) {
    asl_primal_solution_[i] = primal_solution[i];
  }

  // Copy the bound multipliers into the array that is given to the ASL, potentially adjusting the sign
  if (obj_sign == 1) {
    for (int i=0; i<num_variables; ++i) {
      asl_bound_multipliers_[i] = bound_multipliers[i];
    }
  } else {
    for (int i=0; i<num_variables; ++i) {
      asl_bound_multipliers_[i] = -bound_multipliers[i];
    }
  }

  // Give the array to the ASL for the bound_mult suffix
  suf_rput("bound_mult", ASL_Sufkind_var, asl_bound_multipliers_);

  // Copy the activities into the ASL array and give them to the ASL, as long as they
  // are available
  if (bound_activity_status) {
    assert(constraint_activity_status);
    if (asl_bounds_activities_ == nullptr) {
      assert(!asl_constraints_activities_);
      asl_bounds_activities_ = new int [num_variables];
      asl_constraints_activities_ = new int [num_constraints];
    }
    for (int i=0; i<num_variables; ++i) {
      asl_bounds_activities_[i] = (int)bound_activity_status[i];
    }
    for (int i=0; i<num_constraints; ++i) {
      asl_constraints_activities_[i] = (int)constraint_activity_status[i];
    }
    // Give the arrays to the ASL for the activity suffix
    suf_iput("activity", ASL_Sufkind_var, asl_bounds_activities_);
    suf_iput("activity", ASL_Sufkind_con, asl_constraints_activities_);
  }
  else {
    delete [] asl_bounds_activities_;
    asl_bounds_activities_ = nullptr;
    delete [] asl_constraints_activities_;
    asl_constraints_activities_ = nullptr;
  }

  // Give the final value of the penatly parameter to ASL
  final_penalty_parameter_ = stats->final_penalty_parameter_;
  suf_rput("penalty_parameter", ASL_Sufkind_obj, &final_penalty_parameter_);

  // Copy the multipliers for the constraints, potentially adjusting the sign
  if (obj_sign == 1) {
    for (int i=0; i<num_constraints; ++i) {
      asl_constraint_multipliers_[i] = constraint_multipliers[i];
    }
  } else {
    for (int i=0; i<num_constraints; ++i) {
      asl_constraint_multipliers_[i] = -constraint_multipliers[i];
    }
  }

  // Determine output message and exit code
  string message;
  if (status == OPTIMAL)
  {
     message = "Locally Optimal Solution Found";
     solve_result_num = 0;
  }
  else if( status == EXCEED_MAX_ITERATIONS )
  {
     message = "Maximum Number of Iterations Exceeded.";
     solve_result_num = 400;
  }
  else if( status == EXCEED_MAX_CPU_TIME )
  {
    message = "Maximum CPU Time Exceeded.";
    solve_result_num = 401;
  }
  else if( status == TRUST_REGION_TOO_SMALL )
  {
    message = "Search Direction becomes Too Small.";
    solve_result_num = 500;
  }
  /*
  else if( status == STOP_AT_ACCEPTABLE_POINT )
  {
    message = "Solved To Acceptable Level.";
    solve_result_num = 1;
  }
  else if( status == FEASIBLE_POINT_FOUND )
  {
    message = "Found feasible point for square problem.";
    solve_result_num = 2;
  }
  else if( status == LOCAL_INFEASIBILITY )
  {
    message = "Converged to a locally infeasible point. Problem may be infeasible.";
    solve_result_num = 200;
  }
  */
  else if( status == CONVERGE_TO_NONOPTIMAL )
  {
    message = "Congergence to nonoptimal Point.";
    solve_result_num = 501;
  }
  else if( status == PRED_REDUCTION_NEGATIVE )
  {
    message = "Predicted reduction is negative.";
    solve_result_num = 502;
  }
  else if( status == PENALTY_TOO_LARGE )
  {
    message = "Penalty parameter becomes too large.";
    solve_result_num = 503;
  }
  /*
  else if( status == DIVERGING_ITERATES )
  {
    message = "Iterates diverging; problem might be unbounded.";
    solve_result_num = 300;
  }
  */
  else
  {
    message = "Unknown Error with status " + to_string(int(status));
    solve_result_num = 510;
  }

  // Write the .sol file
  message = " \nRestartSqp: " + message;
  //write_solution_file(message.c_str());

  // We need to copy the message into a non-const char array to make
  // it work with the AMPL C function.
  char* cmessage = new char[message.length() + 1];
  strcpy(cmessage, message.c_str());

  write_sol(cmessage, asl_primal_solution_, asl_constraint_multipliers_,
            (Option_Info* )Oinfo_ptr_);

  delete[] cmessage;
}

bool IpoptAmplTNlp::use_initial_working_set()
{
  const double* penalty_parameter = suffix_handler_->GetNumberSuffixValues("penalty_parameter", AmplSuffixHandler::Objective_Source);
  printf("act = %x\n", penalty_parameter);
  if (penalty_parameter != nullptr) {
    return true;
  }
  else {
    return false;
  }
}

bool IpoptAmplTNlp::get_initial_working_sets(int num_variables,
                              ActivityStatus* bounds_working_set,
                              int num_constraints,
                              ActivityStatus* constraints_working_set)
{
  // Get the arrays with the activities from ASL
  const int* bounds_activities =
      suffix_handler_->GetIntegerSuffixValues("activity", AmplSuffixHandler::Variable_Source);
  const int* constraints_activities =
      suffix_handler_->GetIntegerSuffixValues("activity", AmplSuffixHandler::Constraint_Source);

  // Transfer results
  if (bounds_activities) {
    for (int i=0; i<num_variables; ++i) {
      bounds_working_set[i] = (ActivityStatus)bounds_activities[i];
    }
  }
  else {
    for (int i=0; i<num_variables; ++i) {
      bounds_working_set[i] = (ActivityStatus)0;
    }
  }
  if (constraints_activities) {
    for (int i=0; i<num_constraints; ++i) {
      constraints_working_set[i] = (ActivityStatus)constraints_activities[i];
    }
  }
  else {
    for (int i=0; i<num_constraints; ++i) {
      constraints_working_set[i] = (ActivityStatus)0;
    }
  }

  return true;
}

double IpoptAmplTNlp::get_initial_penalty_parameter()
{
  const double* asl_penalty_parameter =
      suffix_handler_->GetNumberSuffixValues("penalty_parameter", AmplSuffixHandler::Objective_Source);
  return *asl_penalty_parameter;
}

}

