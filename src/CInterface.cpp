// Copyright (C) 2004, 2011 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

/// @todo use memcpy for various for-loop-copies

#include "restartsqp/CInterface.h"
#include "restartsqp/CInterfaceTNlp.hpp"
#include "restartsqp/CrossoverSqpSolver.hpp"

using namespace std;
using namespace RestartSqp;
using namespace Ipopt;

struct SqpProblemInfo
{
   shared_ptr<CrossoverSqpSolver> sqp_solver;
   int           num_variables;
   double*         variable_lower_bounds;
   double*         variable_upper_bounds;
   int           num_constraints;
   double*         constraint_lower_bounds;
   double*         constraint_upper_bounds;
   int           num_nonzeros_jacobian;
   int           num_nonzeros_hessian;
   Eval_Obj_CB       eval_obj;
   Eval_Obj_Grad_CB  eval_obj_grad;
   Eval_Constr_CB       eval_constr;
   Eval_Constr_Jac_CB   eval_constr_jac;
   Eval_Lag_Hess_CB       eval_lag_hess;
   string        output_file_name;
};

SqpProblem CreateSqpProblem(int          num_variables,
   double*        variable_lower_bounds,
   double*        variable_upper_bounds,
   int          num_constraints,
   double*        constraint_lower_bounds,
   double*        constraint_upper_bounds,
   int          num_nonzeros_jacobian,
   int          num_nonzeros_hessian,
   Eval_Obj_CB      eval_obj,
   Eval_Obj_Grad_CB eval_obj_grad,
   Eval_Constr_CB      eval_constr,
   Eval_Constr_Jac_CB  eval_constr_jac,
   Eval_Lag_Hess_CB     eval_lag_hess
)
{
#if 0
  NEED TO UPDATE
   // make sure input is Ok
   if( n < 1 || m < 0 || !x_L || !x_U || (m > 0 && (!g_L || !g_U)) || (m == 0 && nele_jac != 0)
       || (m > 0 && nele_jac < 1) || nele_hess < 0 || !eval_f || !eval_grad_f || (m > 0 && (!eval_g || !eval_jac_g)) )
   {
      return NULL;
   }
#endif

   SqpProblem retval = new SqpProblemInfo;

   retval->num_variables   = num_variables;
   retval->variable_lower_bounds = new double[num_variables];
   for( int i = 0; i < num_variables; i++ )
   {
      retval->variable_lower_bounds[i] = variable_lower_bounds[i];
   }

   retval->variable_upper_bounds = new double[num_variables];
   for( int i = 0; i < num_variables; i++ )
   {
      retval->variable_upper_bounds[i] = variable_upper_bounds[i];
   }

   retval->num_constraints = num_constraints;
   if( num_constraints > 0 )
   {
      retval->constraint_lower_bounds = new double[num_constraints];
      for( int i = 0; i < num_constraints; i++ )
      {
         retval->constraint_lower_bounds[i] = constraint_lower_bounds[i];
      }

      retval->constraint_upper_bounds = new double[num_constraints];
      for( int i = 0; i < num_constraints; i++ )
      {
         retval->constraint_upper_bounds[i] = constraint_upper_bounds[i];
      }
   }
   else
   {
      retval->constraint_lower_bounds = NULL;
      retval->constraint_upper_bounds = NULL;
   }

   retval->sqp_solver = make_shared<CrossoverSqpSolver>();

   retval->num_nonzeros_jacobian = num_nonzeros_jacobian;
   retval->num_nonzeros_hessian = num_nonzeros_hessian;
   retval->eval_obj = eval_obj;
   retval->eval_obj_grad = eval_obj_grad;
   retval->eval_constr = eval_constr;
   retval->eval_constr_jac = eval_constr_jac;
   retval->eval_lag_hess = eval_lag_hess;

   return retval;
}

void FreeSqpProblem(
   SqpProblem sqp_problem
)
{
   sqp_problem->sqp_solver.reset();

   delete[] sqp_problem->variable_lower_bounds;
   delete[] sqp_problem->variable_upper_bounds;
   delete[] sqp_problem->constraint_lower_bounds;
   delete[] sqp_problem->constraint_upper_bounds;

   delete sqp_problem;
}

Bool AddSqpStrOption(SqpProblem sqp_problem,
   char*        keyword,
   char*        val
)
{
   return (Bool) sqp_problem->sqp_solver->get_options_list()->SetStringValue(keyword, val);
}

Bool AddSqpNumOption(
   SqpProblem sqp_problem,
   char*        keyword,
   double       val
)
{
   return (Bool) sqp_problem->sqp_solver->get_options_list()->SetNumericValue(keyword, val);
}

Bool AddSqpIntOption(
   SqpProblem sqp_problem,
   char*        keyword,
   int          val
)
{
   return (Bool) sqp_problem->sqp_solver->get_options_list()->SetIntegerValue(keyword, val);
}

Bool SetOutputFileName(
   SqpProblem sqp_problem,
   char*        file_name
)
{
   sqp_problem->output_file_name = file_name;
   return (Bool) true;
}

double* GetVariableLowerBounds(SqpProblem sqp_problem)
{
  return sqp_problem->variable_lower_bounds;
}

double* GetVariableUpperBounds(SqpProblem sqp_problem)
{
  return sqp_problem->variable_upper_bounds;
}

double* GetConstraintLowerBounds(SqpProblem sqp_problem)
{
  return sqp_problem->constraint_lower_bounds;
}

double* GetConstraintUpperBounds(SqpProblem sqp_problem)
{
  return sqp_problem->constraint_upper_bounds;
}

void allocate_initial_arrays_(
    SqpProblem sqp_problem,
    double*      primal_variables,
    double*      constraint_values,
    double*      objective_value,
    double*      bound_multipliers,
    double*      constraint_multipliers,
    char*   bound_activities,
    char* constraint_activities,
    UserDataPtr  user_data,
    double*&      initial_primal_variables,
    double*&      initial_bound_multipliers,
    double*&      initial_constraint_multipliers,
    char*&   initial_bound_activities,
    char*&   initial_constraint_activities,
    shared_ptr<CInterfaceTNlp>& sqp_tnlp
    )
{
  int num_variables = sqp_problem->num_variables;
  int num_constraints = sqp_problem->num_constraints;

   initial_primal_variables = new double[num_variables];
   std::copy(primal_variables, primal_variables+num_variables,
             initial_primal_variables);

   initial_bound_multipliers = NULL;
   if (bound_multipliers) {
     initial_bound_multipliers = new double[num_variables];
     std::copy(bound_multipliers, bound_multipliers+num_variables,
               initial_bound_multipliers);
   }

   initial_constraint_multipliers = NULL;
   if (constraint_multipliers) {
     initial_constraint_multipliers = new double[num_constraints];
     std::copy(constraint_multipliers, constraint_multipliers+num_constraints,
               initial_constraint_multipliers);
   }

   initial_bound_activities = NULL;
   if (bound_activities) {
     initial_bound_activities = new char[num_variables];
     std::copy(bound_activities, bound_activities+num_variables,
               initial_bound_activities);
   }

   initial_constraint_activities = NULL;
   if (constraint_activities) {
     initial_constraint_activities = new char[num_constraints];
     std::copy(constraint_activities, constraint_activities+num_constraints,
               initial_constraint_activities);
   }

   // Create the Cinterface NLP
   sqp_tnlp =
       make_shared<CInterfaceTNlp>(num_variables, sqp_problem->variable_lower_bounds,sqp_problem->variable_upper_bounds,
                                         num_constraints, sqp_problem->constraint_lower_bounds, sqp_problem->constraint_upper_bounds,
                                         sqp_problem->num_nonzeros_jacobian, sqp_problem->num_nonzeros_hessian,
                                         initial_primal_variables, initial_bound_multipliers, initial_constraint_multipliers,
                                         initial_bound_activities, initial_constraint_activities,
                                         sqp_problem->eval_obj, sqp_problem->eval_obj_grad,
                                         sqp_problem->eval_constr, sqp_problem->eval_constr_jac, sqp_problem->eval_lag_hess,
                                         primal_variables, bound_multipliers, constraint_multipliers, constraint_values,
                                         objective_value, bound_activities, constraint_activities, user_data);


}

SqpSolverExitStatus SqpSolve(
   SqpProblem sqp_problem,
   double*      primal_variables,
   double*      constraint_values,
   double*      objective_value,
   double*      bound_multipliers,
   double*      constraint_multipliers,
   char*   bound_activities,
   char* constraint_activities,
   UserDataPtr  user_data
)
{
   if( !primal_variables )
   {
      sqp_problem->sqp_solver->get_jnlst()->Printf(J_ERROR, J_MAIN, "Error: Array x with starting point information is NULL.");
      return SqpSolverExitStatus(INVALID_NLP);
   }

   double* initial_primal_variables;
   double* initial_bound_multipliers;
   double* initial_constraint_multipliers;
   char* initial_bound_activities;
   char* initial_constraint_activities;
   shared_ptr<CInterfaceTNlp> sqp_tnlp;

   allocate_initial_arrays_(sqp_problem,
                            primal_variables, constraint_values, objective_value,
                            bound_multipliers, constraint_multipliers,
                            bound_activities, constraint_activities, user_data,
                            initial_primal_variables, initial_bound_multipliers,
                            initial_constraint_multipliers, initial_bound_activities,
                            initial_constraint_activities,
                            sqp_tnlp);

   sqp_problem->sqp_solver->solve(sqp_tnlp);

   delete[] initial_primal_variables;
   delete[] initial_bound_multipliers;
   delete[] initial_constraint_multipliers;
   delete[] initial_bound_activities;
   delete[] initial_constraint_activities;

   SqpSolverExitStatus retval = sqp_problem->sqp_solver->get_exit_flag();
   printf("retval = %d\n", (int)retval);

   return retval;
}


SqpSolverExitStatus CrossoverSolve(
   SqpProblem sqp_problem,
   double*      primal_variables,
   double*      constraint_values,
   double*      objective_value,
   double*      bound_multipliers,
   double*      constraint_multipliers,
   char*   bound_activities,
   char* constraint_activities,
   UserDataPtr  user_data
)
{
   if( !primal_variables )
   {
      sqp_problem->sqp_solver->get_jnlst()->Printf(J_ERROR, J_MAIN, "Error: Array x with starting point information is NULL.");
      return SqpSolverExitStatus(INVALID_NLP);
   }

   double* initial_primal_variables;
   double* initial_bound_multipliers;
   double* initial_constraint_multipliers;
   char* initial_bound_activities;
   char* initial_constraint_activities;
   shared_ptr<CInterfaceTNlp> sqp_tnlp;

   allocate_initial_arrays_(sqp_problem,
                            primal_variables, constraint_values, objective_value,
                            bound_multipliers, constraint_multipliers,
                            bound_activities, constraint_activities, user_data,
                            initial_primal_variables, initial_bound_multipliers,
                            initial_constraint_multipliers, initial_bound_activities,
                            initial_constraint_activities,
                            sqp_tnlp);

   sqp_problem->sqp_solver->crossover_solve(sqp_tnlp);

   delete[] initial_primal_variables;
   delete[] initial_bound_multipliers;
   delete[] initial_constraint_multipliers;
   delete[] initial_bound_activities;
   delete[] initial_constraint_activities;

   SqpSolverExitStatus retval = sqp_problem->sqp_solver->get_exit_flag();
   printf("retval = %d\n", (int)retval);

   return retval;
}

SqpSolverExitStatus SqpResolve(
   SqpProblem sqp_problem,
   double*      primal_variables,
   double*      constraint_values,
   double*      objective_value,
   double*      bound_multipliers,
   double*      constraint_multipliers,
   char*   bound_activities,
   char* constraint_activities,
   UserDataPtr  user_data
)
{
   double* initial_primal_variables;
   double* initial_bound_multipliers;
   double* initial_constraint_multipliers;
   char* initial_bound_activities;
   char* initial_constraint_activities;
   shared_ptr<CInterfaceTNlp> sqp_tnlp;

   allocate_initial_arrays_(sqp_problem,
                            primal_variables, constraint_values, objective_value,
                            bound_multipliers, constraint_multipliers,
                            bound_activities, constraint_activities, user_data,
                            initial_primal_variables, initial_bound_multipliers,
                            initial_constraint_multipliers, initial_bound_activities,
                            initial_constraint_activities,
                            sqp_tnlp);

   sqp_problem->sqp_solver->resolve(sqp_tnlp);

   delete[] initial_primal_variables;
   delete[] initial_bound_multipliers;
   delete[] initial_constraint_multipliers;
   delete[] initial_bound_activities;
   delete[] initial_constraint_activities;

   SqpSolverExitStatus retval = sqp_problem->sqp_solver->get_exit_flag();
   printf("retval = %d\n", (int)retval);

   return retval;
}
