/* Copyright (C) 2005, 2011 International Business Machines and others.
 * All Rights Reserved.
 * This code is published under the Eclipse Public License.
 *
 * Authors:  Carl Laird, Andreas Waechter     IBM    2005-08-17
 */

#include "restartsqp/CInterface.h"
#include <stdlib.h>
#include <assert.h>
#include <stdio.h>

/* Forward Function Declarations */
Bool eval_obj(
    int         num_variables,
    double*     primal_variables,
    Bool        new_primal_variables,
    double*     objective_value,
    UserDataPtr user_data
);

Bool eval_obj_grad(
    int         num_variables,
    double*     primal_variables,
    Bool        new_priml_variables,
    double*     objective_gradient,
    UserDataPtr user_data
);

Bool eval_constr(
    int         num_variables,
    double*     primal_variables,
    Bool        new_primal_variables,
    int         num_constraints,
    double*     constraint_values,
    UserDataPtr user_data
);

Bool eval_constr_jac(
    int         num_variables,
    double*     primal_variables,
    Bool        new_primal_variables,
    int         num_constraints,
    int         num_nonzeros_jacobian,
    int*        row_indices,
    int*        column_indices,
    double*     nonzero_values,
    UserDataPtr user_data
);

Bool eval_lag_hess(
    int         num_variables,
    double*     primal_variables,
    Bool        new_primal_variables,
    double      objective_scaling_factor,
    int         num_constraints,
    double*     constraint_multipliers,
    Bool        new_constraint_multipliers,
    int         num_nonzeros_hessian,
    int*        row_indices,
    int*        column_indices,
    double*     nonzero_values,
    UserDataPtr user_data
);

/** This is an example how user_data can be used. */
struct MyUserData
{
   double g_offset[2]; /**< This is an offset for the constraints.  */
};

/** Main Program */
/* [MAIN] */
int main()
{
   int num_variables = -1;                    /* number of variables */
   int num_constraints = -1;                 /* number of constraints */
   int num_nonzeros_jacobian;                      /* number of nonzeros in the Jacobian of the constraints */
   int num_nonzeros_hessian;                     /* number of nonzeros in the Hessian of the Lagrangian (lower or upper triangular part only) */
   double* variable_lower_bounds = NULL;                  /* lower bounds on x */
   double* variable_upper_bounds = NULL;                  /* upper bounds on x */
   double* constraint_lower_bounds = NULL;                  /* lower bounds on g */
   double* constraint_upper_bounds = NULL;                  /* upper bounds on g */
   SqpProblem sqp_problem = NULL;             /* IpoptProblem */
   enum SqpSolverExitStatus status; /* Solve return code */
   double* primal_variables = NULL;                    /* starting point and solution vector */
   double* bound_multipliers = NULL;               /* constraint multipliers at the solution */
   double* constraint_multipliers = NULL;             /* lower bound multipliers at the solution */
   char* bound_activities = NULL;             /* upper bound multipliers at the solution */
   char* constraint_activities = NULL;             /* upper bound multipliers at the solution */
   double objective_value;                          /* objective value */
   struct MyUserData user_data;         /* our user data for the function evaluations */
   int i;                             /* generic counter */

   /**************************************************************************/
   /*                  Set up the optimization problem                       */
   /**************************************************************************/

   /* set the number of variables and allocate space for the bounds */
   num_variables = 4;
   variable_lower_bounds = (double*) malloc(sizeof(double) * num_variables);
   variable_upper_bounds = (double*) malloc(sizeof(double) * num_variables);
   /* set the values for the variable bounds */
   for( i = 0; i < num_variables; i++ )
   {
      variable_lower_bounds[i] = 1.0;
      variable_upper_bounds[i] = 5.0;
   }

   /* set the number of constraints and allocate space for the bounds */
   num_constraints = 2;
   constraint_lower_bounds = (double*) malloc(sizeof(double) * num_constraints);
   constraint_upper_bounds = (double*) malloc(sizeof(double) * num_constraints);
   /* set the values of the constraint bounds */
   constraint_lower_bounds[0] = 25;
   constraint_upper_bounds[0] = 2e19;
   constraint_lower_bounds[1] = 40;
   constraint_upper_bounds[1] = 40;

   /* set the number of nonzeros in the Jacobian and Hessian */
   num_nonzeros_jacobian = 8;
   num_nonzeros_hessian = 10;


   /* create the IpoptProblem */
   sqp_problem = CreateSqpProblem(num_variables, variable_lower_bounds, variable_upper_bounds, num_constraints,
                                  constraint_lower_bounds, constraint_upper_bounds, num_nonzeros_jacobian, num_nonzeros_hessian,
                            &eval_obj, &eval_obj_grad, &eval_constr,
                            &eval_constr_jac, &eval_lag_hess);

   /* We can free the memory now - the values for the bounds have been
    * copied internally in CreateIpoptProblem
    */
   free(variable_lower_bounds);
   free(variable_upper_bounds);
   free(constraint_lower_bounds);
   free(constraint_upper_bounds);

   /**************************************************************************/
   /*             Solve the problem with the crossover method.               */
   /**************************************************************************/


   /* Set some options.  Note the following ones are only examples,
    * they might not be suitable for your problem.
    */
   AddSqpNumOption(sqp_problem, "tol", 1e-7);
   AddSqpStrOption(sqp_problem, "mu_strategy", "adaptive");
   AddSqpStrOption(sqp_problem, "output_file", "sqp.out");

   /* Here, we do only provide the primal starting point. */
   AddSqpStrOption(sqp_problem, "starting_mode", "primal");

   /* allocate space for the initial point and set the values */
   primal_variables = (double*) malloc(sizeof(double) * num_variables);
   primal_variables[0] = 1.0;
   primal_variables[1] = 5.0;
   primal_variables[2] = 5.0;
   primal_variables[3] = 1.0;

   /* allocate space to store the multipliers and activitities at the solution */
   bound_multipliers = (double*) malloc(sizeof(double) * num_variables);
   constraint_multipliers = (double*) malloc(sizeof(double) * num_constraints);
   bound_activities = (char*) malloc(sizeof(char) * num_variables);
   constraint_activities = (char*) malloc(sizeof(char) * num_constraints);

   /* Initialize the user data */
   user_data.g_offset[0] = 0.;
   user_data.g_offset[1] = 0.;

   /* solve the problem */
   status = CrossoverSolve(sqp_problem, primal_variables, NULL, &objective_value,
                     bound_multipliers, constraint_multipliers, bound_activities,
                     constraint_activities, &user_data);

   if( status == OPTIMAL )
   {
      printf("\n\nSolution of the primal variables, x\n");
      for( i = 0; i < num_variables; i++ )
      {
         printf("x[%d] = %e\n", i, primal_variables[i]);
      }

      printf("\n\nSolution of the bound multipliers, z\n");
      for( i = 0; i < num_variables; i++ )
      {
         printf("z[%d] = %e\n", i, bound_multipliers[i]);
      }
      printf("\n\nSolution of the constraint multipliers, lambda\n");
      for( i = 0; i < num_constraints; i++ )
      {
         printf("lam[%d] = %e\n", i, constraint_multipliers[i]);
      }

      printf("\n\nObjective value\nf(x*) = %e\n", objective_value);
   }
   else
   {
      printf("\n\nERROR OCCURRED DURING CROSSOVER OPTIMIZATION.\n");
   }

   /**************************************************************************/
   /*        Resolve the problem after data and bounds were changed.         */
   /**************************************************************************/

   /* Now we are going to solve this problem again, but with slightly
    * modified constraints.  We change the constraint offset of the
    * first constraint a bit, and resolve the problem using the warm
    * start option.
    */
   user_data.g_offset[0] = 0.2;
   variable_lower_bounds = GetVariableLowerBounds(sqp_problem);
   variable_lower_bounds[0] -= 1.;

   /* We resolve the problem.  This can obly be done if the structure of the
    * optimization problem has not changed (same order of variables and
    * contraints).   This will use the last iterate from the most recent solve
    * and the last workign set as starting point, ignoring any starting point
    * provided by the user.
    */
   status = SqpResolve(sqp_problem, primal_variables, NULL, &objective_value,
                     bound_multipliers, constraint_multipliers, bound_activities,
                     constraint_activities, &user_data);

   if( status == OPTIMAL )
   {
      printf("\n\nSolution of the primal variables, x\n");
      for( i = 0; i < num_variables; i++ )
      {
         printf("x[%d] = %e\n", i, primal_variables[i]);
      }

      printf("\n\nSolution of the bound multipliers, z\n");
      for( i = 0; i < num_variables; i++ )
      {
         printf("z[%d] = %e\n", i, bound_multipliers[i]);
      }
      printf("\n\nSolution of the constraint multipliers, lambda\n");
      for( i = 0; i < num_constraints; i++ )
      {
         printf("lam[%d] = %e\n", i, constraint_multipliers[i]);
      }

      printf("\n\nObjective value\nf(x*) = %e\n", objective_value);
   }
   else
   {
      printf("\n\nERROR OCCURRED DURING CROSSOVER OPTIMIZATION.\n");
   }

#if 0

   if( status == Solve_Succeeded )
   {
      /* Now resolve with a warmstart. */
      AddIpoptStrOption(nlp, "warm_start_init_point", "yes");
      /* The following option reduce the automatic modification of the
       * starting point done my Ipopt.
       */
      AddIpoptNumOption(nlp, "bound_push", 1e-5);
      AddIpoptNumOption(nlp, "bound_frac", 1e-5);
      status = IpoptSolve(nlp, x, NULL, &obj, mult_g, mult_x_L, mult_x_U, &user_data);

      if( status == Solve_Succeeded )
      {
         printf("\n\nSolution of the primal variables, x\n");
         for( i = 0; i < n; i++ )
         {
            printf("x[%d] = %e\n", i, x[i]);
         }

         printf("\n\nSolution of the constraint multipliers, lambda\n");
         for( i = 0; i < m; i++ )
         {
            printf("lambda[%d] = %e\n", i, mult_g[i]);
         }
         printf("\n\nSolution of the bound multipliers, z_L and z_U\n");
         for( i = 0; i < n; i++ )
         {
            printf("z_L[%d] = %e\n", i, mult_x_L[i]);
         }
         for( i = 0; i < n; i++ )
         {
            printf("z_U[%d] = %e\n", i, mult_x_U[i]);
         }

         printf("\n\nObjective value\nf(x*) = %e\n", obj);
      }
      else
      {
         printf("\n\nERROR OCCURRED DURING IPOPT OPTIMIZATION WITH WARM START.\n");
      }
   }
#endif

   /* free allocated memory */
   FreeSqpProblem(sqp_problem);
   free(primal_variables);
   free(bound_multipliers);
   free(constraint_multipliers);
   free(bound_activities);
   free(constraint_activities);

   return (status == OPTIMAL) ? EXIT_SUCCESS : EXIT_FAILURE;
}
/* [MAIN] */

/* Function Implementations */
Bool eval_obj(int num_variables,
   double* primal_variables,
   Bool        new_primal_variables,
   double* objective_value,
   UserDataPtr user_data
)
{
   assert(num_variables == 4);
   (void) num_variables;

   (void) new_primal_variables;
   (void) user_data;

   *objective_value = primal_variables[0] * primal_variables[3] * (primal_variables[0] + primal_variables[1] + primal_variables[2]) + primal_variables[2];

   return TRUE;
}

Bool eval_obj_grad(int num_variables,
   double* primal_variables,
   Bool        new_priml_variables,
   double* objective_gradient,
   UserDataPtr user_data
)
{
   assert(num_variables == 4);
   (void) num_variables;

   (void) new_priml_variables;
   (void) user_data;

   objective_gradient[0] = primal_variables[0] * primal_variables[3] + primal_variables[3] * (primal_variables[0] + primal_variables[1] + primal_variables[2]);
   objective_gradient[1] = primal_variables[0] * primal_variables[3];
   objective_gradient[2] = primal_variables[0] * primal_variables[3] + 1;
   objective_gradient[3] = primal_variables[0] * (primal_variables[0] + primal_variables[1] + primal_variables[2]);

   return TRUE;
}

Bool eval_constr(int num_variables,
   double* primal_variables,
   Bool        new_primal_variables,
   int num_constraints,
   double* constraint_values,
   UserDataPtr user_data
)
{
   struct MyUserData* my_data = user_data;

   assert(num_variables == 4);
   (void) num_variables;
   assert(num_constraints == 2);
   (void) num_constraints;

   (void) new_primal_variables;

   constraint_values[0] = primal_variables[0] * primal_variables[1] * primal_variables[2] * primal_variables[3] + my_data->g_offset[0];
   constraint_values[1] = primal_variables[0] * primal_variables[0] + primal_variables[1] * primal_variables[1] + primal_variables[2] * primal_variables[2] + primal_variables[3] * primal_variables[3] + my_data->g_offset[1];

   return TRUE;
}

Bool eval_constr_jac(int num_variables,
   double* primal_variables,
   Bool        new_primal_variables,
   int num_constraints,
   int num_nonzeros_jacobian,
   int* row_indices,
   int* column_indices,
   double* nonzero_values,
   UserDataPtr user_data
)
{
   (void) num_variables;
   (void) new_primal_variables;
   (void) num_constraints;
   (void) num_nonzeros_jacobian;
   (void) user_data;

   if( nonzero_values == NULL )
   {
      /* return the structure of the jacobian */

      /* this particular jacobian is dense */
      row_indices[0] = 0;
      column_indices[0] = 0;
      row_indices[1] = 0;
      column_indices[1] = 1;
      row_indices[2] = 0;
      column_indices[2] = 2;
      row_indices[3] = 0;
      column_indices[3] = 3;
      row_indices[4] = 1;
      column_indices[4] = 0;
      row_indices[5] = 1;
      column_indices[5] = 1;
      row_indices[6] = 1;
      column_indices[6] = 2;
      row_indices[7] = 1;
      column_indices[7] = 3;
   }
   else
   {
      /* return the values of the jacobian of the constraints */

      nonzero_values[0] = primal_variables[1] * primal_variables[2] * primal_variables[3]; /* 0,0 */
      nonzero_values[1] = primal_variables[0] * primal_variables[2] * primal_variables[3]; /* 0,1 */
      nonzero_values[2] = primal_variables[0] * primal_variables[1] * primal_variables[3]; /* 0,2 */
      nonzero_values[3] = primal_variables[0] * primal_variables[1] * primal_variables[2]; /* 0,3 */

      nonzero_values[4] = 2 * primal_variables[0]; /* 1,0 */
      nonzero_values[5] = 2 * primal_variables[1]; /* 1,1 */
      nonzero_values[6] = 2 * primal_variables[2]; /* 1,2 */
      nonzero_values[7] = 2 * primal_variables[3]; /* 1,3 */
   }

   return TRUE;
}

Bool eval_lag_hess(int num_variables,
   double* primal_variables,
   Bool        new_primal_variables,
   double objective_scaling_factor,
   int num_constraints,
   double* constraint_multipliers,
   Bool        new_constraint_multipliers,
   int num_nonzeros_hessian,
   int* row_indices,
   int* column_indices,
   double* nonzero_values,
   UserDataPtr user_data
)
{
   int idx = 0; /* nonzero element counter */
   int row = 0; /* row counter for loop */
   int col = 0; /* col counter for loop */

   (void) num_variables;
   (void) new_primal_variables;
   (void) num_constraints;
   (void) new_constraint_multipliers;
   (void) user_data;

   if( nonzero_values == NULL )
   {
      /* return the structure. This is a symmetric matrix, fill the lower left
       * triangle only. */

      /* the hessian for this problem is actually dense */
      idx = 0;
      for( row = 0; row < 4; row++ )
      {
         for( col = 0; col <= row; col++ )
         {
            row_indices[idx] = row;
            column_indices[idx] = col;
            idx++;
         }
      }

      assert(idx == num_nonzeros_hessian);
      (void) num_nonzeros_hessian;
   }
   else
   {
      /* return the values. This is a symmetric matrix, fill the lower left
       * triangle only */

      /* fill the objective portion */
      nonzero_values[0] = objective_scaling_factor * (2 * primal_variables[3]); /* 0,0 */

      nonzero_values[1] = objective_scaling_factor * (primal_variables[3]); /* 1,0 */
      nonzero_values[2] = 0; /* 1,1 */

      nonzero_values[3] = objective_scaling_factor * (primal_variables[3]); /* 2,0 */
      nonzero_values[4] = 0; /* 2,1 */
      nonzero_values[5] = 0; /* 2,2 */

      nonzero_values[6] = objective_scaling_factor * (2 * primal_variables[0] + primal_variables[1] + primal_variables[2]); /* 3,0 */
      nonzero_values[7] = objective_scaling_factor * (primal_variables[0]); /* 3,1 */
      nonzero_values[8] = objective_scaling_factor * (primal_variables[0]); /* 3,2 */
      nonzero_values[9] = 0; /* 3,3 */

      /* add the portion for the first constraint */
      nonzero_values[1] -= constraint_multipliers[0] * (primal_variables[2] * primal_variables[3]); /* 1,0 */

      nonzero_values[3] -= constraint_multipliers[0] * (primal_variables[1] * primal_variables[3]); /* 2,0 */
      nonzero_values[4] -= constraint_multipliers[0] * (primal_variables[0] * primal_variables[3]); /* 2,1 */

      nonzero_values[6] -= constraint_multipliers[0] * (primal_variables[1] * primal_variables[2]); /* 3,0 */
      nonzero_values[7] -= constraint_multipliers[0] * (primal_variables[0] * primal_variables[2]); /* 3,1 */
      nonzero_values[8] -= constraint_multipliers[0] * (primal_variables[0] * primal_variables[1]); /* 3,2 */

      /* add the portion for the second constraint */
      nonzero_values[0] -= constraint_multipliers[1] * 2; /* 0,0 */

      nonzero_values[2] -= constraint_multipliers[1] * 2; /* 1,1 */

      nonzero_values[5] -= constraint_multipliers[1] * 2; /* 2,2 */

      nonzero_values[9] -= constraint_multipliers[1] * 2; /* 3,3 */
   }

   return TRUE;
}
