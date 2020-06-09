// Copyright 2020...
//
// This file is derived from ampl_ipopt.cpp file in the Ipopt project.

#include "SqpAmplTNlp.hpp"
#include "restartsqp/CrossoverSqpSolver.hpp"

using namespace Ipopt;
using namespace std;
using namespace RestartSqp;

int main(
   int argc,
   char** args)
{
#if 0
  // We do not do this yet.

  // Check if executable is run only to print out options documentation
  if( argc == 2 )
  {
    bool print_options = false;
    std::string print_options_mode("text");
    if( !strcmp(args[1], "--print-options=latex") )
    {
      print_options = true;
      print_options_mode = "latex";
    }
    else if( !strcmp(args[1], "--print-options=doxygen") )
    {
      print_options = true;
      print_options_mode = "doxygen";
    }
    else if( !strcmp(args[1], "--print-options") )
    {
      print_options = true;
    }
    else if( !strcmp(args[1], "--print-latex-options") )
    {
      fprintf(stderr, "ampl_ipopt.cpp: Options --print-latex-options has been replaced by --print-options=latex. Please adjust your call.\n");
      exit(-200);
    }
    if( print_options )
    {
      SmartPtr<OptionsList> options = app->Options();
      options->SetStringValue("print_options_documentation", "yes");
      options->SetStringValue("print_options_mode", print_options_mode);
      app->Initialize("");
      return 0;
    }
  }
#endif

  // Create the Crossover SQP algorithm object
  CrossoverSqpSolver sqp_solver;

  // Create a suffix handler for bound multipliers suffix ...
  SmartPtr<AmplSuffixHandler> suffix_handler = new AmplSuffixHandler();
  suffix_handler->AddAvailableSuffix("bound_mult", AmplSuffixHandler::Variable_Source,
                                     AmplSuffixHandler::Number_Type);

  // ... and suffix for the activity status
  //  1: upper bound is active
  // -1: lower bound is active
  // -2: active equality constraints
  //  0: inactice
  suffix_handler->AddAvailableSuffix("activity", AmplSuffixHandler::Variable_Source,
                                     AmplSuffixHandler::Index_Type);
  suffix_handler->AddAvailableSuffix("activity", AmplSuffixHandler::Constraint_Source,
                                     AmplSuffixHandler::Index_Type);

  // Create a suffix for the penalty parameter
  // This is used to carry the parameter from one solve to the next
  suffix_handler->AddAvailableSuffix("penalty_parameter", AmplSuffixHandler::Objective_Source,
                                     AmplSuffixHandler::Number_Type);

  // Create the NLP implementation of the AMPL solver interface, based on the one
  // available from the Ipopt project.  This object has a few specialized methods for
  // SQP-specific tasks
  SmartPtr<IpoptAmplTNlp> ipopt_ampl_tnlp =
      new IpoptAmplTNlp(ConstPtr(sqp_solver.get_jnlst()), sqp_solver.get_ipopt_options_list(), args, suffix_handler);

  // Create the SqpTNlp with the AMPL solver interface
  shared_ptr<SqpAmplTNlp> sqp_ampl_tnlp = make_shared<SqpAmplTNlp>(ipopt_ampl_tnlp);

  // We need to determine whether an initial active has been specified in the model.
  // If not, we call the crossover method, otherwise, we call the regular SQP method
  // (through the CrossoverSqpSolver)

  bool use_initial_working_set = sqp_ampl_tnlp->use_initial_working_set();
  printf("use_initial_working_set = %d\n", use_initial_working_set);

  if (!use_initial_working_set) {
    printf("Calling crossover solver for initial solve.\n");
    sqp_solver.initial_solve(sqp_ampl_tnlp);
  }
  else {
    // Get the last penalty parameter
    double initial_penalty_parameter = sqp_ampl_tnlp->get_initial_penalty_parameter();

    // Set it as an option
    SmartPtr<OptionsList> options = sqp_solver.get_options_list();
    options->SetNumericValue("penalty_parameter_init_value",
                             initial_penalty_parameter);

    printf("Calling crossover solver for resolve.\n");
    sqp_solver.next_solve(sqp_ampl_tnlp);
  }

#if 0
   // Call Initialize the first time to create a journalist, but ignore
   // any options file
   ApplicationReturnStatus retval;
   retval = app->Initialize("");
   if( retval != Solve_Succeeded )
   {
      printf("ampl_ipopt.cpp: Error in first Initialize!!!!\n");
      exit(-100);
   }


   // Call Initialize again to process output related options
   retval = app->Initialize();
   if( retval != Solve_Succeeded )
   {
      printf("ampl_ipopt.cpp: Error in second Initialize!!!!\n");
      exit(-101);
   }

   const int n_loops = 1; // make larger for profiling
   for( Index i = 0; i < n_loops; i++ )
   {
      retval = app->OptimizeTNLP(ampl_tnlp);
   }

   // finalize_solution method in AmplTNLP writes the solution file
#endif

   return 0;
}

