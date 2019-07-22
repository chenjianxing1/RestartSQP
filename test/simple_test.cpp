//#include <iostream>
//#include <memory>
//#include <sqphot/QPhandler.hpp>
//#include <sqphot/MyNLP.hpp>
//#include <sqphot/Algorithm.hpp>
//#include <sqphot/Utils.hpp>
#include <coin/IpTNLP.hpp>
#include <coin/AmplTNLP.hpp>
#include <coin/IpIpoptApplication.hpp>
#include <stdio.h>
#include <string.h>


using namespace Ipopt;
//using namespace SQPhotstart;

int main(int argc, char** args){
    using namespace Ipopt;

    SmartPtr<IpoptApplication> app = IpoptApplicationFactory();
    app->RethrowNonIpoptException(false);

    // Check if executable is run only to print out options documentation
    if (argc == 2) {
        bool print_options = false;
        bool print_latex_options = false;
        if (!strcmp(args[1],"--print-options")) {
            print_options = true;
        }
        else if (!strcmp(args[1],"--print-latex-options")) {
            print_options = true;
            print_latex_options = true;
        }
        if (print_options) {
            SmartPtr<OptionsList> options = app->Options();
            options->SetStringValue("print_options_documentation", "yes");
            if (print_latex_options) {
                options->SetStringValue("print_options_latex_mode", "yes");
            }
            app->Initialize("");
            return 0;
        }
    }

    // Call Initialize the first time to create a journalist, but ignore
    // any options file
    ApplicationReturnStatus retval;
    retval = app->Initialize("");
    if (retval != Solve_Succeeded) {
        printf("ampl_ipopt.cpp: Error in first Initialize!!!!\n");
        exit(-100);
    }

    // Add the suffix handler for scaling
    SmartPtr<AmplSuffixHandler> suffix_handler = new AmplSuffixHandler();
    suffix_handler->AddAvailableSuffix("scaling_factor", AmplSuffixHandler::Variable_Source, AmplSuffixHandler::Number_Type);
    suffix_handler->AddAvailableSuffix("scaling_factor", AmplSuffixHandler::Constraint_Source, AmplSuffixHandler::Number_Type);
    suffix_handler->AddAvailableSuffix("scaling_factor", AmplSuffixHandler::Objective_Source, AmplSuffixHandler::Number_Type);
    // Modified for warm-start from AMPL
    suffix_handler->AddAvailableSuffix("ipopt_zL_out", AmplSuffixHandler::Variable_Source, AmplSuffixHandler::Number_Type);
    suffix_handler->AddAvailableSuffix("ipopt_zU_out", AmplSuffixHandler::Variable_Source, AmplSuffixHandler::Number_Type);
    suffix_handler->AddAvailableSuffix("ipopt_zL_in", AmplSuffixHandler::Variable_Source, AmplSuffixHandler::Number_Type);
    suffix_handler->AddAvailableSuffix("ipopt_zU_in", AmplSuffixHandler::Variable_Source, AmplSuffixHandler::Number_Type);

    SmartPtr<TNLP> ampl_tnlp = new AmplTNLP(ConstPtr(app->Jnlst()),
                                            app->Options(),
                                            args, suffix_handler);

    // Call Initialize again to process output related options
    retval = app->Initialize();
    if (retval != Solve_Succeeded) {
        printf("ampl_ipopt.cpp: Error in second Initialize!!!!\n");
        exit(-101);
    }

    const int n_loops = 1; // make larger for profiling
    for (Index i=0; i<n_loops; i++) {
        retval = app->OptimizeTNLP(ampl_tnlp);
    }

    // finalize_solution method in AmplTNLP writes the solution file

    return 0;


    // finalize_solution method in AmplTNLP writes the solution file


//
//    SmartPtr<MyNLP> nlp= new MyNLP();
//
//    Algorithm alg;
//    alg.Optimize(nlp);

}
