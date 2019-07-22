#include <iostream>
#include <memory>
#include <sqphot/QPhandler.hpp>
#include <sqphot/MyNLP.hpp>
#include <sqphot/Algorithm.hpp>
#include <sqphot/Utils.hpp>
#include <coin/IpTNLP.hpp>
#include <coin/AmplTNLP.hpp>
#include <coin/IpIpoptApplication.hpp>
#include <stdio.h>
#include <string.h>


using namespace Ipopt;
using namespace SQPhotstart;

int main(int argc, char** args){



        SmartPtr<IpoptApplication> app = IpoptApplicationFactory();


//    app->RethrowNonIpoptException(false);

    // Check if executable is run only to print out options documentation

    // Call Initialize the first time to create a journalist, but ignore
    // any options file
    ApplicationReturnStatus retval;
//    retval = app->Initialize("");
//    if (retval != Solve_Succeeded) {
//        printf("ampl_ipopt.cpp: Error in first Initialize!!!!\n");
//        exit(-100);
//    }

    Algorithm alg;

    SmartPtr<TNLP> ampl_tnlp = new AmplTNLP(ConstPtr(alg.getJnlst()),
                                            alg.getRoptions2(),
                                            args);

    alg.Optimize(ampl_tnlp);
    // Call Initialize again to process output related options

    // finalize_solution method in AmplTNLP writes the solution file

    return 0;





//    SmartPtr<MyNLP> nlp= new MyNLP();
//
//
}
