#include <iostream>
#include <memory>
#include <sqphot/QPhandler.hpp>
#include <sqphot/Algorithm.hpp>
#include <sqphot/Utils.hpp>
#include <IpTNLP.hpp>
#include <AmplTNLP.hpp>
#include <IpIpoptApplication.hpp>
#include <stdio.h>
#include <string.h>


using namespace Ipopt;
using namespace SQPhotstart;

int main(int argc, char** args) {
    Algorithm alg;
    //printf("\n=====================================================================================\n");
    //printf( "	Solving Problem ");
    //std::cout << args[1]<<std::endl;
    //printf("=====================================================================================\n");

        SmartPtr<TNLP> ampl_tnlp = new AmplTNLP(ConstPtr(alg.getJnlst()),
                                                alg.getRoptions2(),
                                                args);
        alg.initialization(ampl_tnlp);
        alg.Optimize();
	if(alg.getExitFlag()==OPTIMAL)
 		std::cout<<args[1]<<std::endl;
	else
		std::cout<<std::endl;



    return 0;






//
//
}
