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

int main(int argc, char** args) {
    Algorithm alg;
    //printf("\n=====================================================================================\n");
    //printf( "	Solving Problem ");
    //std::cout << args[1]<<std::endl;
    //printf("=====================================================================================\n");
    SmartPtr<MyNLP> nlp= new MyNLP();
    try {
        SmartPtr<TNLP> ampl_tnlp = new AmplTNLP(ConstPtr(alg.getJnlst()),
                                                alg.getRoptions2(),
                                                args);
        alg.Optimize(ampl_tnlp);
	if(alg.getExitFlag()==OPTIMAL)
 		std::cout<<args[1]<<std::endl;
	else
		std::cout<<std::endl;

    }
    catch(...) {
        printf( "WARNING, the NLP is invalid!");

    }


    return 0;






//
//
}
