#include <iostream>
#include <memory>
#include <sqphot/QPhandler.hpp>
#include <sqphot/MyNLP.hpp>
#include <sqphot/Algorithm.hpp>
#include <sqphot/Utils.hpp>
#include <coin/IpTNLP.hpp>


using namespace SQPhotstart;
int main(){
   // USING_NAMESPACE_QPOASES

    SmartPtr<MyNLP> nlp= new MyNLP();
    Algorithm alg;
    alg.Optimize(nlp);
}
