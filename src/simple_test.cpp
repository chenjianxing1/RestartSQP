#include <iostream>
#include <memory>
#include "QPhandler.hpp"
#include "MyNLP.hpp"
#include "Algorithm.hpp"
#include "Utils.hpp"
#include "IpTNLP.hpp"


using namespace SQPhotstart;
int main(){
   // USING_NAMESPACE_QPOASES

    SmartPtr<MyNLP> nlp= new MyNLP();
    Algorithm alg;
    alg.Optimize(nlp);
}
