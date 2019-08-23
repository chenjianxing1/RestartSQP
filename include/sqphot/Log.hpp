/* Copyright (C) 2019
* All Rights Reserved.
*
* Authors: Xinyi Luo
* Date:2019-06
*/
//This is a simple log_ file for output temporarily
#ifndef SQPHOTSTART_LOG_HPP
#define SQPHOTSTART_LOG_HPP

#include <sqphot/Utils.hpp>

namespace SQPhotstart {
class Log {
public:
    /** Default constructor*/
    Log() {}

    /** Default destructor*/
    ~Log() {}

    void print_header() {
        printf("\n=====================================================================================\n");
        printf("%6s   %23s%9s %9s %9s %9s\n", "iter", "f", "||p_k||", "||c_k||",
               "Delta", "rho");
        printf("=====================================================================================\n");
    }


    void print_main_iter(int iter,
                         double obj_value,
                         double norm_p_k,
                         double infea_measure,
                         double delta,
                         double rho) {
        printf("%6i %23e %9.2e  %9.2e %9.2e %9.2e\n",
               iter, obj_value, norm_p_k, infea_measure, delta, rho);
    }

    void print_penalty_update(int iter, double rho_trial, double infea_measure_model,
                              double infea_measure_infty) {
//            if(iter%10==0){
//            printf("\n=====================================================================================\n");
//            printf("%6s %18s %18s %18s \n", "iter", "rho_trial", "infea_measure_model",
//                    "infea_measure_infty");
//            printf("=====================================================================================\n");
//            }
//            printf("%6i %18e %18e %18e \n",iter,rho_trial,infea_measure_model,
//                    infea_measure_infty);

    }

};

}
#endif //CPP_LOG_HPP

