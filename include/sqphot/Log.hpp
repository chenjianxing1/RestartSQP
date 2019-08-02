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
        printf("%6s   %23s%9s %9s %9s %9s\n","iter", "f", "||p_k||","||c_k||", "Delta","rho");
        printf("=====================================================================================\n");
    }



    void print_main_iter(int iter,
                         double obj_value,
                         double norm_p_k,
                         double infea_measure,
                         double delta,
                         double rho) {
        printf("%6i %23e %9.2e  %9.2e %9.2e %9.2e\n",
               iter, obj_value, norm_p_k, infea_measure,delta, rho);
    }

    void print_final(int iter,
                     int qp_iter,
                     double obj_value,
                     double norm_p_k,
                     double infea_measure,
                     int exitflag) {
        printf("\n=====================================================================================\n");
        switch(exitflag) {
        case OPTIMAL:
            printf("Exitflag:                                                   %23s\n","OPTIMAL");
            break;
        case INVALID_NLP :
            printf("Exitflag:                                                   %23s\n","INVALID_NLP");
            break;
        case EXCEED_MAX_ITER :
            printf("Exitflag:                                                   %23s\n","EXCEED_MAX_ITER");
            break;
        case QPERROR_INTERNAL_ERROR :
            printf("Exitflag:                                                   %23s\n","QP_INTERNAL_ERROR");
            break;
        case QPERROR_INFEASIBLE :
            printf("Exitflag:                                                   %23s\n","QP_INFEASIBLE");
            break;
        case QPERROR_UNBOUNDED :
            printf("Exitflag:                                                   %23s\n","QP_UNBOUNDED");
            break;
        case QPERROR_EXCEED_MAX_ITER :
            printf("Exitflag:                                                   %23s\n","QP_EXCEED_MAX_ITER");
            break;
        case QPERROR_NOTINITIALISED :
            printf("Exitflag:                                                   %23s\n","QP_NOTINITIALISED");
            break;
        case AUXINPUT_NOT_OPTIMAL :
            printf("Exitflag:                                                   %23s\n","AUXINPUT_NOT_OPTIMAL");
            break;
        case CONVERGE_TO_NONOPTIMAL :
            printf("Exitflag:                                                   %23s\n","CONVERGE_TO_NONOPTIMAL");
            break;
        case UNKNOWN :
            printf("Exitflag:                                                   %23s\n","UNKNOWN ERROR");

            break;
        }

        printf("Iterations:                                                 %23i\n", iter);
        printf("QP Solver Iterations:                                       %23i\n", qp_iter);
        printf("Final Objectives:                                           %23e\n",obj_value);
        printf("||p_k||                                                     %23e\n",norm_p_k);
        printf("||c_k||                                                     %23e\n",infea_measure);

        printf("=====================================================================================\n\n");

    }
};

}
#endif //CPP_LOG_HPP

