//This is a simple log file for output temporarily
#ifndef SQPHOTSTART_LOG_HPP
#define SQPHOTSTART_LOG_HPP

#include <sqphot/Utils.hpp>
namespace SQPhotstart{
    class Log{
    public:
        /** Default constructor*/
        Log(){}
        /** Default destructor*/
        ~Log(){}

        void print_header(){
            printf("\n=====================================================================================\n");
            printf("%6s %23s %9s %9s %9s %9s\n","iter", "f", "||p_k||","||c_k||", "Delta","rho");
            printf("=====================================================================================\n");
        }

        void print_main_iter(int iter,
                             double obj_value,
                             double norm_p_k,
                             double infea_measure,
                             double delta,
                             double rho){
            printf("%6i %23e %9.2e  %9.2e %9.2e %9.2e\n",
          iter, obj_value, norm_p_k, infea_measure,delta, rho);
        }

        void print_final(int iter,
                          int qp_iter,
                          double obj_value,
                          double norm_p_k,
                          double infea_measure,
                          int exitflag){
            printf("\n=====================================================================================\n");
            printf("Exitflag:                                                   %23i\n", exitflag);
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
