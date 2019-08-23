#ifndef _SQPHOTSTART_MESSAGE_HANDLE_
#define _SQPHOTSTART_MESSAGE_HANDLE_

#define DOUBLE_DIVIDER "===============================================\n"
#define SINGLE_DIVIDER "-----------------------------------------------\n"
#define DOUBLE_LONG_DIVIDER "=====================================================================================\n"
#define STANDARD_HEADER "%6s   %23s%9s %9s %9s %9s\n", "iter", "f", "||p_k||", "||c_k||", "Delta","rho"

#define STANDARD_OUTPUT "%6i %23e %9.2e  %9.2e %9.2e %9.2e\n", stats_->iter,obj_value_,norm_p_k_, infea_measure_, delta_, rho_
#define PENALTY_UPDATE_HEADER "%6s %18s %18s %18s \n", "iter", "rho_trial", "infea_measure_model"
#define PENALTY_UPDATE_OUTPUT "%6i %18e %18e %18e \n",iter,rho_trial,infea_measure_model,infea_measure_infty
#endif
