#ifndef _SQPHOTSTART_MESSAGE_HANDLE_
#define _SQPHOTSTART_MESSAGE_HANDLE_

//OUTPUTS
#define DOUBLE_DIVIDER "===============================================================\n"
#define SINGLE_DIVIDER "---------------------------------------------------------------\n"
#define DOUBLE_LONG_DIVIDER "=====================================================================================\n"

#define STANDARD_HEADER "%6s  %23s%9.2s %9.2s %9.2s %9.2s\n", "iter", "f", "||p_k||", "||c_k||", "Delta","rho"
#define STANDARD_OUTPUT "%6i  %23e%9.2e %9.2e %9.2e %9.2e\n", stats_->iter,obj_value_,norm_p_k_, infea_measure_, delta_, rho_
#define PENALTY_UPDATE_HEADER "%6s %18s %18s %18s \n", "iter", "rho_trial", "infea_measure_model" "infea_measure_inf"
#define PENALTY_UPDATE_OUTPUT "%6i %18e %18e %18e \n",iter,rho_trial,infea_measure_model,infea_measure_infty
//WARNING MESSAGES
#define INVALID_WORKING_SET_MSG "Warning! The working set index in QP is not well defined!"
#define INVALID_RETURN_TYPE_MSG " The return type is invalid for QOREInterface"
#endif
