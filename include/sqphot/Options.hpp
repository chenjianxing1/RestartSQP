/* Copyright (C) 2019
* All Rights Reserved.
*
* Authors: Xinyi Luo
* Date:2019-05
*/
#ifndef SQPHOTSTART_OPTIONS_HPP
#define SQPHOTSTART_OPTIONS_HPP

namespace SQPhotstart{
    class Options{
    public:
        /**Default Constructor */
        Options();
        
        /**Destructor**/
        ~Options();
        
        /** Sets all options to the default values */
        int setToDefault();
        
    public:
        /** Public Member Variables*/
        int iter_max;
        int printLevel;
        int qpPrintLevel;
        //int update_method
        bool second_order_correction;
        double eta_c;
        double eta_s;
        double eta_e;
        double gamma_c;
        double gamma_e;
        double delta;
        double delta_max;
        double opt_tol;
        double opt_compl_tol;
        double opt_dual_fea_tol;
        double opt_prim_fea_tol;
        double opt_second_tol;
        double tol;
        double rho;
        int qp_maxiter;
        //double penalty_tol;
        double increase_parm;
        //double decrease_parm;
        //double decrease_tol;
        int rho_max;
        int penalty_iter_max ;
        double eps1;//FIXME: it may need to be changed inside the main loop
        double eps2;
        bool EnablePertubation;
    };//ENDCLASS
}
#endif /* SQPHOTSTART_OPTIONS_HPP */
