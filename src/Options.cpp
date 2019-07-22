/* Copyright (C) 2019
 * All Rights Reserved.
 *
 * Authors: Xinyi Luo
 * Date:2019-07
 */
#include <sqphot/Options.hpp>

namespace SQPhotstart{
    
    Options::Options(){
        setToDefault();
    }
    
    Options::~Options(){
    }
    
    
    int Options::setToDefault()
    {
        iter_max = 200;
        printLevel = 2;
        qpPrintLevel = 0;       //does not print anything
        QPsolverChoice = "qpOASES";
        LPsolverChoice = "qpOASES";
        second_order_correction = true;
        penalty_update = true;
        eta_c = 0.25;
        eta_s = 1.0e-8;
        eta_e = 0.75;
        gamma_c = 0.5;
        gamma_e = 2;
        delta = 1;
        delta_max =1.0e8;
        opt_tol = 1.0e-5;
        opt_compl_tol = 1.0e-6;
        opt_dual_fea_tol = 1.0e-6;
        opt_prim_fea_tol = 1.0e-5;
        opt_second_tol = 1.0e-8;
        tol = 1.0e-8;
        penalty_update_tol = 1.0e-8;
        rho = 1;
        qp_maxiter = 1000;
        //penalty_tol = 1.0e-8;
        increase_parm = 10.0;
        rho_max = 1.0e6;
        penalty_iter_max = 1000;
        eps1 = 0.3;
        eps2 = 1.0e-6;
        EnablePertubation = false;
        lp_maxiter = 100;
        return 0;
        
    }
}

