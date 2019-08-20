/* Copyright (C) 2019
 * All Rights Reserved.
 *
 * Authors: Xinyi Luo
 * Date:2019-07
 */
#include <sqphot/Options.hpp>

namespace SQPhotstart {

Options::Options() {
    setToDefault();
}

Options::~Options() {
}


int Options::setToDefault() {
    iter_max =50;
    printLevel = 2;
    qpPrintLevel = 0;       //does not print anything
    QPsolverChoice = QORE_QP;
    LPsolverChoice = QORE_LP;
    second_order_correction = false;
    penalty_update = false;
    eta_c = 0.25;
    eta_s = 1.0e-10;
    eta_e = 0.75;
    gamma_c = 0.5;
    gamma_e = 2;
    delta = 1;
    delta_min = 1.0e-16;
    delta_max = 1.0e10;
    active_set_tol = 1.0e-5;
    opt_tol = 1.0e-3;
    opt_compl_tol = 1.0e-4;
    opt_dual_fea_tol = 1.0e-5;
    opt_prim_fea_tol = 1.0e-5;
    opt_second_tol = 1.0e-8;
    tol = 1.0e-8;
    penalty_update_tol = 1.0e-8;
    rho = 100;
    qp_maxiter = 1000;
    //penalty_tol = 1.0e-8;
    increase_parm = 2;
    rho_max = 1.0e8;
    penalty_iter_max = 200;
    eps1 = 0.3;
    eps1_change_parm = 0.1;
    eps2 = 1.0e-6;
    EnablePertubation = false;
    lp_maxiter = 100;
    return 0;

}
}


