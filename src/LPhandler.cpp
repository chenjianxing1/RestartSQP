/* Copyright (C) 2019
 * All Rights Reserved.
 *
 * Authors: Xinyi Luo
 * Date:2019-07
 */

#include <sqphot/LPhandler.hpp>

namespace SQPhotstart {

LPhandler::LPhandler(Index_info nlp_info, shared_ptr<const Options> options,
                     Ipopt::SmartPtr<Ipopt::Journalist> jnlst) :
    QPhandler(nlp_info, options, jnlst),//TODO: modify it
    nlp_info_(nlp_info) {
    solverInterface_ = make_shared<qpOASESInterface>(nlp_info, LP);
}

/**
 *Default destructor
 */
LPhandler::~LPhandler() {
}


/**
 * Get the optimal solution from the LPhandler_interface
 *
 *This is only an interface for user to avoid call interface directly.
 * @param p_k       the pointer to an empty array with the length equal to the size of the QP subproblem
 */
double* LPhandler::GetOptimalSolution() {
    return solverInterface_->get_optimal_solution();
}


/**
 * @brief solve the LP problem with bounds and objectives defined by the class
 * members
 */

void LPhandler::solveLP(shared_ptr<SQPhotstart::Stats> stats,
                        shared_ptr<Options> options) {
    solverInterface_->optimizeLP(stats, options);
}

void LPhandler::set_A(shared_ptr<const SpTripletMat> jacobian) {
    if (!isAinitialised) {
        I_info_A.irow1 = 1;
        I_info_A.irow2 = 1;
        I_info_A.jcol1 = nlp_info_.nVar + 1;
        I_info_A.jcol2 = nlp_info_.nVar + nlp_info_.nCon + 1;
        I_info_A.size = nlp_info_.nCon;
        solverInterface_->set_A_structure(jacobian, I_info_A);
        isAinitialised = true;
    };
    solverInterface_->set_A_values(jacobian, I_info_A);
}

void LPhandler::set_g(double rho) {
    for (int i = nlp_info_.nVar; i < nlp_info_.nVar + 2 * nlp_info_.nCon; i++)
        solverInterface_->set_g(i, rho);
}

void LPhandler::set_bounds(double delta, shared_ptr<const Vector> x_k,
                           shared_ptr<const Vector> x_l,
                           shared_ptr<const Vector> x_u,
                           shared_ptr<const Vector> c_k,
                           shared_ptr<const Vector> c_l,
                           shared_ptr<const Vector> c_u) {
    for (int i = 0; i < nlp_info_.nCon; i++) {
        solverInterface_->set_lbA(i, c_l->values()[i] - c_k->values()[i]);
        solverInterface_->set_ubA(i, c_u->values()[i] - c_k->values()[i]);
    }

    for (int i = 0; i < nlp_info_.nVar; i++) {
        solverInterface_->set_lb(i, std::max(
                                     x_l->values()[i] - x_k->values()[i], -delta));
        solverInterface_->set_ub(i, std::min(
                                     x_u->values()[i] - x_k->values()[i], delta));
    }
    /**
     * only set the upper bound for the last 2*nCon entries (those are slack variables).
     * The lower bounds are initialized as 0
     */
    for (int i = 0; i < nlp_info_.nCon * 2; i++)
        solverInterface_->set_ub(nlp_info_.nVar + i, INF);

}

const shared_ptr<QPSolverInterface> &LPhandler::getSolverInterface() const {
    return solverInterface_;
}


} // namespace SQPhotstart


