/* Copyright (C) 2019
 * All Rights Reserved.
 *
 * Authors: Xinyi Luo
 * Date:2019-07
 */

#include <sqphot/LPhandler.hpp>

namespace SQPhotstart {

LPhandler::LPhandler(NLPInfo nlp_info, shared_ptr<const Options> options,
                     Ipopt::SmartPtr<Ipopt::Journalist> jnlst) :
    QPhandler(nlp_info, options,jnlst),
    nlp_info_(nlp_info),
    jnlst_(jnlst),
    LPsolverChoice_(options->LPsolverChoice) {

#if NEW_FORMULATION
    nConstr_LP_ = nlp_info.nCon+nlp_info.nVar;
    nVar_LP_ = nlp_info.nVar*3+2*nlp_info.nCon;
    I_info_A_.length = 5;
    I_info_A_.irow = new int[5];
    I_info_A_.jcol = new int[5];
    I_info_A_.size = new int[5];
    I_info_A_.value = new double[5];
    I_info_A_.irow[0] = I_info_A_.irow[1] = 1;
    I_info_A_.irow[2] = I_info_A_.irow[3] = I_info_A_.irow[4] = nlp_info.nCon+1;
    I_info_A_.jcol[0] = nlp_info.nVar+1;
    I_info_A_.jcol[1] = nlp_info_.nVar+nlp_info.nCon+1;
    I_info_A_.jcol[2] = 1;
    I_info_A_.jcol[3] = nlp_info_.nVar+nlp_info.nCon*2+1;
    I_info_A_.jcol[4] = nlp_info_.nVar*2+nlp_info.nCon*2+1;
    I_info_A_.size[0] = I_info_A_.size[1] = nlp_info.nCon;
    I_info_A_.size[2] = I_info_A_.size[3] = I_info_A_.size[4] = nlp_info.nVar;
    I_info_A_.value[0] = I_info_A_.value[2] = I_info_A_.value[3] = 1.0;
    I_info_A_.value[1] = I_info_A_.value[4] = -1.0;
#else
    nConstr_LP_ = nlp_info.nCon;
    nVar_LP_ = nlp_info.nVar+2*nlp_info.nCon;
    I_info_A_.length = 2;
    I_info_A_.irow = new int[2];
    I_info_A_.jcol = new int[2];
    I_info_A_.size = new int[2];
    I_info_A_.value = new double[2];
    I_info_A_.irow[0] = I_info_A_.irow[1] = 1;
    I_info_A_.jcol[0] = nlp_info.nVar+1;
    I_info_A_.jcol[1] = nlp_info_.nVar+nlp_info.nCon+1;
    I_info_A_.size[0] = I_info_A_.size[1] = nlp_info.nCon;
    I_info_A_.value[0] = 1.0;
    I_info_A_.value[1] = -1.0;
#endif
    switch(LPsolverChoice_) {
    case QPOASES:
        solverInterface_ = make_shared<qpOASESInterface>(nlp_info, LP,options,jnlst);
        break;
    case QORE:
        solverInterface_ = make_shared<QOREInterface>(nlp_info,LP,options,jnlst);
        break;
    case GUROBI:
#ifdef USE_GUROBI
        solverInterface_ = make_shared<GurobiInterface>(nlp_info,LP,options,jnlst);
#endif
        break;
    case CPLEX:
#ifdef USE_CPLEX
        solverInterface_ = make_shared<CplexInterface>(nlp_info,LP,options,jnlst);
#endif
        break;
    }

}


/**
*Default destructor
*/
LPhandler::~LPhandler() {
    delete[] I_info_A_.irow;
    I_info_A_.irow = NULL;
    delete[] I_info_A_.jcol;
    I_info_A_.jcol = NULL;
    delete[] I_info_A_.size;
    I_info_A_.size = NULL;
    delete[] I_info_A_.value;
    I_info_A_.value = NULL;
}


/**
 * Get the optimal solution from the LPhandler_interface
 *
 *This is only an interface for user to avoid call interface directly.
 * @param p_k       the pointer to an empty array with the length equal to the size of the QP subproblem
 */
//double* LPhandler::get_optimal_solution() {
//    return solverInterface_->get_optimal_solution();
//
//}


/**
 * @brief solve the LP problem with bounds and objectives defined by the class
 * members
 */

void LPhandler::solveLP(shared_ptr<SQPhotstart::Stats> stats) {
    solverInterface_->optimizeLP(stats);
}

//void LPhandler::set_A(shared_ptr<const SpTripletMat> jacobian) {
//    if (!isinitialised) {
//        solverInterface_->set_A_structure(jacobian, I_info_A_);
//        isinitialised = true;
//    }
//    solverInterface_->set_A_values(jacobian, I_info_A_);
//}

void LPhandler::set_g(double rho) {
    for (int i = nlp_info_.nVar; i < nlp_info_.nVar + 2 * nlp_info_.nCon; i++)
        solverInterface_->set_g(i, rho);
}
//
//void LPhandler::set_bounds(double delta, shared_ptr<const Vector> x_l,
//                           shared_ptr<const Vector> x_u, shared_ptr<const Vector> x_k,
//                           shared_ptr<const Vector> c_l, shared_ptr<const Vector> c_u,
//                           shared_ptr<const Vector> c_k) {
//
//#if DEBUG
//#if COMPARE_QP_SOLVER
//    set_bounds_debug(delta, x_l, x_u, x_k, c_l, c_u, c_k);
//#endif
//#endif
//    /*-------------------------------------------------------------*/
//    /* Set lbA, ubA as well as lb and ub as qpOASES differentiates */
//    /*the bound constraints from the linear constraints            */
//    /*-------------------------------------------------------------*/
//    if(LPsolverChoice_ != QORE) {
//        for (int i = 0; i < nConstr_LP_; i++) {
//            solverInterface_->set_lbA(i, c_l->values()[i] - c_k->values()[i]);
//            solverInterface_->set_ubA(i, c_u->values()[i] - c_k->values()[i]);
//        }
//
//        for (int i = 0; i < nlp_info_.nVar; i++) {
//            solverInterface_->set_lb(i, std::max(
//                                         x_l->values()[i] - x_k->values()[i], -delta));
//            solverInterface_->set_ub(i, std::min(
//                                         x_u->values()[i] - x_k->values()[i], delta));
//        }
//        /**
//         * only set the upper bound for the last 2*nCon entries (those are slack variables).
//         * The lower bounds are initialized as 0
//         */
//        for (int i = 0; i < nlp_info_.nCon * 2; i++)
//            solverInterface_->set_ub(nlp_info_.nVar + i, INF);
//
//    }
//
//    /*-------------------------------------------------------------*/
//    /* Only set lb and ub, where lb = [lbx;lbA]; and ub=[ubx; ubA] */
//    /*-------------------------------------------------------------*/
//    else {
//        for (int i = 0; i < nlp_info_.nVar; i++) {
//            solverInterface_->set_lb(i, std::max(
//                                         x_l->values()[i] - x_k->values()[i], -delta));
//            solverInterface_->set_ub(i, std::min(
//                                         x_u->values()[i] - x_k->values()[i], delta));
//
//        }
//
//        for (int i = 0; i < nlp_info_.nCon * 2; i++)
//            solverInterface_->set_ub(nlp_info_.nVar + i, INF);
//
//        for (int i = 0; i < nlp_info_.nCon; i++) {
//            solverInterface_->set_lb(nlp_info_.nVar +2*nlp_info_.nCon+i, c_l->values()[i]
//                                     - c_k->values()[i]);
//            solverInterface_->set_ub(nlp_info_.nVar +2*nlp_info_.nCon+i, c_u->values()[i]
//                                     - c_k->values()[i]);
//        }
//    }
//
//}




} // namespace SQPhotstart


