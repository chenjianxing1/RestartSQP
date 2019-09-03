/** Copyright (C) 2019
 * All Rights Reserved.
 *
 * Authors: Xinyi Luo
 * Date:2019-07
 */
#include <sqphot/QPhandler.hpp>


namespace SQPhotstart {
using namespace std;

QPhandler::QPhandler(Index_info nlp_info, shared_ptr<const Options> options,
                     Ipopt::SmartPtr<Ipopt::Journalist> jnlst) :
    nlp_info_(nlp_info),
    jnlst_(jnlst),
    nConstr_QP_(nlp_info.nCon),
    nVar_QP_(nlp_info.nVar +2*nlp_info.nCon),
    QPsolverChoice_(options->QPsolverChoice) {
#if DEBUG
#if COMPARE_QP_SOLVER
    qpOASESInterface_ = make_shared<qpOASESInterface>(nlp_info, QP,options);
    QOREInterface_= make_shared<QOREInterface>(nlp_info,QP,options,jnlst);
    W_b_qpOASES_ = new ActiveType[nVar_QP_];
    W_c_qpOASES_ = new ActiveType[nConstr_QP_];
    W_b_qore_ = new ActiveType[nVar_QP_];
    W_c_qore_ = new ActiveType[nConstr_QP_];
#endif
#else
    switch(QPsolverChoice_) {
    case QPOASES_QP:
        solverInterface_ = make_shared<qpOASESInterface>(nlp_info, QP,options);
        break;
    case QORE_QP:
        solverInterface_ = make_shared<QOREInterface>(nlp_info,QP,options,jnlst);
        break;
    default:
        THROW_EXCEPTION(INVALID_QP_SOLVER_CHOICE,"The QP solver choice is invalid!")
    }
#endif
}




/**
 *Default destructor
 */
QPhandler::~QPhandler() {
#if DEBUG
#if COMPARE_QP_SOLVER
    delete[] W_b_qpOASES_;
    delete[] W_c_qpOASES_;
    delete[] W_b_qore_;
    delete[] W_c_qore_;
#endif
#endif

}


/**
 * Get the optimal solution from the QPhandler_interface
 *
 *This is only an interface for user to avoid call interface directly.
 * @param p_k       the pointer to an empty array with the length equal to the size
 *                  of the QP subproblem
 */
double* QPhandler::GetOptimalSolution() {
#if DEBUG
    return qpOASESInterface_->get_optimal_solution();
#else
    return solverInterface_->get_optimal_solution();
#endif
}


/**
 *Get the multipliers from the QPhandler_interface
 *
 *This is only an interface for user to avoid call interface directly.
 *
 * @param y_k       the pointer to an empty array with the length equal to the size of
 * multipliers of the QP subproblem
 */
double*  QPhandler::GetMultipliers() {
#if DEBUG

    return qpOASESInterface_->get_multipliers();
#else
    return solverInterface_->get_multipliers();
#endif

}


/**
 * Setup the bounds for the QP subproblems according to the information from current
 * iterate. We have
 * 	c_l -c_k <=J_p+ u-v<=c_u-c_k
 * The bound is formulated as
 *   max(-delta, x_l-x_k)<=p<=min(delta,x_u-x_k)
 * and  u,v>=0
 * @param delta      trust region radius
 * @param x_k        current iterate point
 * @param c_k        current constraint value evaluated at x_k
 * @param x_l        the lower bounds for variables
 * @param x_u        the upper bounds for variables
 * @param c_l        the lower bounds for constraints
 * @param c_u        the upper bounds for constraints
 */

void QPhandler::set_bounds(double delta, shared_ptr<const Vector> x_l,
                           shared_ptr<const Vector> x_u, shared_ptr<const Vector> x_k,
                           shared_ptr<const Vector> c_l, shared_ptr<const Vector> c_u,
                           shared_ptr<const Vector> c_k) {


#if DEBUG
#if COMPARE_QP_SOLVER
    set_bounds_debug(delta, x_l, x_u, x_k, c_l, c_u, c_k);
#endif
#else
    /*-------------------------------------------------------------*/
    /* Set lbA, ubA as well as lb and ub as qpOASES differentiates */
    /*the bound constraints from the linear constraints            */
    /*-------------------------------------------------------------*/
    if(QPsolverChoice_==QPOASES_QP) {
        for (int i = 0; i < nConstr_QP_; i++) {
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

    /*-------------------------------------------------------------*/
    /* Only set lb and ub, where lb = [lbx;lbA]; and ub=[ubx; ubA] */
    /*-------------------------------------------------------------*/
    else if(QPsolverChoice_==QORE_QP) {
        for (int i = 0; i < nlp_info_.nVar; i++) {
            solverInterface_->set_lb(i, std::max(
                                         x_l->values()[i] - x_k->values()[i], -delta));
            solverInterface_->set_ub(i, std::min(
                                         x_u->values()[i] - x_k->values()[i], delta));

        }

        for (int i = 0; i < nlp_info_.nCon * 2; i++)
            solverInterface_->set_ub(nlp_info_.nVar + i, INF);

        for (int i = 0; i < nlp_info_.nCon; i++) {
            solverInterface_->set_lb(nlp_info_.nVar +2*nlp_info_.nCon+i, c_l->values()[i]
                                     - c_k->values()[i]);
            solverInterface_->set_ub(nlp_info_.nVar +2*nlp_info_.nCon+i, c_u->values()[i]
                                     - c_k->values()[i]);
        }
    }
#endif
}


/**
 * This function sets up the object vector g of the QP problem
 * The (2*nCon+nVar) vector g_^T in QP problem will be the same as
 * [grad_f^T, rho* e^T], where the unit vector is of length (2*nCon).
 * @param grad      Gradient vector from nlp class
 * @param rho       Penalty Parameter
 */

void QPhandler::set_g(shared_ptr<const Vector> grad, double rho) {

#if DEBUG
#if COMPARE_QP_SOLVER
    for (int i = 0; i < nlp_info_.nVar + 2 * nlp_info_.nCon; i++)
        if (i < nlp_info_.nVar) {
            qpOASESInterface_->set_g(i, grad->values()[i]);
            QOREInterface_->set_g(i, grad->values()[i]);
        }
        else {
            qpOASESInterface_->set_g(i, rho);
            QOREInterface_->set_g(i,rho);
        }

#endif
#else
    if(QPsolverChoice_==QPOASES_QP||QPsolverChoice_==QORE_QP) {
        for (int i = 0; i < nlp_info_.nVar + 2 * nlp_info_.nCon; i++)
            if (i < nlp_info_.nVar)
                solverInterface_->set_g(i, grad->values()[i]);

            else
                solverInterface_->set_g(i, rho);

    }
#endif

}

/**
 * @brief Set up the H for the first time in the QP problem.
 * It will be concatenated as [H_k 0]
 *          		         [0   0]
 * where H_k is the Lagragian hessian evaluated at x_k and lambda_k.
 *
 * This method should only be called for once.
 *
 * @param hessian the Lagragian hessian evaluated at x_k and lambda_k from nlp
 * readers.
 */
void QPhandler::set_H(shared_ptr<const SpTripletMat> hessian) {
#if DEBUG
#if COMPARE_QP_SOLVER
    qpOASESInterface_->set_H_structure(hessian);
    qpOASESInterface_->set_H_values(hessian);
    QOREInterface_->set_H_structure(hessian);
    QOREInterface_->set_H_values(hessian);
#endif
#else
    solverInterface_->set_H_structure(hessian);
    solverInterface_->set_H_values(hessian);
#endif
}


/**
 * @brief This function sets up the matrix A in the QP subproblem
 * The matrix A in QP problem will be concatenate as [J I -I]
 * @param jacobian  the Matrix object for Jacobian from c(x)
 */
void QPhandler::set_A(shared_ptr<const SpTripletMat> jacobian) {
    I_info_A.irow1 = 1;
    I_info_A.irow2 = 1;
    I_info_A.jcol1 = nlp_info_.nVar + 1;
    I_info_A.jcol2 = nlp_info_.nVar + nlp_info_.nCon + 1;
    I_info_A.size = nlp_info_.nCon;
#if DEBUG
#if COMPARE_QP_SOLVER
    qpOASESInterface_->set_A_structure(jacobian, I_info_A);
    qpOASESInterface_->set_A_values(jacobian, I_info_A);
    QOREInterface_->set_A_structure(jacobian, I_info_A);
    QOREInterface_->set_A_values(jacobian, I_info_A);
#endif
#else
    solverInterface_->set_A_structure(jacobian, I_info_A);
    solverInterface_->set_A_values(jacobian, I_info_A);
#endif
}



/**
 * @brief This function updates the constraint if there is any changes to
 * the iterates
 */
void QPhandler::update_bounds(double delta, shared_ptr<const Vector> x_l,
                              shared_ptr<const Vector> x_u,
                              shared_ptr<const Vector> x_k,
                              shared_ptr<const Vector> c_l,
                              shared_ptr<const Vector> c_u,
                              shared_ptr<const Vector> c_k) {
#if DEBUG
#if COMPARE_QP_SOLVER
    set_bounds_debug(delta, x_l, x_u, x_k, c_l, c_u, c_k);
#endif
#else
    if(QPsolverChoice_==QPOASES_QP) {
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
    }
    else if(QPsolverChoice_== QORE_QP) {
        for (int i = 0; i < nlp_info_.nVar; i++) {
            solverInterface_->set_lb(i, std::max(
                                         x_l->values()[i] - x_k->values()[i], -delta));
            solverInterface_->set_ub(i, std::min(
                                         x_u->values()[i] - x_k->values()[i], delta));

        }
        for (int i = 0; i < nlp_info_.nCon; i++) {
            solverInterface_->set_lb(nlp_info_.nVar +2*nlp_info_.nCon+i, c_l->values()[i]
                                     - c_k->values()[i]);
            solverInterface_->set_ub(nlp_info_.nVar +2*nlp_info_.nCon+i, c_u->values()[i]
                                     - c_k->values()[i]);
        }
    }
#endif
}


/**
 * @brief This function updates the vector g in the QP subproblem when there are any
 * change to the values of penalty parameter
 *
 * @param rho               penalty parameter
 * @param nVar              number of variables in NLP
 */

void QPhandler::update_penalty(double rho) {
#if DEBUG
#if COMPARE_QP_SOLVER
    for (int i = nlp_info_.nVar; i < nlp_info_.nVar + nlp_info_.nCon * 2; i++) {
        qpOASESInterface_->set_g(i, rho);
        QOREInterface_->set_g(i, rho);
    }
#endif
#else
    if(QPsolverChoice_==QPOASES_QP||QPsolverChoice_== QORE_QP) {
        for (int i = nlp_info_.nVar; i < nVar_QP_; i++)
            solverInterface_->set_g(i, rho);
    }
#endif
}


/**
 * @brief This function updates the vector g in the QP subproblem when there are any
 * change to the values of gradient in NLP
 *
 * @param grad              the gradient vector from NLP
 */
void QPhandler::update_grad(shared_ptr<const Vector> grad) {

#if DEBUG
#if COMPARE_QP_SOLVER

    for (int i = 0; i < nlp_info_.nVar; i++) {
        qpOASESInterface_->set_g(i, grad->values()[i]);
        QOREInterface_->set_g(i, grad->values()[i]);
    }
#endif
#else
    if(QPsolverChoice_==QPOASES_QP||QPsolverChoice_== QORE_QP) {
        for (int i = 0; i < nlp_info_.nVar; i++)
            solverInterface_->set_g(i, grad->values()[i]);
    }
#endif
}



/**
 *@brief Solve the QP with objective and constraints defined by its class members
 */
void QPhandler::solveQP(shared_ptr<SQPhotstart::Stats> stats,
                        shared_ptr<Options> options) {

#if DEBUG
#if COMPARE_QP_SOLVER
    QOREInterface_->optimizeQP(stats);
    qpOASESInterface_->optimizeQP(stats);
    bool qpOASES_optimal = OptimalityTest(QPOASES_QP);
    bool qore_optimal = OptimalityTest(QORE_QP);
//        OptimalityTest(QPOASES_QP, , nullptr, shared_ptr<const
//                Vector>(),
//                       shared_ptr<const Vector>(), shared_ptr<const Vector>(),
//                       shared_ptr<const Vector>(), 0, nullptr, shared_ptr<const Vector>(),
//                       shared_ptr<const Vector>());

    if(!qpOASES_optimal||!qore_optimal)    
        testQPsolverDifference();


#endif
#else
    solverInterface_->optimizeQP(stats);
#endif
}


double QPhandler::GetObjective() {

#if DEBUG
#if COMPARE_QP_SOLVER

    return qpOASESInterface_->get_obj_value();
#endif
#else
    return solverInterface_->get_obj_value();
#endif
}


void QPhandler::update_H(shared_ptr<const SpTripletMat> Hessian) {

#if DEBUG
#if COMPARE_QP_SOLVER
    QOREInterface_->set_H_values(Hessian);
    qpOASESInterface_->set_H_values(Hessian);
#endif
#else
    solverInterface_->set_H_values(Hessian);
#endif
}


void QPhandler::update_A(shared_ptr<const SpTripletMat> Jacobian) {

#if DEBUG
#if COMPARE_QP_SOLVER
    QOREInterface_->set_A_values(Jacobian, I_info_A);
    qpOASESInterface_->set_A_values(Jacobian, I_info_A);
#endif
#else
    solverInterface_->set_A_values(Jacobian, I_info_A);
#endif
}


    void QPhandler::update_delta(double delta, shared_ptr<const Vector> x_l,
                             shared_ptr<const Vector> x_u,
                             shared_ptr<const Vector> x_k) {

#if DEBUG
#if COMPARE_QP_SOLVER

    for (int i = 0; i < nlp_info_.nVar; i++) {
        qpOASESInterface_->set_lb(i, std::max(
                                      x_l->values()[i] - x_k->values()[i], -delta));
        qpOASESInterface_->set_ub(i, std::min(
                                      x_u->values()[i] - x_k->values()[i], delta));
        QOREInterface_->set_lb(i, std::max(
                                   x_l->values()[i] - x_k->values()[i], -delta));
        QOREInterface_->set_ub(i, std::min(
                                   x_u->values()[i] - x_k->values()[i], delta));
    }
#endif
#else
    if(QPsolverChoice_==QPOASES_QP||QPsolverChoice_== QORE_QP) {
        for (int i = 0; i < nlp_info_.nVar; i++) {
            solverInterface_->set_lb(i, std::max(
                                         x_l->values()[i] - x_k->values()[i], -delta));
            solverInterface_->set_ub(i, std::min(
                                         x_u->values()[i] - x_k->values()[i], delta));
        }
    }
#endif

}


#if DEBUG
#if COMPARE_QP_SOLVER
void QPhandler::set_bounds_debug(double delta, shared_ptr<const Vector> x_l,
                                 shared_ptr<const Vector> x_u,
                                 shared_ptr<const Vector> x_k,
                                 shared_ptr<const Vector> c_l,
                                 shared_ptr<const Vector> c_u,
                                 shared_ptr<const Vector> c_k) {

    for (int i = 0; i < nlp_info_.nCon; i++) {
        qpOASESInterface_->set_lbA(i, c_l->values()[i] - c_k->values()[i]);
        qpOASESInterface_->set_ubA(i, c_u->values()[i] - c_k->values()[i]);
    }

    for (int i = 0; i < nlp_info_.nVar; i++) {
        qpOASESInterface_->set_lb(i, std::max(
                                      x_l->values()[i] - x_k->values()[i], -delta));
        qpOASESInterface_->set_ub(i, std::min(
                                      x_u->values()[i] - x_k->values()[i], delta));
    }
    /**
     * only set the upper bound for the last 2*nCon entries (those are slack variables).
     * The lower bounds are initialized as 0
     */
    for (int i = 0; i < nlp_info_.nCon * 2; i++)
        qpOASESInterface_->set_ub(nlp_info_.nVar + i, INF);


    /*-------------------------------------------------------------*/
    /* Only set lb and ub, where lb = [lbx;lbA]; and ub=[ubx; ubA] */
    /*-------------------------------------------------------------*/
    for (int i = 0; i < nlp_info_.nVar; i++) {
        QOREInterface_->set_lb(i, std::max(
                                   x_l->values()[i] - x_k->values()[i], -delta));
        QOREInterface_->set_ub(i, std::min(
                                   x_u->values()[i] - x_k->values()[i], delta));

    }

    for (int i = 0; i < nConstr_QP_ * 2; i++)
        QOREInterface_->set_ub(nlp_info_.nVar + i, INF);

    for (int i = 0; i < nConstr_QP_; i++) {
        QOREInterface_->set_lb(nVar_QP_+i, c_l->values()[i]
                               - c_k->values()[i]);
        QOREInterface_->set_ub(nVar_QP_+i, c_u->values()[i]
                               - c_k->values()[i]);
    }

}

bool QPhandler::testQPsolverDifference() {
    shared_ptr<Vector> qpOASESsol = make_shared<Vector>(nlp_info_.nVar);
    shared_ptr<Vector> QOREsol = make_shared<Vector>(nlp_info_.nVar);
    shared_ptr<Vector> difference = make_shared<Vector>(nlp_info_.nVar);
    qpOASESsol->print("qpOASESsol",jnlst_);
    QOREsol->print("QOREsol",jnlst_);
    qpOASESsol->copy_vector(qpOASESInterface_->get_optimal_solution());
    QOREsol->copy_vector(QOREInterface_->get_optimal_solution());
    difference->copy_vector(qpOASESsol);
    difference->subtract_vector(QOREsol->values());
    double diff_norm = difference->getOneNorm();
    if(diff_norm>1.0e-8) {
        printf("difference is %10e\n",diff_norm);
        qpOASESsol->print("qpOASESsol");
        QOREsol->print("QOREsol");
        QOREInterface_->WriteQPDataToFile(jnlst_,J_ALL,J_DBG);
    }
    assert(diff_norm<1.0e-8);
    return true;

}

bool QPhandler::OptimalityTest(QPSolver qpSolver) {

    int i;
//create local variables and set all violation values to be 0
    double primal_violation = 0;
    double dual_violation = 0;
    double compl_violation = 0;
    double statioanrity_violation = 0;

    //create two temporary vector for storage some data if they are needed
    shared_ptr<Vector> Ax = make_shared<Vector>(nConstr_QP_);
    shared_ptr<Vector> stationary_gap = make_shared<Vector>(nVar_QP_);
    shared_ptr<Vector> multiplier_constr=make_shared<Vector> (nConstr_QP_);
    shared_ptr<Vector> multiplier_bounds=make_shared<Vector> (nVar_QP_);
    shared_ptr<Vector> x = make_shared<Vector>(nVar_QP_);
    if (qpSolver == QPOASES_QP) {
        auto lb = qpOASESInterface_->getLb();
        auto ub = qpOASESInterface_->getUb();
        auto lbA =qpOASESInterface_->getLbA();
        auto ubA =qpOASESInterface_->getUbA();
        auto g=qpOASESInterface_->getG();

        shared_ptr<const SpTripletMat> A = qpOASESInterface_->getA();
        shared_ptr<const SpTripletMat> H = qpOASESInterface_->getH();
        x->copy_vector(qpOASESInterface_->get_optimal_solution());
        multiplier_bounds->copy_vector(qpOASESInterface_->get_multipliers());
        multiplier_constr->copy_vector(qpOASESInterface_->get_multipliers()+nVar_QP_);
        qpOASESInterface_->GetWorkingSet(W_c_qpOASES_, W_b_qpOASES_);
        /**-------------------------------------------------------**/
        /**                    primal feasibility                 **/
        /**-------------------------------------------------------**/
        assert(lbA != nullptr && ubA != nullptr);
        for (i = 0; i < nVar_QP_; i++) {
            primal_violation += max(0.0, (lb->values()[i] - x->values()[i]));
            primal_violation += -min(0.0, (ub->values()[i] - x->values()[i]));
        }
        if (A != nullptr) {
            A->times(x, Ax); //tmp_vec_nCon=A*x
            for (i = 0; i < nConstr_QP_; i++) {
                primal_violation += max(0.0, (lbA->values()[i] -
                                              Ax->values()[i]));
                primal_violation += -min(0.0, (ubA->values()[i] -
                                               Ax->values()[i]));
            }
        }
        /**-------------------------------------------------------**/
        /**                    dual feasibility                   **/
        /**-------------------------------------------------------**/
        for (i = 0; i < nVar_QP_; i++) {
            switch (W_b_qpOASES_[i]) {
            case INACTIVE://the constraint is inactive, then the dual multiplier
                // should be 0
                dual_violation += fabs(multiplier_bounds->values()[i]);
                break;
            case ACTIVE_BELOW://the constraint is active at the lower bound, so the
                // multiplier should be positive
                dual_violation += -min(0.0, multiplier_bounds->values()[i]);
                break;
            case ACTIVE_ABOVE: //the contraint is active at the upper bounds, so the
                // multiplier should be negavie
                dual_violation += max(0.0, multiplier_bounds->values()[i]);
                break;
            default:
                printf("failed in dual fea test, the working set  for qpOASES at "
                       "the "
                       "bounds is %i", W_b_qore_[i]);
                    THROW_EXCEPTION(INVALID_WORKING_SET, INVALID_WORKING_SET_MSG);
            }
        }
        if(A!= nullptr) {
            for (i = 0; i < nConstr_QP_; i++) {
                switch (W_c_qpOASES_[i]) {
                case INACTIVE://the constraint is inactive, then the dual multiplier
                    // should be 0
                    dual_violation += fabs(multiplier_constr->values()[i]);
                    break;
                case ACTIVE_BELOW://the constraint is active at the lower bound, so the
                    // multiplier should be positive
                    dual_violation += -min(0.0, multiplier_constr->values()[i]);
                    break;
                case ACTIVE_ABOVE: //the contraint is active at the upper bounds, so the
                    // multiplier should be negavie
                    dual_violation += max(0.0, multiplier_constr->values()[i]);
                    break;
                default:
                    printf("failed in dual fea test, the working set  for qpOASES at "
                           "the "
                           "constr is %i", W_c_qore_[i]);
                    THROW_EXCEPTION(INVALID_WORKING_SET, INVALID_WORKING_SET_MSG);
                }
            }
        }
        /**-------------------------------------------------------**/
        /**                   stationarity                        **/
        /**-------------------------------------------------------**/
        //calculate A'*y+lambda-(g+Hx)
        if(A!=nullptr) {
            A->transposed_times(multiplier_constr, stationary_gap);
        }
        shared_ptr<Vector> Hx = make_shared<Vector>(nVar_QP_);
        H->times(x,Hx);

        stationary_gap->add_vector(multiplier_bounds->values());
        stationary_gap->subtract_vector(g->values());
        stationary_gap->subtract_vector(Hx->values());
        statioanrity_violation += stationary_gap->getOneNorm();


        /**-------------------------------------------------------**/
        /**                    Complemtarity                      **/
        /**-------------------------------------------------------**/

        for (i = 0; i < nVar_QP_; i++) {
            switch (W_b_qpOASES_[i]) {
            case INACTIVE: //constraint is inactive, multiplier should be 0
                compl_violation += abs(multiplier_bounds->values()[i]);
                break;
            case ACTIVE_BELOW://the constraint is active at the lower bound
                compl_violation += abs(multiplier_bounds->values()[i] *
                                       (x->values()[i] - lb->values()[i]));
                break;
            case ACTIVE_ABOVE: //the contraint is active at the upper bounds, so the
                // multiplier should be negavie
                compl_violation += abs(multiplier_bounds->values()[i] *
                                       (ub->values()[i] - x->values()[i]));
                break;
            default:
                printf("failed in compl test, the working set  for qpOASES at "
                       "the "
                       "bounds is %i", W_b_qore_[i]);
                    THROW_EXCEPTION(INVALID_WORKING_SET, INVALID_WORKING_SET_MSG);
            }
        }
        if (A != nullptr) {
            for (i = 0; i < nConstr_QP_; i++) {
                switch (W_c_qpOASES_[i]) {
                case INACTIVE: //constraint is inactive, multiplier should be 0
                    compl_violation += abs(multiplier_constr->values()[i]);
                    break;
                case ACTIVE_BELOW://the constraint is active at the lower bound
                    compl_violation += abs(multiplier_constr->values()[i] *
                                           (Ax->values()[i] - lbA->values()[i]));
                    break;
                case ACTIVE_ABOVE: //the contraint is active at the upper bounds, so the
                    // multiplier should be negavie
                    compl_violation += abs(multiplier_constr->values()[i] *
                                           (ubA->values()[i] - Ax->values()[i]));
                    break;
                default:
                    printf("failed in compl test, the working set  for qpOASES at "
                           "the "
                           "constr is %i", W_c_qore_[i]);
                        THROW_EXCEPTION(INVALID_WORKING_SET, INVALID_WORKING_SET_MSG);
                }
            }
        }
    }

    else if(qpSolver == QORE_QP) {
        auto lb = QOREInterface_->getLb();
        auto ub = QOREInterface_->getUb();
        auto g=QOREInterface_->getG();

        x->copy_vector(QOREInterface_->get_optimal_solution());
        shared_ptr<const SpTripletMat> A = QOREInterface_->getA();
        shared_ptr<const SpTripletMat> H = QOREInterface_->getH();
        multiplier_bounds->copy_vector(QOREInterface_->get_multipliers());
        multiplier_constr->copy_vector(QOREInterface_->get_multipliers()+nVar_QP_);
        QOREInterface_->GetWorkingSet(W_c_qore_, W_b_qore_);

        /**-------------------------------------------------------**/
        /**                    primal feasibility                 **/
        /**-------------------------------------------------------**/
        for (i = 0; i < nVar_QP_; i++) {
            primal_violation += max(0.0, (lb->values()[i] - x->values()[i]));
            primal_violation += -min(0.0, (ub->values()[i] - x->values()[i]));
        }
        if (A != nullptr) {
            A->times(x, Ax); //tmp_vec_nCon=A*x
            for (i = 0; i < nConstr_QP_; i++) {
                primal_violation += max(0.0, (lb->values()[i + nVar_QP_] -
                                              Ax->values()[i]));
                primal_violation += -min(0.0, (ub->values()[i + nVar_QP_] -
                                               Ax->values()[i]));
            }
        }
        /**-------------------------------------------------------**/
        /**                    dual feasibility                   **/
        /**-------------------------------------------------------**/
        for (i = 0; i < nVar_QP_; i++) {
            switch (W_b_qpOASES_[i]) {
            case INACTIVE://the constraint is inactive, then the dual multiplier
                // should be 0
                dual_violation += fabs(multiplier_bounds->values()[i]);
                break;
            case ACTIVE_BELOW://the constraint is active at the lower bound, so the
                // multiplier should be positive
                dual_violation += -min(0.0, multiplier_bounds->values()[i]);
                break;
            case ACTIVE_ABOVE: //the contraint is active at the upper bounds, so the
                // multiplier should be negavie
                dual_violation += max(0.0, multiplier_bounds->values()[i]);
                break;
            default:
                printf("failed in dual fea test, the working set  for qore at the "
                       "bounds is %i", W_b_qpOASES_[i]);
                THROW_EXCEPTION(INVALID_WORKING_SET, INVALID_WORKING_SET_MSG);
            }
        }
        if (A != nullptr) {
            for (i = 0; i < nConstr_QP_; i++) {
                switch (W_c_qpOASES_[i]) {
                case INACTIVE://the constraint is inactive, then the dual multiplier
                    // should be 0
                    dual_violation += fabs(multiplier_constr->values()[i]);
                    break;
                case ACTIVE_BELOW://the constraint is active at the lower bound, so the
                    // multiplier should be positive
                    dual_violation += -min(0.0, multiplier_constr->values()[i]);
                    break;
                case ACTIVE_ABOVE: //the contraint is active at the upper bounds, so the
                    // multiplier should be negavie
                    dual_violation += max(0.0, multiplier_constr->values()[i]);
                    break;
                default:
                    printf("failed in dual fea test, the working set  for qore at the "
                           "constr is %i", W_c_qpOASES_[i]);
                    THROW_EXCEPTION(INVALID_WORKING_SET, INVALID_WORKING_SET_MSG);
                }
            }
        }
        /**-------------------------------------------------------**/
        /**                   stationarity                        **/
        /**-------------------------------------------------------**/
        //calculate A'*y+lambda-g-Hx
        if (A != nullptr) {
            A->transposed_times(multiplier_constr, stationary_gap);
        }
        shared_ptr<Vector> Hx = make_shared<Vector>(nVar_QP_);
        H->times(x, Hx);

        stationary_gap->add_vector(multiplier_bounds->values());
        stationary_gap->subtract_vector(g->values());
        stationary_gap->subtract_vector(Hx->values());
        statioanrity_violation += stationary_gap->getOneNorm();


        /**-------------------------------------------------------**/
        /**                    Complemtarity                      **/
        /**-------------------------------------------------------**/
        for (i = 0; i < nVar_QP_; i++) {
            switch (W_b_qpOASES_[i]) {
            case INACTIVE: //constraint is inactive, multiplier should be 0
                compl_violation += abs(multiplier_bounds->values()[i]);
                break;
            case ACTIVE_BELOW://the constraint is active at the lower bound
                compl_violation += abs(multiplier_bounds->values()[i] *
                                       (x->values()[i] - lb->values()[i]));
                break;
            case ACTIVE_ABOVE: //the contraint is active at the upper bounds, so the
                // multiplier should be negavie
                compl_violation += abs(multiplier_bounds->values()[i] *
                                       (ub->values()[i] - x->values()[i]));
                break;
            default:
                printf("failed in compl test, the working set  for qore at "
                       "the "
                       "bounds is %i", W_b_qpOASES_[i]);
                    THROW_EXCEPTION(INVALID_WORKING_SET, INVALID_WORKING_SET_MSG);
            }
        }
        if (A != nullptr) {
            for (i = 0; i < nConstr_QP_; i++) {
                switch (W_c_qpOASES_[i]) {
                case INACTIVE: //constraint is inactive, multiplier should be 0
                    compl_violation += abs(multiplier_constr->values()[i]);
                    break;
                case ACTIVE_BELOW://the constraint is active at the lower bound
                    compl_violation += abs(multiplier_constr->values()[i] *
                                           (Ax->values()[i] - lb->values()
                                            [i + nVar_QP_]));
                    break;
                case ACTIVE_ABOVE: //the contraint is active at the upper bounds, so the
                    // multiplier should be negavie
                    compl_violation += abs(multiplier_constr->values()[i] *
                                           (ub->values()[i + nVar_QP_] -
                                            Ax->values()
                                            [i]));
                    break;
                default:
                    printf("failed in compl test, the working set  for qpOASES at "
                           "the "
                           "constr is %i", W_c_qpOASES_[i]);
                        THROW_EXCEPTION(INVALID_WORKING_SET, INVALID_WORKING_SET_MSG);
                }
            }
        }
    }
    /**-------------------------------------------------------**/
    /**             Decide if x_k is optimal                  **/
    /**-------------------------------------------------------**/

    qpOptimalStatus_.compl_violation=compl_violation;
    qpOptimalStatus_.statioanrity_violation=statioanrity_violation;
    qpOptimalStatus_.dual_violation=dual_violation;
    qpOptimalStatus_.primal_feasibility=primal_violation;

if(qpOptimalStatus_.KKT_error>1.0e-6) {
    printf("comp_violation %10e\n", compl_violation);
    printf("stat_violation %10e\n", statioanrity_violation);
    printf("prim_violation %10e\n", primal_violation);
    printf("dual_violation %10e\n", dual_violation);
    qpOptimalStatus_.KKT_error =
            compl_violation + statioanrity_violation + dual_violation + primal_violation;
    assert(qpOptimalStatus_.KKT_error < 1.0e-4);
    return false;
}

return true;
}


#endif
#endif
} // namespace SQPhotstart







