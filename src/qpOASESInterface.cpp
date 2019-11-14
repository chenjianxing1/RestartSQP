/* Copyright (C) 2019
 * All Rights Reserved.
 *
 * Authors: Xinyi Luo
 * Date:    2019-07
 */
#include "sqphot/qpOASESInterface.hpp"

using namespace std;

namespace SQPhotstart {


/**
 * @brief Constructor which also initializes the qpOASES SQProblem objects
 *
 * It will either formulate the SQP subproblem as
 * 		\minimize 1/2*p^T*H_k*p + g_k^T*p + rho*(e^T u_1 + e^Tu_2 +e^T v_1 e^T v_2)
 * 		\subject  c_l-c_k <= J_k*p + u_1 - u_2 <= c_u-c_k,
 * 			  x_l-x_k <= p + v_1 -v_2 <= x_u -x_k, *
 * 			  ||p||<=\delta,
 * 			  u_1,u_2,v_1,v_2>=0.
 * Or
 * 		\minimize 1/2*p^T*H_k*p + g_k^T*p + rho*(e^T u_1 + e^Tu_2)
 * 		\subject  c_l-c_k <= J_k*p + u_1 - u_2 <= c_u-c_k,
 * 			  x_l-x_k <= p <= x_u -x_k, *
 * 			  ||p||<=\delta,
 * 			  u_1,u_2.
 *
 * depending on user's choice.
 *
 * @param nlp_info the struct that stores simple nlp dimension info
 * @param qptype  is the problem to be solved QP or LP or SOC?
 */
qpOASESInterface::qpOASESInterface(NLPInfo nlp_info, QPType qptype,
                                   shared_ptr<const Options> options,
                                   Ipopt::SmartPtr<Ipopt::Journalist> jnlst):
    jnlst_(jnlst),
    options_(options)
{

#if NEW_FORMULATION
    nConstr_QP_ = nlp_info.nCon+nlp_info.nVar;
    nVar_QP_ = nlp_info.nVar*3+2*nlp_info.nCon;
#else
    nConstr_QP_ = nlp_info.nCon;
    nVar_QP_ = nlp_info.nVar+2*nlp_info.nCon;
#endif
    allocate_memory(nlp_info, qptype);
}



qpOASESInterface::qpOASESInterface(shared_ptr<SpHbMat> H,
                                   shared_ptr<SpHbMat> A,
                                   shared_ptr<Vector> g,
                                   shared_ptr<Vector> lb,
                                   shared_ptr<Vector> ub,
                                   shared_ptr<Vector> lbA,
                                   shared_ptr<Vector> ubA,
                                   shared_ptr<Options> options):
    nVar_QP_(A->ColNum()),
    nConstr_QP_(A->RowNum()),
    H_(H),
    A_(A),
    g_(g),
    lb_(lb),
    ub_(ub),
    lbA_(lbA),
    ubA_(ubA),
    options_(options)
{

    x_qp_ = make_shared<Vector>(nVar_QP_);
    y_qp_ = make_shared<Vector>(nConstr_QP_+nVar_QP_);
    solver_ = std::make_shared<qpOASES::SQProblem>((qpOASES::int_t) nVar_QP_,
              (qpOASES::int_t) nConstr_QP_);
    H_qpOASES_ = std::make_shared<qpOASES::SymSparseMat>(H_->RowNum(),
                 H_->ColNum(),
                 H_->RowIndex(),
                 H_->ColIndex(),
                 H_->MatVal());
    H_qpOASES_->createDiagInfo();

    A_qpOASES_ = std::make_shared<qpOASES::SymSparseMat>(A_->RowNum(),
                 A_->ColNum(),
                 A_->RowIndex(),
                 A_->ColIndex(),
                 A_->MatVal());



}

/**Default destructor*/
qpOASESInterface::~qpOASESInterface() = default;


/**
 * @brief Allocate memory for the class members
 * @param nlp_info  the struct that stores simple nlp dimension info
 * @param qptype is the problem to be solved QP or LP or SOC?
 * @return
 */
void qpOASESInterface::allocate_memory(NLPInfo nlp_info, QPType qptype) {

#if NEW_FORMULATION
    int nnz_g_QP = nlp_info.nnz_jac_g+2*nlp_info.nCon+3*nlp_info.nVar;
#else
    int nnz_g_QP = nlp_info.nnz_jac_g + 2 * nlp_info.nCon;
#endif
    lbA_ = make_shared<Vector>(nConstr_QP_);
    ubA_ = make_shared<Vector>(nConstr_QP_);
    lb_ = make_shared<Vector>(nVar_QP_);
    ub_ = make_shared<Vector>(nVar_QP_);
    g_ = make_shared<Vector>(nVar_QP_);
    A_ = make_shared<SpHbMat>(nnz_g_QP, nConstr_QP_, nVar_QP_,false);
    x_qp_ = make_shared<Vector>(nVar_QP_);
    y_qp_ = make_shared<Vector>(nConstr_QP_+nVar_QP_);

    if (qptype != LP) {
        H_ = make_shared<SpHbMat>(nVar_QP_, nVar_QP_, false);
    }

    solver_ = std::make_shared<qpOASES::SQProblem>((qpOASES::int_t) nVar_QP_,
              (qpOASES::int_t) nConstr_QP_);
}

/**
 * @brief This method solves the QP problem specified in the data, with given
 * options.
 * After the QP being solved, it updates the stats, adding the iteration
 * number used to solve the QP to the qp_iter in object stats
 */
void
qpOASESInterface::optimizeQP(shared_ptr<Stats> stats) {

    qpOASES::int_t nWSR = options_->qp_maxiter;

    if (!firstQPsolved_) {//if haven't solve any QP before then initialize the first QP
        set_solver_options();

//@{
//for debugging
//        H_qpOASES_->print("H_qp_oases");
//        A_qpOASES_->print("A_qpoases");
//        g_->print("g");
//        lbA_->print("LbA");
//        ubA_->print("ubA");
//        lb_->print("lb");
//        ub_->print("ub");
        //@}

        solver_->init(H_qpOASES_.get(), g_->values(), A_qpOASES_.get(), lb_->values(),
                      ub_->values(), lbA_->values(), ubA_->values(), nWSR);

        if (solver_->isSolved()) {
            firstQPsolved_ = true;
        }
        else {
            handle_error(QP, stats);
        }
    }
    else {
//for debugging
//@{
//                  H_qpOASES_->print("H_qp_oases");
//                  A_qpOASES_->print("A_qpoases");
//                  g_->print("g");
//                  lbA_->print("LbA");
//                  ubA_->print("ubA");
//                  lb_->print("lb");
//                  ub_->print("ub");
        //@}
        get_Matrix_change_status();
        if (new_QP_matrix_status_ == UNDEFINED) {
            assert(old_QP_matrix_status_ != UNDEFINED);
            if (old_QP_matrix_status_ == FIXED)
                solver_->hotstart(g_->values(), lb_->values(), ub_->values(),
                                  lbA_->values(),
                                  ubA_->values(), nWSR);
            else {
                solver_->hotstart(H_qpOASES_.get(), g_->values(), A_qpOASES_.get(),
                                  lb_->values(), ub_->values(), lbA_->values(),
                                  ubA_->values(), nWSR);
            }
        }
        else {
            if (new_QP_matrix_status_ == FIXED && old_QP_matrix_status_ == FIXED) {
                solver_->hotstart(g_->values(), lb_->values(), ub_->values(),
                                  lbA_->values(),
                                  ubA_->values(), nWSR);
            }
            else if (new_QP_matrix_status_ == VARIED &&
                     old_QP_matrix_status_ == VARIED) {
                solver_->hotstart(H_qpOASES_.get(), g_->values(), A_qpOASES_.get(),
                                  lb_->values(), ub_->values(), lbA_->values(),
                                  ubA_->values(), nWSR);
            }
            else if (new_QP_matrix_status_ != old_QP_matrix_status_) {
                qpOASES::Bounds tmp_bounds;
                solver_->getBounds(tmp_bounds);
                solver_->init(H_qpOASES_.get(), g_->values(), A_qpOASES_.get(),
                              lb_->values(),
                              ub_->values(), lbA_->values(), ubA_->values(), nWSR,0,         x_qp_->values(),y_qp_->values(),&tmp_bounds);
                new_QP_matrix_status_ = old_QP_matrix_status_ = UNDEFINED;

            }
        }
    }

    reset_flags();

    if(stats!=nullptr)
        stats->qp_iter_addValue((int) nWSR);

    if (!solver_->isSolved()) {
        handle_error(QP, stats);
    }
    solver_->getPrimalSolution(x_qp_->values());
    solver_->getDualSolution(y_qp_->values());

}


void qpOASESInterface::optimizeLP(shared_ptr<Stats> stats) {
    qpOASES::int_t nWSR = options_->lp_maxiter;//TODO modify it
    if (!firstQPsolved_) {
        set_solver_options();
        solver_->init(0, g_->values(), A_qpOASES_.get(), lb_->values(),
                      ub_->values(), lbA_->values(), ubA_->values(), nWSR);
        if (solver_->isSolved()) {
            firstQPsolved_ = true;
        }
        else {
            handle_error(LP, stats);
        }
    }
    else {
        get_Matrix_change_status();
        if (new_QP_matrix_status_ == UNDEFINED) {
            if (old_QP_matrix_status_ == FIXED)
                solver_->hotstart(g_->values(), lb_->values(), ub_->values(),
                                  lbA_->values(),
                                  ubA_->values(), nWSR);
            else {
                solver_->hotstart(0, g_->values(), A_qpOASES_.get(),
                                  lb_->values(), ub_->values(), lbA_->values(),
                                  ubA_->values(), nWSR);
            }
        }

        else {
            if (new_QP_matrix_status_ == FIXED && old_QP_matrix_status_ == FIXED) {
                solver_->hotstart(g_->values(), lb_->values(), ub_->values(),
                                  lbA_->values(),
                                  ubA_->values(), nWSR);
            }
            else if (new_QP_matrix_status_ == VARIED &&
                     old_QP_matrix_status_ == VARIED) {
                solver_->hotstart(0, g_->values(), A_qpOASES_.get(),
                                  lb_->values(), ub_->values(), lbA_->values(),
                                  ubA_->values(), nWSR);
            }
            else if (new_QP_matrix_status_ != old_QP_matrix_status_) {
                solver_->init(0, g_->values(), A_qpOASES_.get(), lb_->values(),
                              ub_->values(), lbA_->values(), ubA_->values(), nWSR);
                new_QP_matrix_status_ = old_QP_matrix_status_ = UNDEFINED;
            }
        }
        reset_flags();

        if (!solver_->isSolved()) {
            handle_error(LP, stats);
        }
    }
    //get primal and dual solutions
    if(stats!=nullptr)
        stats->qp_iter_addValue((int) nWSR);
    solver_->getPrimalSolution(x_qp_->values());
    solver_->getDualSolution(y_qp_->values());

}


/**
 * @brief get the pointer to the multipliers to the bounds constraints.
 */
double* qpOASESInterface::get_multipliers_bounds() {
    return y_qp_->values();
}


/**
 * @brief get the pointer to the multipliers to the regular constraints.
 */
double* qpOASESInterface::get_multipliers_constr() {
    return y_qp_->values()+nVar_QP_;
};


/**
 * @brief copy the optimal solution of the QP to the input pointer
 */
double* qpOASESInterface::get_optimal_solution() {
    return x_qp_->values();
}


/**
 *@brief get the objective value from the QP solvers
 *
 * @return the objective function value of the QP problem
 */


inline double qpOASESInterface::get_obj_value() {

    return (double) (solver_->getObjVal());
}




Exitflag qpOASESInterface::get_status() {
    qpOASES::QProblemStatus finalStatus = solver_->getStatus();

    if (solver_->isInfeasible()) {
        return   QPERROR_INFEASIBLE;
    }
    else if (solver_->isUnbounded()) {
        return  QPERROR_UNBOUNDED;
    }
    else if(solver_->isSolved()) {
        return QP_OPTIMAL;
    }
    else
        switch (finalStatus) {
        case qpOASES::QPS_NOTINITIALISED:
            return  QPERROR_NOTINITIALISED;
        case qpOASES::QPS_PREPARINGAUXILIARYQP:
            return  QPERROR_PREPARINGAUXILIARYQP;
        case qpOASES::QPS_AUXILIARYQPSOLVED:
            return  QPERROR_AUXILIARYQPSOLVED;
        case qpOASES::QPS_PERFORMINGHOMOTOPY:
            return  QPERROR_PERFORMINGHOMOTOPY;
        case qpOASES::QPS_HOMOTOPYQPSOLVED:
            return  QPERROR_HOMOTOPYQPSOLVED;
        }
}


//@}
void qpOASESInterface::set_lb(int location, double value) {
    if (firstQPsolved_ && !data_change_flags_.Update_bounds)
        data_change_flags_.Update_bounds = true;
    lb_->setValueAt(location, value);
}


void qpOASESInterface::set_ub(int location, double value) {

    if (firstQPsolved_ && !data_change_flags_.Update_bounds)
        data_change_flags_.Update_bounds = true;
    ub_->setValueAt(location, value);
}


void qpOASESInterface::set_lbA(int location, double value) {

    if (firstQPsolved_ && !data_change_flags_.Update_bounds)
        data_change_flags_.Update_bounds = true;
    lbA_->setValueAt(location, value);
}


void qpOASESInterface::set_ubA(int location, double value) {
    if (firstQPsolved_ && !data_change_flags_.Update_bounds)
        data_change_flags_.Update_bounds = true;
    ubA_->setValueAt(location, value);
}


void qpOASESInterface::set_g(int location, double value) {
    if (firstQPsolved_ && !data_change_flags_.Update_g)
        data_change_flags_.Update_g = true;
    g_->setValueAt(location, value);
}




void qpOASESInterface::set_H(shared_ptr<const SpTripletMat> rhs) {
    //@for debugging
    //@{
    //	H_->print("H");
    //	rhs->print("rhs");
    //@}

    if (firstQPsolved_ && !data_change_flags_.Update_H) {
        data_change_flags_.Update_H = true;
    }
    if(!H_->isinitialized()) {
        H_->setStructure(rhs);//TODO: move to somewhere else?
        H_qpOASES_ = std::make_shared<qpOASES::SymSparseMat>(H_->RowNum(),
                     H_->ColNum(),
                     H_->RowIndex(),
                     H_->ColIndex(),
                     H_->MatVal());
    }
    else {
        H_->setMatVal(rhs);
        H_qpOASES_->setVal(H_->MatVal());
    }
    H_qpOASES_->createDiagInfo();
}


void qpOASESInterface::set_A(shared_ptr<const SQPhotstart::SpTripletMat> rhs, IdentityInfo I_info) {
    if (firstQPsolved_ && !data_change_flags_.Update_A) {
        data_change_flags_.Update_A = true;
    }
    if(!A_->isinitialized()) {
        A_->setStructure(rhs, I_info);
        A_qpOASES_ = std::make_shared<qpOASES::SparseMatrix>(A_->RowNum(),
                     A_->ColNum(),
                     A_->RowIndex(),
                     A_->ColIndex(),
                     A_->MatVal());
    }
    else {
        A_->setMatVal(rhs, I_info);
        A_qpOASES_->setVal(A_->MatVal());
    }
}


void qpOASESInterface::set_ub(shared_ptr<const Vector> rhs) {

    if (firstQPsolved_ && !data_change_flags_.Update_bounds)
        data_change_flags_.Update_bounds = true;
    ub_->copy_vector(rhs->values());
}


void qpOASESInterface::set_lb(shared_ptr<const Vector> rhs) {

    if (firstQPsolved_ && !data_change_flags_.Update_bounds)
        data_change_flags_.Update_bounds = true;
    lb_->copy_vector(rhs->values());

}


void qpOASESInterface::set_lbA(shared_ptr<const Vector> rhs) {

    if (firstQPsolved_ && !data_change_flags_.Update_bounds)
        data_change_flags_.Update_bounds = true;
    lbA_->copy_vector(rhs->values());
}


void qpOASESInterface::set_ubA(shared_ptr<const Vector> rhs) {

    if (firstQPsolved_ && !data_change_flags_.Update_bounds)
        data_change_flags_.Update_bounds = true;
    ubA_->copy_vector(rhs->values());

}


void qpOASESInterface::set_g(shared_ptr<const Vector> rhs) {

    if (firstQPsolved_ && !data_change_flags_.Update_g)
        data_change_flags_.Update_g = true;
    g_->copy_vector(rhs->values());
}



void qpOASESInterface::reset_flags() {

    data_change_flags_.Update_A = false;
    data_change_flags_.Update_delta = false;
    data_change_flags_.Update_H = false;
    data_change_flags_.Update_penalty = false;
    data_change_flags_.Update_g = false;
    data_change_flags_.Update_bounds = false;
}

bool qpOASESInterface::test_optimality(ActiveType* W_c, ActiveType* W_b) {

    int i;
    //create local variables and set all violation values to be 0
    double primal_violation = 0.0;
    double dual_violation = 0.0;
    double compl_violation = 0.0;
    double statioanrity_violation = 0.0;
    shared_ptr<Vector> Ax = make_shared<Vector>(nConstr_QP_);


    if(W_c==NULL&& W_b==NULL) {
        W_c = new ActiveType[nConstr_QP_];
        W_b = new ActiveType[nVar_QP_];
    }
    get_working_set(W_c,W_b);

    /**-------------------------------------------------------**/
    /**                    primal feasibility                 **/
    /**-------------------------------------------------------**/
    for (i = 0; i < nVar_QP_; i++) {
        primal_violation += max(0.0, (lb_->values(i) - x_qp_->values(i)));
        primal_violation += -min(0.0, (ub_->values(i) - x_qp_->values(i)));
    }
    if(A_ != nullptr) {
        A_->times(x_qp_, Ax); //tmp_vec_nCon=A*x
        for (i = 0; i < nConstr_QP_; i++) {
            primal_violation += max(0.0, (lbA_->values(i) -Ax->values(i)));
            primal_violation += -min(0.0, (ubA_->values(i) -Ax->values(i)));
        }
    }

    /**-------------------------------------------------------**/
    /**                    dual feasibility                   **/
    /**-------------------------------------------------------**/
    for (i = 0; i < nVar_QP_; i++) {
        switch (W_b[i]) {
        case INACTIVE://the constraint is inactive, then the dual multiplier
            // should be 0
            dual_violation += fabs(y_qp_->values(i));
            break;
        case ACTIVE_BELOW://the constraint is active at the lower bound, so the
            // multiplier should be positive
            dual_violation += -min(0.0, y_qp_->values(i));
            break;
        case ACTIVE_ABOVE: //the contraint is active at the upper bounds, so the
            // multiplier should be negavie
            dual_violation += max(0.0, y_qp_->values(i));
            break;
        case ACTIVE_BOTH_SIDE:
            break;
        default:
            printf("failed in dual fea test, the working set  for qpOASES at "
                   "the bounds is %i", W_b[i]);
            THROW_EXCEPTION(INVALID_WORKING_SET, INVALID_WORKING_SET_MSG);
        }
    }
    if (A_ != nullptr) {
        for (i = 0; i < nConstr_QP_; i++) {
            switch (W_c[i]) {
            case INACTIVE://the constraint is inactive, then the dual multiplier
                // should be 0
                dual_violation += fabs(y_qp_->values(i+nVar_QP_));
                break;
            case ACTIVE_BELOW://the constraint is active at the lower bound, so the
                // multiplier should be positive
                dual_violation += -min(0.0, y_qp_->values(i+nVar_QP_));
                break;
            case ACTIVE_ABOVE: //the contraint is active at the upper bounds, so the
                // multiplier should be negavie
                dual_violation += max(0.0, y_qp_->values(i+nVar_QP_));
                break;
            case ACTIVE_BOTH_SIDE:
                break;
            default:
                printf("failed in dual fea test, the working set  for qpOASES at "
                       "the constr is %i", W_c[i]);
                THROW_EXCEPTION(INVALID_WORKING_SET, INVALID_WORKING_SET_MSG);
            }
        }
    }


    /**-------------------------------------------------------**/
    /**                   stationarity                        **/
    /**-------------------------------------------------------**/
    //calculate A'*y+lambda-(g+Hx)


    shared_ptr<Vector> stationary_gap = make_shared<Vector>(nVar_QP_);
//    A_->print("A");
//    H_->print("H");
//    x_qp_->print("x_qp");
//    lb_->print("lb");
//    ub_->print("ub");
//    y_qp_->print("y_qp");

    if (A_ != nullptr) {
        A_->transposed_times(y_qp_->values()+nVar_QP_, stationary_gap->values());
    }
    shared_ptr<Vector> Hx = make_shared<Vector>(nVar_QP_);
    H_->times(x_qp_, Hx);

    stationary_gap->add_vector(y_qp_->values());
    stationary_gap->subtract_vector(g_->values());
    stationary_gap->subtract_vector(Hx->values());
    statioanrity_violation = stationary_gap->getOneNorm();


    /**-------------------------------------------------------**/
    /**                    Complemtarity                      **/
    /**-------------------------------------------------------**/

    for (i = 0; i < nVar_QP_; i++) {
        switch (W_b[i]) {
        case INACTIVE: //constraint is inactive, multiplier should be 0
            compl_violation += abs(y_qp_->values(i));
            break;
        case ACTIVE_BELOW://the constraint is active at the lower bound
            compl_violation += abs(y_qp_->values(i) *
                                   (x_qp_->values(i) - lb_->values(i)));
            break;
        case ACTIVE_ABOVE: //the contraint is active at the upper bounds, so the
            // multiplier should be negavie
            compl_violation += abs(y_qp_->values(i) *
                                   (ub_->values(i) - x_qp_->values(i)));
            break;
        case ACTIVE_BOTH_SIDE:
            break;
        default:
            printf("failed in compl test, the working set  for qpOASES at "
                   "the "
                   "bounds is %i", W_b[i]);
            THROW_EXCEPTION(INVALID_WORKING_SET, INVALID_WORKING_SET_MSG);
        }
    }
    if (A_ != nullptr) {
        for (i = 0; i < nConstr_QP_; i++) {
            switch (W_c[i]) {
            case INACTIVE: //constraint is inactive, multiplier should be 0
                compl_violation += abs(y_qp_->values(i+nVar_QP_));
                break;
            case ACTIVE_BELOW://the constraint is active at the lower bound
                compl_violation += abs(y_qp_->values(i+ nVar_QP_) *
                                       (Ax->values(i) - lbA_->values(i)));
                break;
            case ACTIVE_ABOVE: //the contraint is active at the upper bounds, so the
                // multiplier should be negavie
                compl_violation += abs(y_qp_->values(i+ nVar_QP_) *
                                       (ubA_->values(i) - Ax->values(i)));
                break;
            case ACTIVE_BOTH_SIDE:
                break;
            default:
                printf("failed in compl test, the working set  for qpOASES at "
                       "the "
                       "constr is %i", W_c[i]);
                THROW_EXCEPTION(INVALID_WORKING_SET, INVALID_WORKING_SET_MSG);
            }
        }
    }

    qpOptimalStatus_.compl_violation = compl_violation;
    qpOptimalStatus_.stationarity_violation = statioanrity_violation;
    qpOptimalStatus_.dual_violation = dual_violation;
    qpOptimalStatus_.primal_violation = primal_violation;
    qpOptimalStatus_.KKT_error =
        compl_violation + statioanrity_violation + dual_violation + primal_violation;


    if(W_c==NULL&&W_b==NULL) {
        delete [] W_c;
        delete [] W_b;
    }

    if(qpOptimalStatus_.KKT_error>1.0e-6) {
//        printf("comp_violation %10e\n", compl_violation);
//        printf("stat_violation %10e\n", statioanrity_violation);
//        printf("prim_violation %10e\n", primal_violation);
//        printf("dual_violation %10e\n", dual_violation);
//        printf("KKT_error %10e\n", qpOptimalStatus_.KKT_error);
        return false;
    }


    return true;
}

void qpOASESInterface::handle_error(QPType qptype, shared_ptr<Stats> stats) {

    if (qptype == LP) {
        qpOASES::int_t nWSR = options_->lp_maxiter;//TODO modify it
        if (solver_->isInfeasible()) {
            shared_ptr<Vector> x_0 = make_shared<Vector>(nVar_QP_);
            shared_ptr<Vector> Ax = make_shared<Vector>(nConstr_QP_);
            x_0->copy_vector(x_qp_);
            A_->times(x_0,Ax);

            for(int i=0; i<nConstr_QP_; i++) {
                x_0->setValueAt(i+nVar_QP_-2*nConstr_QP_,max(0.0,lbA_->values(i)));
                x_0->setValueAt(i+nVar_QP_-nConstr_QP_,-min(0.0,ubA_->values(i)));
            }
            solver_->init(H_qpOASES_.get(), g_->values(), A_qpOASES_.get(),
                          lb_->values(), ub_->values(), lbA_->values(),
                          ubA_->values(), nWSR,NULL, x_0->values());
        }
        else {
            solver_->init(0, g_->values(), A_qpOASES_.get(),
                          lb_->values(), ub_->values(), lbA_->values(),
                          ubA_->values(), nWSR);

        }
        old_QP_matrix_status_ = new_QP_matrix_status_ = UNDEFINED;
        if(stats!=nullptr)
            stats->qp_iter_addValue((int) nWSR);

        if (!solver_->isSolved()) {
            THROW_EXCEPTION(LP_NOT_OPTIMAL,LP_NOT_OPTIMAL_MSG);
        }
    }
    else {
        qpOASES::int_t nWSR = options_->qp_maxiter;//TODO modify it
        if (solver_->isInfeasible()) {
            shared_ptr<Vector> x_0 = make_shared<Vector>(nVar_QP_);
            for(int i=0; i<nConstr_QP_; i++) {
                x_0->setValueAt(i+nVar_QP_-2*nConstr_QP_,max(0.0,lbA_->values(i)));
                x_0->setValueAt(i+nVar_QP_-nConstr_QP_,-min(0.0,ubA_->values(i)));
            }
            solver_->init(H_qpOASES_.get(), g_->values(), A_qpOASES_.get(),
                          lb_->values(), ub_->values(), lbA_->values(),
                          ubA_->values(), nWSR,NULL, x_0->values());
            //for debugging
            //@{
//            H_qpOASES_->print("H_qp_oases");
//            A_qpOASES_->print("A_qpoases");
//            g_->print("g");
//            lbA_->print("LbA");
//            ubA_->print("ubA");
//            lb_->print("lb");
//            ub_->print("ub");
//            x_0->print("x_0");
//            shared_ptr<Vector> Ax= make_shared<Vector>(nConstr_QP_);
//
//            A_->times(x_0, Ax);
//            Ax->print("Ax");

            //@}
        }
        else {
            solver_->init(H_qpOASES_.get(), g_->values(), A_qpOASES_.get(),
                          lb_->values(), ub_->values(), lbA_->values(),
                          ubA_->values(), nWSR);
        }
        old_QP_matrix_status_ = new_QP_matrix_status_ = UNDEFINED;
        if(stats!=nullptr)
            stats->qp_iter_addValue((int) nWSR);
        if (!solver_->isSolved()) {
            THROW_EXCEPTION(QP_NOT_OPTIMAL,QP_NOT_OPTIMAL_MSG);
        }
    }
}


void qpOASESInterface::set_solver_options() {

    qpOASES::Options qp_options;

    qp_options.setToReliable();
    //setup the printlevel of q
    switch (options_->qpPrintLevel) {
    case 0:
        qp_options.printLevel = qpOASES::PL_NONE;
        break;
    case 1:
        qp_options.printLevel = qpOASES::PL_TABULAR;
        break;
    case 2:
        qp_options.printLevel = qpOASES::PL_LOW;
        break;
    case 3:
        qp_options.printLevel = qpOASES::PL_MEDIUM;
        break;
    case 4:
        qp_options.printLevel = qpOASES::PL_HIGH;
        break;
    case -2:
        qp_options.printLevel = qpOASES::PL_DEBUG_ITER;
        break;
    }
    solver_->setOptions(qp_options);
}


void qpOASESInterface::WriteQPDataToFile(Ipopt::EJournalLevel level,
        Ipopt::EJournalCategory category,
        const string filename) {
#if DEBUG
#if PRINT_OUT_QP_WITH_ERROR
    jnlst_->DeleteAllJournals();
    Ipopt::SmartPtr<Ipopt::Journal> QPdata_jrnl= jnlst_->AddFileJournal("QPdata",
            "qpOASES"+filename,Ipopt::J_WARNING);
    QPdata_jrnl->SetAllPrintLevels(level);
    QPdata_jrnl->SetPrintLevel(category,level);

    lb_->write_to_file("lb",jnlst_,level,category,QPOASES);
    lbA_->write_to_file("lbA",jnlst_,level,category,QPOASES);
    ub_->write_to_file("ub",jnlst_,level,category,QPOASES);
    ubA_->write_to_file("ubA",jnlst_,level,category,QPOASES);

    g_->write_to_file("g",jnlst_,level,category,QPOASES);
    A_->write_to_file("A",jnlst_,level,category,QPOASES);
    H_->write_to_file("H",jnlst_,level,category,QPOASES);
    jnlst_->DeleteAllJournals();
#endif
#endif

}


void qpOASESInterface::get_Matrix_change_status() {

    if (old_QP_matrix_status_ == UNDEFINED) {
        old_QP_matrix_status_ = data_change_flags_.Update_A ||
                                data_change_flags_.Update_H ? VARIED
                                : FIXED;
    }
    else {
        if (new_QP_matrix_status_ != UNDEFINED)
            old_QP_matrix_status_ = new_QP_matrix_status_;
        new_QP_matrix_status_ = data_change_flags_.Update_A ||
                                data_change_flags_.Update_H ? VARIED
                                : FIXED;
    }


}

void qpOASESInterface::get_working_set(SQPhotstart::ActiveType* W_constr,
                                       SQPhotstart::ActiveType* W_bounds) {
    int* tmp_W_c = new int[nConstr_QP_];
    int* tmp_W_b = new int[nVar_QP_];


    assert(nConstr_QP_==solver_->getNC());
    assert(nVar_QP_==solver_->getNV());
    solver_->getWorkingSetConstraints(tmp_W_c);
    solver_->getWorkingSetBounds(tmp_W_b);

    for (int i = 0; i < nVar_QP_; i++) {
        switch((int)tmp_W_b[i]) {
        case 1:
            if(fabs(x_qp_->values(i)-lb_->values(i))<sqrt_m_eps)
                W_bounds[i] = ACTIVE_BOTH_SIDE;
            else
                W_bounds[i] = ACTIVE_ABOVE;
            break;
        case -1:
            if(fabs(x_qp_->values(i)-ub_->values(i))<sqrt_m_eps)
                W_bounds[i] = ACTIVE_BOTH_SIDE;
            else
                W_bounds[i] = ACTIVE_BELOW;

            break;
        case 0:
            W_bounds[i] = INACTIVE;
            break;
        default:
            printf("invalud workingset for qpoases1");
            THROW_EXCEPTION(INVALID_WORKING_SET,INVALID_WORKING_SET_MSG);
        }
    }
    auto Ax = make_shared<Vector>(nConstr_QP_);
    A_->times(x_qp_, Ax); //tmp_vec_nCon=A*x
    for (int i = 0; i < nConstr_QP_; i++) {
        switch((int)tmp_W_c[i]) {
        case 1:
            if(fabs(Ax->values(i)-lbA_->values(i)<sqrt_m_eps))
                W_constr[i] = ACTIVE_BOTH_SIDE;
            else
                W_constr[i] = ACTIVE_ABOVE;
            break;
        case -1:
            if(fabs(Ax->values(i)-ubA_->values(i)<sqrt_m_eps))
                W_constr[i] = ACTIVE_BOTH_SIDE;
            else
                W_constr[i] = ACTIVE_BELOW;
            break;
        case 0:
            W_constr[i] = INACTIVE;
            break;
        default:
            printf("invalud workingset for qpoases2;");
            THROW_EXCEPTION(INVALID_WORKING_SET,INVALID_WORKING_SET_MSG);
        }
    }
    delete[] tmp_W_b;
    delete[] tmp_W_c;
}

void qpOASESInterface::reset_constraints() {
    lb_->set_zeros();
    ub_->set_zeros();
    lbA_->set_zeros();
    ubA_->set_zeros();
}




}//SQPHOTSTART

