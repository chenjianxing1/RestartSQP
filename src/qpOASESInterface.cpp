/* Copyright (C) 2019
 * All Rights Reserved.
 *
 * Authors: Xinyi Luo
 * Date:    2019-07
 */
#include <sqphot/qpOASESInterface.hpp>



namespace SQPhotstart {


/**
 * @brief Constructor which also initializes the qpOASES SQProblem objects
 * @param nlp_index_info the struct that stores simple nlp dimension info
 * @param qptype  is the problem to be solved QP or LP or SOC?
 */
qpOASESInterface::qpOASESInterface(NLPInfo nlp_index_info, QPType qptype,
                                   shared_ptr<const Options> options,
                                   Ipopt::SmartPtr<Ipopt::Journalist> jnlst):
    nConstr_QP_(nlp_index_info.nCon),
    jnlst_(jnlst),
    nVar_QP_(nlp_index_info.nVar+2*nlp_index_info.nCon),
    options_(options)
{
    allocate_memory(nlp_index_info, qptype);
}


/**Default destructor*/
qpOASESInterface::~qpOASESInterface() = default;


/**
 * @brief Allocate memory for the class members
 * @param nlp_index_info  the struct that stores simple nlp dimension info
 * @param qptype is the problem to be solved QP or LP or SOC?
 * @return
 */
void qpOASESInterface::allocate_memory(NLPInfo nlp_index_info, QPType qptype) {
    lbA_ = make_shared<Vector>(nConstr_QP_);
    ubA_ = make_shared<Vector>(nConstr_QP_);
    lb_ = make_shared<Vector>(nVar_QP_);
    ub_ = make_shared<Vector>(nVar_QP_);
    g_ = make_shared<Vector>(nVar_QP_);
    A_ = make_shared<SpHbMat>(
             nlp_index_info.nnz_jac_g + 2 * nlp_index_info.nCon, nConstr_QP_, nVar_QP_,false);
    x_qp_ = make_shared<Vector>(nVar_QP_);
    y_qp_ = make_shared<Vector>(nConstr_QP_+nVar_QP_);

    if (qptype != LP) {
        H_ = make_shared<SpHbMat>(nVar_QP_, nVar_QP_, false);
    }

    //@{
//TODO: for debugging
    A_triplet_ = make_shared<SpTripletMat>(nlp_index_info
                                           .nnz_jac_g+2*nlp_index_info.nCon,nConstr_QP_,
                                           nVar_QP_,false);
    H_triplet_ = make_shared<SpTripletMat>(nlp_index_info.nnz_h_lag,nConstr_QP_,
                                           nVar_QP_,true);
    //@}
    //FIXME: the qpOASES does not accept any extra input
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
        //          H_qpOASES_->print("H_qp_oases");
        //          A_qpOASES_->print("A_qpoases");
        //          g_->print("g");
        //          lbA_->print("LbA");
        //          ubA_->print("ubA");
        //          lb_->print("lb");
        //          ub_->print("ub");
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
                solver_->init(H_qpOASES_.get(), g_->values(), A_qpOASES_.get(),
                              lb_->values(),
                              ub_->values(), lbA_->values(), ubA_->values(), nWSR);
                new_QP_matrix_status_ = old_QP_matrix_status_ = UNDEFINED;
            }
        }
    }

    reset_flags();

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
//    x_qp_->print("x_qp_");
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


void qpOASESInterface::set_H_structure(shared_ptr<const SpTripletMat> rhs) {

    H_->setStructure(rhs);//TODO: move to somewhere else?
    H_qpOASES_ = std::make_shared<qpOASES::SymSparseMat>(H_->RowNum(),
                 H_->ColNum(),
                 H_->RowIndex(),
                 H_->ColIndex(),
                 H_->MatVal());
}


void qpOASESInterface::set_H_values(shared_ptr<const SpTripletMat> rhs) {
    //@for debugging
    //@{
//	H_->print("H");
//	rhs->print("rhs");
    //@}
    if (firstQPsolved_ && !data_change_flags_.Update_H) {
        data_change_flags_.Update_H = true;
    }
    H_->setMatVal(rhs);

    H_qpOASES_->setVal(H_->MatVal());
    H_qpOASES_->createDiagInfo();
}


void qpOASESInterface::set_A_structure(shared_ptr<const SpTripletMat> rhs,
                                       Identity2Info I_info) {

    A_->setStructure(rhs, I_info);
    A_qpOASES_ = std::make_shared<qpOASES::SparseMatrix>(A_->RowNum(),
                 A_->ColNum(),
                 A_->RowIndex(),
                 A_->ColIndex(),
                 A_->MatVal());
}


void qpOASESInterface::set_A_values(

    shared_ptr<const SQPhotstart::SpTripletMat> rhs, Identity2Info I_info) {

    if (firstQPsolved_ && !data_change_flags_.Update_A) {
        data_change_flags_.Update_A = true;
    }
    A_->setMatVal(rhs, I_info);
    A_qpOASES_->setVal(A_->MatVal());
#if DEBUG
#if COMPARE_QP_SOLVER
    A_triplet_->convert2Triplet(A_);
#endif
#endif
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
            W_bounds[i] = ACTIVE_ABOVE;
            break;
        case -1:
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
    for (int i = 0; i < nConstr_QP_; i++) {
        switch((int)tmp_W_c[i]) {
        case 1:
            W_constr[i] = ACTIVE_ABOVE;
            break;
        case -1:
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

}

}//SQPHOTSTART
