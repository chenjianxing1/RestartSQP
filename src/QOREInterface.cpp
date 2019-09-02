/* Copyright (C) 2019
 *
 * Authors: Xinyi Luo
 * Date:    2019-07
 * Date:    2019-08-15
 */

#include <sqphot/QOREInterface.hpp>

namespace SQPhotstart {

QOREInterface::QOREInterface(Index_info nlp_info,
                             QPType qptype,
                             shared_ptr<const Options> options,
                             Ipopt::SmartPtr<Ipopt::Journalist> jnlst) :
    jnlst_(jnlst),
    firstQPsolved_(false),
    solver_(0),
    nVar_QP_(nlp_info.nVar + nlp_info.nCon * 2),
    nConstr_QP_(nlp_info.nCon) {
    allocate_memory(nlp_info, qptype);


}


QOREInterface::~QOREInterface() {
#if DEBUG
#if COMPARE_QP_SOLVER
    delete[] working_set_;
#endif
#endif
    QPFree(&solver_);
}

void QOREInterface::optimizeQP(shared_ptr<Stats> stats) {
    rv_ = QPSetData(solver_, nVar_QP_, nConstr_QP_, A_->RowIndex(), A_->ColIndex(),
                    A_->MatVal(), H_->RowIndex(), H_->ColIndex(), H_->MatVal());
#if DEBUG
#if CHECK_QP_INFEASIBILITY
    A_->print();
    H_->print();
    lb_->print("lb_");
    ub_->print("ub_");
    g_->print("g_");
#endif
#endif
    assert(rv_ == QPSOLVER_OK);
    QPSetInt(solver_, "prtfreq", -1);
    if (!firstQPsolved_) {
        rv_ = QPOptimize(solver_, lb_->values(), ub_->values(), g_->values(), 0, 0);//
        firstQPsolved_ = true;
        assert(rv_ == QPSOLVER_OK);
        rv_ = QPGetInt(solver_, "status", &status_);
        if (status_ != QPSOLVER_OPTIMAL) {

#if DEBUG
#if PRINT_OUT_QP_WITH_ERROR
            WriteQPDataToFile(jnlst_,Ipopt::J_LAST_LEVEL, Ipopt::J_USER1);
#endif
#endif
            THROW_EXCEPTION(QP_NOT_OPTIMAL,
                            "the QP problem didn't solved to optimality\n")
        }
    }
    else {
        //TODO:monitor the data changes.
        rv_ = QPOptimize(solver_, lb_->values(), ub_->values(), g_->values(),
                         0, 0);//TODO: does not update primal-dual sol here
        assert(rv_ == QPSOLVER_OK);
        rv_ = QPGetInt(solver_, "status", &status_);
        if (status_ != QPSOLVER_OPTIMAL) {
#if DEBUG
#if PRINT_OUT_QP_WITH_ERROR
            WriteQPDataToFile(jnlst_,Ipopt::J_LAST_LEVEL, Ipopt::J_USER1);
#endif
#endif
            THROW_EXCEPTION(QP_NOT_OPTIMAL,
                            "the QP problem didn't solved to optimality\n")
        }


    }
    assert(rv_ == QPSOLVER_OK);

}


void QOREInterface::optimizeLP(shared_ptr<Stats> stats) {
    rv_ = QPSetData(solver_, nVar_QP_, nConstr_QP_, A_->RowIndex(), A_->ColIndex(),
                    +                   A_->MatVal(), NULL, NULL, NULL);
    assert(rv_ == QPSOLVER_OK);
    //set the print level.
    //TODO: adjust it based on user option
    QPSetInt(solver_, "prtfreq", -1);

    if (!firstQPsolved_) {
        rv_ = QPOptimize(solver_, lb_->values(), ub_->values(), g_->values(), 0, 0);//
        assert(rv_ == QPSOLVER_OK);
        rv_ = QPGetInt(solver_, "status", &status_);
        if (status_ != QPSOLVER_OPTIMAL)
            THROW_EXCEPTION(LP_NOT_OPTIMAL,
                            "the LP problem didn't solved to optimality\n")

            firstQPsolved_ = true;
    } else {
        rv_ = QPOptimize(solver_, lb_->values(), ub_->values(), g_->values(),
                         0, 0);//TODO: does not update primal-dual sol here
        rv_ = QPGetInt(solver_, "status", &status_);
        if (status_ != QPSOLVER_OPTIMAL)
            THROW_EXCEPTION(LP_NOT_OPTIMAL,
                            "the LP problem didn't solved to optimality\n")
        }
    assert(rv_ == QPSOLVER_OK);
}


void QOREInterface::allocate_memory(Index_info nlp_info, QPType qptype) {

    int nnz_g_QP = nlp_info.nnz_jac_g +
                   2 * nlp_info.nCon;//number of nonzero variables in jacobian
    //The Jacobian has the structure [J I -I], so it will contains extra 2*number_constr nonzero elements
    lb_ = make_shared<Vector>(nVar_QP_ + nConstr_QP_);
    ub_ = make_shared<Vector>(nVar_QP_ + nConstr_QP_);
    g_ = make_shared<Vector>(nVar_QP_);
    A_ = make_shared<SpHbMat>(nnz_g_QP, nConstr_QP_, nVar_QP_,true);
    x_qp_ = make_shared<Vector>(nVar_QP_ + nConstr_QP_);
    y_qp_ = make_shared<Vector>(nConstr_QP_ + nVar_QP_);
    if (qptype != LP) {
        rv_ = QPNew(&solver_, nVar_QP_, nConstr_QP_, nnz_g_QP,
                    nlp_info.nnz_h_lag);
        H_ = make_shared<SpHbMat>(nVar_QP_, nVar_QP_, true,true);
    } else {
        //if we are solving an LP, the number of nonzero in the Hessian is 0

        rv_ = QPNew(&solver_, nVar_QP_, nConstr_QP_, nnz_g_QP, 0);
    }

#if DEBUG
#if COMPARE_QP_SOLVER
    working_set_  =  new int[nConstr_QP_+nVar_QP_];
    A_triplet_ = make_shared<SpTripletMat>(nnz_g_QP,nConstr_QP_,
                                           nVar_QP_);
    H_triplet_ = make_shared<SpTripletMat>(nlp_info.nnz_h_lag,nConstr_QP_,
                                           nVar_QP_,true);
#endif
#endif
    assert(rv_ == QPSOLVER_OK);
}

double* QOREInterface::get_optimal_solution() {
    rv_ = QPGetDblVector(solver_, "primalsol", x_qp_->values());
    assert(rv_ == QPSOLVER_OK);
#if DEBUG
#if CHECK_QP_INFEASIBILITY
    x_qp_->print("x_qp_");
#endif
#endif
    return x_qp_->values();
}

double QOREInterface::get_obj_value() {
    //FIXME: QORE does not have existing emthod to return the qp obj, now it is calculated
    // in Algorithm class for QOREInterface...
}


double* QOREInterface::get_multipliers() {
    rv_ = QPGetDblVector(solver_, "dualsol", y_qp_->values());
    assert(rv_ == QPSOLVER_OK);

#if DEBUG
#if CHECK_QP_INFEASIBILITY
    y_qp_->print("y_qp_");
#endif
#endif

    return y_qp_->values();
}

void QOREInterface::set_ub(int location, double value) {
    value = value < INF ? value : INF;
    ub_->setValueAt(location, value);

}

void QOREInterface::set_g(int location, double value) {
    value = value < INF ? value : INF;
    g_->setValueAt(location, value);
}

void QOREInterface::set_lb(int location, double value) {
    value = value > -INF ? value : -INF;
    lb_->setValueAt(location, value);
}

void QOREInterface::set_A_structure(shared_ptr<const SpTripletMat> rhs,
                                    Identity2Info I_info) {

#if DEBUG
#if GET_QP_INTERFACE_MEMBERS or COMPARE_QP_SOLVER
    A_triplet_->copy(rhs, false);
#endif
#endif
    A_->setStructure(rhs, I_info);
}

void
QOREInterface::set_A_values(shared_ptr<const SpTripletMat> rhs,
                            Identity2Info I_info) {

    A_->setMatVal(rhs, I_info);
//    A_->print("A",jnlst_);
}

QPReturnType QOREInterface::get_status() {
    switch(status_) {
    case QPSOLVER_ITER_LIMIT:
        return QP_EXCEED_MAX_ITER;
    case QPSOLVER_OPTIMAL:
        return QP_OPTIMAL;
    case QPSOLVER_INFEASIBLE:
        return QP_INFEASIBLE;
    case QPSOLVER_UNBOUNDED:
        return QP_UNBOUNDED;
    default:
        return QP_UNKNOWN_ERROR;

    }
}



void QOREInterface::WriteQPDataToFile(Ipopt::SmartPtr<Ipopt::Journalist> jnlst,
                                      Ipopt::EJournalLevel level,
                                      Ipopt::EJournalCategory category) {


#if DEBUG

    Ipopt::SmartPtr<Ipopt::Journal> QPdata_jrnl = jnlst->GetJournal("QPdata");
    if (IsNull(QPdata_jrnl)) {
        QPdata_jrnl= jnlst->AddFileJournal("QPdata", "qpdata.out",
                                           Ipopt::J_WARNING);
    }
    QPdata_jrnl->SetAllPrintLevels(level);
    QPdata_jrnl->SetPrintLevel(category,level);
#if PRINT_OUT_QP_WITH_ERROR
#if PRINT_QP_IN_CPP
    jnlst->Printf(level, category, "#include <stdio.h>\n"
                  "#include <assert.h>\n"
                  "#include <stdlib.h>\n"
                  "#include <string.h>\n"
                  "#include <math.h>\n"
                  "#include <matrixconversion.h>\n"
                  "#include <qpsolver.h>\n"
                  "\n"
                  "int main(){\n"
                  "#define NV %i\n#define NC %i\n", nVar_QP_, nConstr_QP_);
#else
    jnlst->Printf(level, category, "%d\n", nVar_QP_);
    jnlst->Printf(level, category, "%d\n", nCon_QP_);
    jnlst->Printf(level, category, "%d\n", A_->EntryNum();
                  jnlst->Printf(level, category, "%d\n", H_->EntryNum();
#endif

    lb_->write_to_file("lb",jnlst,level,category,QORE_QP);
    ub_->write_to_file("ub",jnlst,level,category,QORE_QP);
    g_->write_to_file("g",jnlst,level,category,QORE_QP);
    A_->write_to_file("A",jnlst,level,category,QORE_QP);
    H_->write_to_file("H",jnlst,level,category,QORE_QP);

#if PRINT_QP_IN_CPP
    jnlst->Printf(level,category,"QoreProblem * qp = 0;\n");
    jnlst->Printf(level,category,"qp_int rv = QPNew( &qp, NV, NC, %i, %i );\n",
                  A_->EntryNum(), H_->EntryNum());
    jnlst->Printf(level,category,"assert( rv == QPSOLVER_OK );\n");
    jnlst->Printf(level,category,"assert( qp!= 0 );\n");
    jnlst->Printf(level,category,"QPSetInt(qp, \"prtfreq\", 0); \n");
    // pass problem data to solver
    jnlst->Printf(level, category, "rv = QPSetData( qp, NV, NC, A_jc, A_ir, A_val, H_jc,"
                  " H_ir, H_val );\n");

    jnlst->Printf(level,category,"assert( rv == QPSOLVER_OK );\n");
    // solve first QP
    jnlst->Printf(level, category,"rv = QPOptimize( qp, lb, ub, g, 0, 0 );\n");

    jnlst->Printf(level,category,"qp_int status;\n");
    jnlst->Printf(level,category,"QPGetInt(qp, \"status\", &status);\n");
    jnlst->Printf(level,category,"if(status!=QPSOLVER_OPTIMAL){\n)");
    jnlst->Printf(level,category,"    printf(\"Warning! The QP is not solved to optimality!\");\n");
    jnlst->Printf(level,category,"}\n");
    jnlst->Printf(level,category,"assert( rv == QPSOLVER_OK );\n");
    // get and print primal solution

    jnlst->Printf(level, category, "QPFree(&qp);\n"
                  "\n return 0; \n"
                  "}\n");
#endif
#endif
    jnlst->DeleteAllJournals();
#endif

}


void QOREInterface::GetWorkingSet(ActiveType* W_constr, ActiveType* W_bounds) {

    QPGetIntVector(solver_,"workingset",working_set_);
    for(int i=0; i<nConstr_QP_+nVar_QP_; i++) {
        if(i<nVar_QP_) {
            switch (working_set_[i]) {
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
                printf("invalud workingset for qore1;");
                THROW_EXCEPTION(INVALID_WORKING_SET,INVALID_WORKING_SET_MSG);
            }
        }
        else {
            switch (working_set_[i]) {
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
                printf("invalud workingset for qore2, the working set is %i;",
                        W_constr[i]);
                THROW_EXCEPTION(INVALID_WORKING_SET,INVALID_WORKING_SET_MSG);
            }
        }
    }

}

void QOREInterface::set_H_structure(shared_ptr<const SpTripletMat> rhs) {
#if DEBUG
#if GET_QP_INTERFACE_MEMBERS or COMPARE_QP_SOLVER
    H_triplet_->copy(rhs, false);
#endif
#endif
    H_->setStructure(rhs);
}

void QOREInterface::set_H_values(shared_ptr<const SpTripletMat> rhs) {
    H_->setMatVal(rhs);
//        H_->print("H",jnlst_);
}


}//SQP_HOTSTART

