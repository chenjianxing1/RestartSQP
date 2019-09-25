/* Copyright (C) 2019
 *
 * Authors: Xinyi Luo
 * Date:    2019-08-15
 */

#include <sqphot/QOREInterface.hpp>
using namespace std;
namespace SQPhotstart {
/**
 * @brief Constructor
 * @param nlp_info the index information for NLP
 * @param qptype QP or LP
 * @param options object stored user-defined parameter values
 * @param jnlst Ipopt Jourlist object, for printing out log files
 */
QOREInterface::QOREInterface(Index_info nlp_info,
                             QPType qptype,
                             shared_ptr<const Options> options,
                             Ipopt::SmartPtr<Ipopt::Journalist> jnlst) :
    jnlst_(jnlst),
    firstQPsolved_(false),
    solver_(0),
    nVar_QP_(nlp_info.nVar + nlp_info.nCon * 2),
    nConstr_QP_(nlp_info.nCon) {
    qpiter_[0] = 0;
    allocate_memory(nlp_info, qptype);
    set_solver_options(options);
}

/**
 * @brief Destructor
 */
QOREInterface::~QOREInterface() {
    delete[] working_set_;
    QPFree(&solver_);
}

/**
 * @brief Solve a regular QP with given data and options.
 */

void QOREInterface::optimizeQP(shared_ptr<Stats> stats) {
    /**-------------------------------------------------------**/
    /**                   Set Data and Optimize QP            **/
    /**-------------------------------------------------------**/
    rv_ = QPSetData(solver_, nVar_QP_, nConstr_QP_, A_->RowIndex(), A_->ColIndex(),
                    A_->MatVal(), H_->RowIndex(), H_->ColIndex(), H_->MatVal());

    assert(rv_ == QPSOLVER_OK);
    rv_ = QPOptimize(solver_, lb_->values(), ub_->values(), g_->values(), 0, 0);//
    firstQPsolved_ = true;
    assert(rv_ == QPSOLVER_OK);
    rv_ = QPGetInt(solver_, "status", &status_);
    assert(rv_ == QPSOLVER_OK);
    handle_error(QP);
    if (status_ != QPSOLVER_OPTIMAL) {
        THROW_EXCEPTION(QP_NOT_OPTIMAL,QP_NOT_OPTIMAL_MSG);
#if DEBUG
#if PRINT_OUT_QP_WITH_ERROR
        WriteQPDataToFile(jnlst_,Ipopt::J_LAST_LEVEL, Ipopt::J_USER1);
#endif
#endif
    }

    /**-------------------------------------------------------**/
    /**               Get Primal and Dual Solution            **/
    /**-------------------------------------------------------**/
    rv_ = QPGetDblVector(solver_, "primalsol", x_qp_->values());
    assert(rv_ == QPSOLVER_OK);
    rv_ = QPGetDblVector(solver_, "dualsol", y_qp_->values());
    assert(rv_ == QPSOLVER_OK);

    /**-------------------------------------------------------**/
    /**                     Update Stats                      **/
    /**-------------------------------------------------------**/
    QPGetInt(solver_, "itercount", qpiter_);
    stats->qp_iter_addValue(qpiter_[0]);

#if DEBUG
#if CHECK_QP_INFEASIBILITY
    A_->print("A_");
    H_->print("H_");
    lb_->print("lb_");
    ub_->print("ub_");
    g_->print("g_");
#endif
#endif
}


/**
 * @brief Solve a regular LP with given data and options
 */
void QOREInterface::optimizeLP(shared_ptr<Stats> stats) {
    /**-------------------------------------------------------**/
    /**                   Set Data and Optimize LP            **/
    /**-------------------------------------------------------**/
    rv_ = QPSetData(solver_, nVar_QP_, nConstr_QP_, A_->RowIndex(), A_->ColIndex(),
                    A_->MatVal(), NULL, NULL, NULL);
    assert(rv_ == QPSOLVER_OK);
    rv_ = QPOptimize(solver_, lb_->values(), ub_->values(), g_->values(), 0, 0);//
    assert(rv_ == QPSOLVER_OK);

    rv_ = QPGetInt(solver_, "status", &status_);
    assert(rv_ == QPSOLVER_OK);
    handle_error(LP);
    if (status_ != QPSOLVER_OPTIMAL)
        THROW_EXCEPTION(LP_NOT_OPTIMAL, LP_NOT_OPTIMAL_MSG);


    /**-------------------------------------------------------**/
    /**               Get Primal and Dual Solution            **/
    /**-------------------------------------------------------**/
    rv_ = QPGetDblVector(solver_, "primalsol", x_qp_->values());
    assert(rv_ == QPSOLVER_OK);
    rv_ = QPGetDblVector(solver_, "dualsol", y_qp_->values());
    assert(rv_ == QPSOLVER_OK);
    /**-------------------------------------------------------**/
    /**                     Update Stats                      **/
    /**-------------------------------------------------------**/
    QPGetInt(solver_, "itercount", qpiter_);
    stats->qp_iter_addValue(qpiter_[0]);
}


/**
 * @brief Allocate memory for the class members
 * @param nlp_index_info  the struct that stores simple nlp dimension info
 * @param qptype is the problem to be solved QP or LP?
 */

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
    working_set_  =  new int[nConstr_QP_+nVar_QP_];
    if (qptype != LP) {
        rv_ = QPNew(&solver_, nVar_QP_, nConstr_QP_, nnz_g_QP,
                    nlp_info.nnz_h_lag);
        H_ = make_shared<SpHbMat>(nVar_QP_, nVar_QP_, true,true);
    } else {
        //if we are solving an LP, the number of nonzero in the Hessian is 0
        rv_ = QPNew(&solver_, nVar_QP_, nConstr_QP_, nnz_g_QP, 0);
    }


    assert(rv_ == QPSOLVER_OK);


    //TODO: for debugging use only
    A_triplet_ = make_shared<SpTripletMat>(nnz_g_QP,nConstr_QP_,
                                           nVar_QP_);
    H_triplet_ = make_shared<SpTripletMat>(nlp_info.nnz_h_lag,nConstr_QP_,
                                           nVar_QP_,true);
}



/**@name Getters*/
//@{

double QOREInterface::get_obj_value() {
    shared_ptr<Vector> Hx = make_shared<Vector>(nVar_QP_);
    H_->times(x_qp_,Hx);
    return (Hx->times(x_qp_)*0.5+g_->times(x_qp_));
}


Exitflag QOREInterface::get_status() {
    switch(status_) {
    case QPSOLVER_ITER_LIMIT:
        return QPERROR_EXCEED_MAX_ITER;
    case QPSOLVER_INFEASIBLE:
        return QPERROR_INFEASIBLE;
    case QPSOLVER_UNBOUNDED:
        return QPERROR_UNBOUNDED;
    default:
        return QPERROR_UNKNOWN;

    }
}

void QOREInterface::get_working_set(ActiveType* W_constr, ActiveType* W_bounds) {
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
                W_constr[i-nVar_QP_] = ACTIVE_ABOVE;
                break;
            case -1:
                W_constr[i-nVar_QP_] = ACTIVE_BELOW;
                break;
            case 0:
                W_constr[i-nVar_QP_] = INACTIVE;
                break;
            default:
                printf("invalud workingset for qore2, the working set is %i;",
                       W_constr[i]);
                THROW_EXCEPTION(INVALID_WORKING_SET,INVALID_WORKING_SET_MSG);
            }
        }
    }

}



//@}



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

    lb_->write_to_file("lb",jnlst,level,category,QORE);
    ub_->write_to_file("ub",jnlst,level,category,QORE);
    g_->write_to_file("g",jnlst,level,category,QORE);
    A_->write_to_file("A",jnlst,level,category,QORE);
    H_->write_to_file("H",jnlst,level,category,QORE);

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

void QOREInterface::handle_error(QPType qptype) {
    switch(status_) {
    case QPSOLVER_OPTIMAL:
        //do nothing here
        break;
    case QPSOLVER_INFEASIBLE:
        shared_ptr<Vector> x_0 = make_shared<Vector>(nVar_QP_);
        shared_ptr<Vector> Ax = make_shared<Vector>(nConstr_QP_);
        x_0->copy_vector(x_qp_);
        A_->times(x_0,Ax);
//            Ax->print("Ax");
//            lb_->print("lb_");
//            ub_->print("ub_");
        //setup the slack variables to satisfy the bound constraints
        for(int i=0; i<nConstr_QP_; i++) {
            x_0->setValueAt(i+nVar_QP_-2*nConstr_QP_,max(0.0,lb_->values(nVar_QP_+i)));
            x_0->setValueAt(i+nVar_QP_-nConstr_QP_,-min(0.0,ub_->values(nVar_QP_+i)));
        }

//            x_0->print("x_0");
//
//            A_->times(x_0,Ax);
//
//            Ax->print("Ax");
//            A_->print("A");


        rv_ = QPOptimize(solver_, lb_->values(), ub_->values(), g_->values(),
                         x_0->values(), NULL);//
        assert(rv_ == QPSOLVER_OK);
        break;
    }

}


void QOREInterface::set_solver_options(shared_ptr<const Options> options) {


    if(options->qpPrintLevel==0) {
        //does not print anything
        QPSetInt(solver_, "prtfreq", -1);
    }
}

}//SQP_HOTSTART

