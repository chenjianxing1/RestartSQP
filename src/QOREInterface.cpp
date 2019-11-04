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
QOREInterface::QOREInterface(NLPInfo nlp_info,
                             QPType qptype,
                             shared_ptr<const Options> options,
                             Ipopt::SmartPtr<Ipopt::Journalist> jnlst) :
    jnlst_(jnlst),
    firstQPsolved_(false),
    solver_(0) {
#if NEW_FORMULATION
    nConstr_QP_ = nlp_info.nCon+nlp_info.nVar;
    nVar_QP_ = nlp_info.nVar*3+2*nlp_info.nCon;
#else
    nConstr_QP_ = nlp_info.nCon;
    nVar_QP_ = nlp_info.nVar+2*nlp_info.nCon;
#endif
    qpiter_[0] = 0;
    allocate_memory(nlp_info, qptype);
    set_solver_options(options);
}

QOREInterface::QOREInterface(shared_ptr<SpHbMat> H,
                             shared_ptr<SpHbMat> A,
                             shared_ptr<Vector> g,
                             shared_ptr<Vector> lb,
                             shared_ptr<Vector> ub,
                             shared_ptr<const Options> options):
    A_(A),
    H_(H),
    g_(g),
    lb_(lb),
    ub_(ub),
    nVar_QP_(A->ColNum()),
    nConstr_QP_(A->RowNum()),
    firstQPsolved_(false),
    solver_(0) {
    x_qp_= make_shared<Vector>(nVar_QP_+ nConstr_QP_);
    y_qp_= make_shared<Vector>(nVar_QP_+ nConstr_QP_);
    working_set_ = new int[nVar_QP_+ nConstr_QP_]();
    qpiter_[0] = 0;
    rv_ = QPNew(&solver_, nVar_QP_, nConstr_QP_, A->EntryNum(),H->EntryNum());
    assert(rv_ == QPSOLVER_OK);
    if(options!=nullptr)
        set_solver_options(options);
}




/**
 * @brief Destructor
 */
QOREInterface::~QOREInterface() {
    delete[] working_set_;
    working_set_ = NULL;
    QPFree(&solver_);
}

/**
 * @brief Solve a regular QP with given data and options.
 */

void QOREInterface::optimizeQP(shared_ptr<Stats> stats) {


    //@{ For debug
//        printf("QP_var is %d, QP_con is %d\n\n",nVar_QP_,nConstr_QP_);
//        A_->print("A");
//        H_->print_full("H");

    //@}
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
    handle_error(QP,stats);
    if (status_ != QPSOLVER_OPTIMAL) {
        THROW_EXCEPTION(QP_NOT_OPTIMAL,QP_NOT_OPTIMAL_MSG);
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
    if(stats!=nullptr)
        stats->qp_iter_addValue(qpiter_[0]);
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
    handle_error(LP,stats);
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

void QOREInterface::allocate_memory(NLPInfo nlp_info, QPType qptype) {

#if NEW_FORMULATION
    int nnz_g_QP = nlp_info.nnz_jac_g+2*nlp_info.nCon+2*nlp_info.nVar;
#else
    int nnz_g_QP = nlp_info.nnz_jac_g +
                   2 * nlp_info.nCon;//number of nonzero variables in jacobian
    //The Jacobian has the structure [J I -I], so it will contains extra 2*number_constr
    //nonzero elements
#endif

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
        H_ = make_shared<SpHbMat>(nVar_QP_, nVar_QP_, true);
    } else {
        //if we are solving an LP, the number of nonzero in the Hessian is 0
        rv_ = QPNew(&solver_, nVar_QP_, nConstr_QP_, nnz_g_QP, 0);
    }

    assert(rv_ == QPSOLVER_OK);
}


bool QOREInterface::test_optimality(ActiveType* W_c, ActiveType* W_b) {
    int i;
    //create local variables and set all violation values to be 0
    double primal_violation = 0.0;
    double dual_violation = 0.0;
    double compl_violation = 0.0;
    double statioanrity_violation = 0.0;

    shared_ptr<Vector> Ax = make_shared<Vector>(nConstr_QP_);
    shared_ptr<Vector> stationary_gap = make_shared<Vector>(nVar_QP_);


    if(W_c==NULL&&W_b==NULL) {
        W_c = new ActiveType[nConstr_QP_];
        W_b = new ActiveType[nVar_QP_];
    }

    get_working_set(W_c,W_b);



    /**-------------------------------------------------------**/
    /**                    primal feasibility                 **/
    /**-------------------------------------------------------**/
    for (i = 0; i < nVar_QP_+nConstr_QP_; i++) {
        primal_violation += max(0.0, (lb_->values(i) - x_qp_->values(i)));
        primal_violation += -min(0.0, (ub_->values(i) - x_qp_->values(i)));
//        printf("primal_violation[%d]=%f\n",i,primal_violation);
    }
    /**-------------------------------------------------------**/
    /**                    dual feasibility                   **/
    /**-------------------------------------------------------**/


    for (i = 0; i < nVar_QP_; i++) {
        switch (W_b[i]) {
        case INACTIVE://the constraint is inactive, then the dual multiplier
            // should be 0
            //printf("0\n");
            dual_violation += fabs(y_qp_->values(i));
            break;
        case ACTIVE_BELOW://the constraint is active at the lower bound, so the
            // multiplier should be positive
//            printf("-1\n");
            dual_violation += -min(0.0, y_qp_->values(i));
            break;
        case ACTIVE_ABOVE: //the contraint is active at the upper bounds, so the
            // multiplier should be negavie
//            printf("1\n");
            dual_violation += max(0.0, y_qp_->values(i));
            break;
        case ACTIVE_BOTH_SIDE:
//           printf("99\n");
            break;
        default:
            printf("failed in dual fea test, the working set  for qore at the "
                   "bounds is %i", W_b[i]);
            THROW_EXCEPTION(INVALID_WORKING_SET, INVALID_WORKING_SET_MSG);
        }
    }

    if (A_ != nullptr) {
        for (i = 0; i < nConstr_QP_; i++) {
            switch (W_c[i]) {
            case INACTIVE://the constraint is inactive, then the dual multiplier
                // should be 0
                dual_violation += fabs(y_qp_->values(i+ nVar_QP_));
                //                  printf("0\n");
                break;
            case ACTIVE_BELOW://the constraint is active at the lower bound, so the
                // multiplier should be positive
                dual_violation += -min(0.0, y_qp_->values(i+ nVar_QP_));
                //                 printf("-1\n");
                break;
            case ACTIVE_ABOVE: //the contraint is active at the upper bounds, so the
                // multiplier should be negavie
                dual_violation += max(0.0, y_qp_->values(i+ nVar_QP_));
                //                printf("1\n");
                break;
            case ACTIVE_BOTH_SIDE:
                //       printf("99\n");
                break;
            default:
                printf("failed in dual fea test, the working set  for qore at the "
                       "constr is %i", W_c[i]);
                THROW_EXCEPTION(INVALID_WORKING_SET, INVALID_WORKING_SET_MSG);
            }
        }
    }
    /**-------------------------------------------------------**/
    /**                   stationarity                        **/
    /**-------------------------------------------------------**/
    //calculate A'*y+lambda-g-Hx
    A_->transposed_times(y_qp_->values()+nVar_QP_, stationary_gap->values());

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
            printf("failed in compl test, the working set  for qore at "
                   "the bounds is %i", W_b[i]);
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
                compl_violation += abs(y_qp_->values(i+nVar_QP_) *
                                       (x_qp_->values(i+nVar_QP_) -
                                        lb_->values(i + nVar_QP_)));
                break;
            case ACTIVE_ABOVE: //the contraint is active at the upper bounds, so the
                // multiplier should be negavie
                compl_violation += abs(y_qp_->values(i+nVar_QP_) *
                                       (ub_->values(i + nVar_QP_) -
                                        x_qp_->values(i+nVar_QP_)));
                break;
            case ACTIVE_BOTH_SIDE:
                break;
            default:
                printf("failed in compl test, the working set  for qpOASES at "
                       "the constr is %i", W_c[i]);
                THROW_EXCEPTION(INVALID_WORKING_SET, INVALID_WORKING_SET_MSG);
            }
        }
    }

    /**-------------------------------------------------------**/
    /**             Decide if x_k is optimal                  **/
    /**-------------------------------------------------------**/

    qpOptimalStatus_.compl_violation = compl_violation;
    qpOptimalStatus_.stationarity_violation = statioanrity_violation;
    qpOptimalStatus_.dual_violation = dual_violation;
    qpOptimalStatus_.primal_violation = primal_violation;

    qpOptimalStatus_.KKT_error =
        compl_violation + statioanrity_violation + dual_violation + primal_violation;

    if(W_c==NULL&&W_b == NULL) {
        delete[] W_c;
        delete[] W_b;
    }

    double tol= 1.0e-5;
//    if(A!= nullptr)
//        tol =  (std::max(H->oneNorm(),A->oneNorm())+1)*1.0e-6;
//    else
//        tol =  (H->oneNorm()+1)*1.0e-6;
//
    if(qpOptimalStatus_.KKT_error>tol) {
        //    printf("comp_violation %10e\n", compl_violation);
        //    printf("stat_violation %10e\n", statioanrity_violation);
        //    printf("prim_violation %10e\n", primal_violation);
        //    printf("dual_violation %10e\n", dual_violation);
        //    printf("KKT_error %10e\n", qpOptimalStatus_.KKT_error);
        return false;

    }

    return true;


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
    case QPSOLVER_OPTIMAL:
        return QP_OPTIMAL;
    default:
        return QPERROR_UNKNOWN;
    }
}

void QOREInterface::get_working_set(ActiveType* W_constr, ActiveType* W_bounds) {
    QPGetIntVector(solver_,"workingset",working_set_);
    for(int i=0; i<nConstr_QP_+nVar_QP_; i++) {
        if(i<nVar_QP_) {
            switch (working_set_[i]) {
            case -1:
                if(fabs(x_qp_->values(i)-lb_->values(i))<sqrt_m_eps)
                    W_bounds[i] = ACTIVE_BOTH_SIDE;
                else
                    W_bounds[i] = ACTIVE_ABOVE;
                break;
            case 1:
                if(fabs(x_qp_->values(i)-ub_->values(i))<sqrt_m_eps)
                    W_bounds[i] = ACTIVE_BOTH_SIDE;
                else
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
            case -1:
                if(fabs(x_qp_->values(i)-lb_->values(i)<sqrt_m_eps))
                    W_constr[i-nVar_QP_] = ACTIVE_BOTH_SIDE;
                else
                    W_constr[i-nVar_QP_] = ACTIVE_ABOVE;
                break;
            case 1:
                if(fabs(x_qp_->values(i)-ub_->values(i)<sqrt_m_eps))
                    W_constr[i-nVar_QP_] =ACTIVE_BOTH_SIDE;
                else
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



void QOREInterface::WriteQPDataToFile(Ipopt::EJournalLevel level,
                                      Ipopt::EJournalCategory category,
                                      const string filename) {


#if DEBUG
#if PRINT_OUT_QP_WITH_ERROR
#if PRINT_QP_IN_CPP
    jnlst_->DeleteAllJournals();
    Ipopt::SmartPtr<Ipopt::Journal> QPdata_jrnl= jnlst_->AddFileJournal("QPdata",
            "QORE_cpp_"+filename,Ipopt::J_WARNING);
    QPdata_jrnl->SetAllPrintLevels(level);
    QPdata_jrnl->SetPrintLevel(category,level);
    jnlst_->Printf(level, category, "#include <stdio.h>\n"
                   "#include <assert.h>\n"
                   "#include <stdlib.h>\n"
                   "#include <string.h>\n"
                   "#include <math.h>\n"
                   "#include <matrixconversion.h>\n"
                   "#include <qpsolver.h>\n"
                   "\n"
                   "int main(){\n"
                   "#define NV %i\n#define NC %i\n", nVar_QP_, nConstr_QP_);
    lb_->write_to_file("lb",jnlst_,level,category,QORE);
    ub_->write_to_file("ub",jnlst_,level,category,QORE);
    g_->write_to_file("g",jnlst_,level,category,QORE);
    A_->write_to_file("A",jnlst_,level,category,QORE);
    H_->write_to_file("H",jnlst_,level,category,QORE);

    jnlst_->Printf(level,category,"QoreProblem * qp = 0;\n");
    jnlst_->Printf(level,category,"qp_int rv = QPNew( &qp, NV, NC, %i, %i );\n",
                   A_->EntryNum(), H_->EntryNum());
    jnlst_->Printf(level,category,"assert( rv == QPSOLVER_OK );\n");
    jnlst_->Printf(level,category,"assert( qp!= 0 );\n");
    jnlst_->Printf(level,category,"QPSetInt(qp, \"prtfreq\", 0); \n");
    // pass problem data to solver
    jnlst_->Printf(level, category, "rv = QPSetData( qp, NV, NC, A_jc, A_ir, A_val, H_jc,"
                   " H_ir, H_val );\n");

    jnlst_->Printf(level,category,"assert( rv == QPSOLVER_OK );\n");
    // solve first QP
    jnlst_->Printf(level, category,"rv = QPOptimize( qp, lb, ub, g, 0, 0 );\n");

    jnlst_->Printf(level,category,"qp_int status;\n");
    jnlst_->Printf(level,category,"QPGetInt(qp, \"status\", &status);\n");
    jnlst_->Printf(level,category,"if(status!=QPSOLVER_OPTIMAL){\n)");
    jnlst_->Printf(level,category,"    printf(\"Warning! The QP is not solved to optimality!\");\n");
    jnlst_->Printf(level,category,"}\n");
    jnlst_->Printf(level,category,"assert( rv == QPSOLVER_OK );\n");
    // get and print primal solution

    jnlst_->Printf(level, category, "QPFree(&qp);\n"
                   "\n return 0; \n"
                   "}\n");

#else
#if PRINT_DATA_FOR_QPOASES
    jnlst_->DeleteAllJournals();
    Ipopt::SmartPtr<Ipopt::Journal> QPdata_qpOASES_jrnl= jnlst_->AddFileJournal("QPdata",
            "QPOASES_"+filename,Ipopt::J_WARNING);
    QPdata_qpOASES_jrnl->SetAllPrintLevels(level);
    QPdata_qpOASES_jrnl->SetPrintLevel(category,level);

    jnlst_->Printf(level, category, "%d\n", nVar_QP_);
    jnlst_->Printf(level, category, "%d\n", nConstr_QP_);
    jnlst_->Printf(level, category, "%d\n", A_->EntryNum());
    jnlst_->Printf(level, category, "%d\n", H_->EntryNum());
    lb_->write_to_file("lb",jnlst_,level,category,QORE);
    ub_->write_to_file("ub",jnlst_,level,category,QORE);
    g_->write_to_file("g",jnlst_,level,category,QORE);

    shared_ptr<Vector> A_dense = make_shared<Vector>(nVar_QP_*nConstr_QP_);
    shared_ptr<Vector> H_dense = make_shared<Vector>(nVar_QP_*nVar_QP_);

    A_->get_dense_matrix(A_dense->values());
    shared_ptr<SpHbMat> A_out = make_shared<SpHbMat>(A_dense->values(),nConstr_QP_,nVar_QP_,true,false);
    A_out->write_to_file("A",jnlst_,level,category,QPOASES);

    H_->get_dense_matrix(H_dense->values());

    shared_ptr<SpHbMat> H_out = make_shared<SpHbMat>(H_dense->values(),nVar_QP_,nVar_QP_,true,false);
    H_out->write_to_file("H",jnlst_,level,category,QPOASES);
#endif

#if PRINT_DATA_FOR_QORE
    jnlst_->DeleteAllJournals();
    Ipopt::SmartPtr<Ipopt::Journal> QPdata_qore_jrnl= jnlst_->AddFileJournal("QPdata",
            "QORE_"+filename,Ipopt::J_WARNING);
    QPdata_qore_jrnl->SetAllPrintLevels(level);
    QPdata_qore_jrnl->SetPrintLevel(category,level);

    jnlst_->Printf(level, category, "%d\n", nVar_QP_);
    jnlst_->Printf(level, category, "%d\n", nConstr_QP_);
    jnlst_->Printf(level, category, "%d\n", A_->EntryNum());
    jnlst_->Printf(level, category, "%d\n", H_->EntryNum());
    lb_->write_to_file("lb",jnlst_,level,category,QORE);
    ub_->write_to_file("ub",jnlst_,level,category,QORE);
    g_->write_to_file("g",jnlst_,level,category,QORE);
    A_->write_to_file("A",jnlst_,level,category,QORE);
    H_->write_to_file("H",jnlst_,level,category,QORE);

#endif
#endif
#endif
    jnlst_->DeleteAllJournals();
#endif

}

void QOREInterface::handle_error(QPType qptype, shared_ptr<Stats> stats) {
    switch(status_) {
    case QPSOLVER_OPTIMAL:
        //do nothing here
        break;
    case QPSOLVER_INFEASIBLE:
        shared_ptr<Vector> x_0 = make_shared<Vector>(nVar_QP_);
        //setup the slack variables to satisfy the bound constraints
        for(int i=0; i<nConstr_QP_; i++) {
            x_0->setValueAt(i+nVar_QP_-2*nConstr_QP_,max(0.0,lb_->values(nVar_QP_+i)));
            x_0->setValueAt(i+nVar_QP_-nConstr_QP_,-min(0.0,ub_->values(nVar_QP_+i)));
        }

        rv_ = QPOptimize(solver_, lb_->values(), ub_->values(), g_->values(),
                         x_0->values(), NULL);//
        QPGetInt(solver_, "itercount", qpiter_);
        if(stats!=nullptr)
            stats->qp_iter_addValue(qpiter_[0]);

        assert(rv_ == QPSOLVER_OK);
        break;
    }
}


void QOREInterface::set_solver_options(shared_ptr<const Options> options) {

    if(options->qpPrintLevel==0) {
        //does not print anything
        QPSetInt(solver_, "prtfreq", -1);
    }
    QPSetInt(solver_,"maxiter",options->qp_maxiter);

}
}//SQP_HOTSTART

