/* Copyright (C) 2019
 *
 * Authors: Xinyi Luo
 * Date:    2019-07
 * Date:    2019-08-15
 */

#include <sqphot/QOREInterface.hpp>

namespace SQPhotstart {

QOREInterface::QOREInterface(Index_info nlp_info, QPType qptype):
    firstQPsolved_(false),
    qp_(0),
    nVar_QP_(nlp_info.nVar+nlp_info.nCon*2),
    nConstr_QP_(nlp_info.nCon+nlp_info.nVar)
{

    rv_ = allocate(nlp_info, qptype);
    assert(rv_==QPSOLVER_OK);

}


QOREInterface::~QOREInterface()
{
    QPFree(&qp_);
}

void QOREInterface::optimizeQP(shared_ptr<Stats> stats, shared_ptr<Options>)
{
    rv_= QPSetData(qp_, nVar_QP_, nConstr_QP_, A_->ColIndex(), A_->RowIndex(),
                   A_->MatVal(), H_->ColIndex(), H_->RowIndex(), H_->MatVal());
    assert(rv_==QPSOLVER_OK);
    QPSetInt(qp_,"prtfreq", 0);
    if(!firstQPsolved_) {
        rv_=QPOptimize(qp_,lb_->values(),ub_->values(),g_->values(),0,0);// FIXME:the last two should be the pointer to x and y(optimal sol and multiplier)
        firstQPsolved_ = true;
    }
    else
        rv_=QPOptimize(qp_,lb_->values(),ub_->values(),g_->values(),0,0);// FIXME:the last two should be the pointer to x and y(optimal sol and multiplier)
    assert(rv_==QPSOLVER_OK);

}


void QOREInterface::optimizeLP(shared_ptr<Stats> stats, shared_ptr<Options>)
{
}



int QOREInterface::allocate(Index_info nlp_info, QPType qptype) {

    int nnz_g_QP = nlp_info.nnz_jac_g+2*nlp_info.nCon;//number of nonzero variables in jacobian
    //The Jacobian has the structure [J I -I], so it will contains extra 2*number_constr nonzero elements
    lb_ = make_shared<Vector>(nVar_QP_+nConstr_QP_);
    ub_ = make_shared<Vector>(nVar_QP_+nConstr_QP_);
    g_ = make_shared<Vector>(nVar_QP_);
    A_= make_shared<SpHbMat>(nnz_g_QP,nConstr_QP_,nVar_QP_);
    H_= make_shared<SpHbMat>(nnz_g_QP,nConstr_QP_,nVar_QP_);
    int returnvalue;
    if(qptype!=LP) {
        returnvalue =QPNew(&qp_, nVar_QP_,nConstr_QP_,nnz_g_QP,
                           nlp_info.nnz_h_lag);
        H_= make_shared<SpHbMat>(nVar_QP_, nVar_QP_, true);
    }
    else {
        returnvalue =QPNew(&qp_, nVar_QP_,nConstr_QP_,nnz_g_QP, 0);
    }
    return returnvalue;
}

void QOREInterface::get_optimal_solution(double* x_optimal) {
    rv_ = QPGetDblVector(qp_,"primalsol", x_optimal);
    assert(rv_==QPSOLVER_OK);
}

double QOREInterface::get_obj_value() {

}


void QOREInterface::get_multipliers(double* y_optimal) {
    rv_ = QPGetDblVector(qp_,"dualsol", y_optimal);
    assert(rv_==QPSOLVER_OK);
}

}//SQP_HOTSTART

