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
        rv_=QPOptimize(qp_,lb_->values(),ub_->values(),g_->values(),x_qp_->values(),
                       y_qp_->values());//
        firstQPsolved_ = true;
    }
    else
        rv_=QPOptimize(qp_,lb_->values(),ub_->values(),g_->values(),x_qp_->values(),
                       y_qp_->values());
    lb_->print("lb_");
    ub_->print("ub_");
    g_->print("g");
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
    x_qp_ = make_shared<Vector>(nVar_QP_);
    y_qp_ = make_shared<Vector>(nConstr_QP_);
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

double* QOREInterface::get_optimal_solution() {

    return x_qp_->values();
}

double QOREInterface::get_obj_value() {
//    std::cout<<(&qp_)<<std::endl;
//    shared_ptr<Vector> Hx = make_shared<Vector>(nVar_QP_);
//    H_->times(x_qp_,Hx);
//    return 0.5*(x_qp_->times(Hx))+g_->times(x_qp_);
}


double* QOREInterface::get_multipliers() {
    return x_qp_->values();
}


}//SQP_HOTSTART

