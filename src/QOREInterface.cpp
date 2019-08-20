/* Copyright (C) 2019
 *
 * Authors: Xinyi Luo
 * Date:    2019-07
 * Date:    2019-08-15
 */

#include <sqphot/QOREInterface.hpp>

namespace SQPhotstart {

QOREInterface::QOREInterface(Index_info nlp_info, QPType qptype) :
    firstQPsolved_(false),
    qp_(0),
    nVar_QP_(nlp_info.nVar + nlp_info.nCon * 2),
    nConstr_QP_(nlp_info.nCon + nlp_info.nVar) {
    rv_ = allocate(nlp_info, qptype);
    assert(rv_ == QPSOLVER_OK);

}


QOREInterface::~QOREInterface() {
    QPFree(&qp_);
}

void QOREInterface::optimizeQP(shared_ptr<Stats> stats, shared_ptr<Options>) {
    rv_ = QPSetData(qp_, nVar_QP_, nConstr_QP_, A_->ColIndex(), A_->RowIndex(),
                    A_->MatVal(), H_->ColIndex(), H_->RowIndex(), H_->MatVal());
    A_->print();
    H_->print();
    assert(rv_ == QPSOLVER_OK);
    QPSetInt(qp_, "prtfreq", 0);
    if (!firstQPsolved_) {
        rv_ = QPOptimize(qp_, lb_->values(), ub_->values(), g_->values(),0,0);//
        firstQPsolved_ = true;
    } else
        rv_ = QPOptimize(qp_, lb_->values(), ub_->values(), g_->values(),
                         0,0);//TODO: does not update primal-dual sol here
    assert(rv_ == QPSOLVER_OK);

}


void QOREInterface::optimizeLP(shared_ptr<Stats> stats, shared_ptr<Options>) {
}


int QOREInterface::allocate(Index_info nlp_info, QPType qptype) {

    int nnz_g_QP = nlp_info.nnz_jac_g +
                   2 * nlp_info.nCon;//number of nonzero variables in jacobian
    //The Jacobian has the structure [J I -I], so it will contains extra 2*number_constr nonzero elements
    lb_ = make_shared<Vector>(nVar_QP_ + nConstr_QP_);
    ub_ = make_shared<Vector>(nVar_QP_ + nConstr_QP_);
    g_ = make_shared<Vector>(nVar_QP_);
    A_ = make_shared<SpHbMat>(nnz_g_QP, nConstr_QP_, nVar_QP_);
    H_ = make_shared<SpHbMat>(nnz_g_QP, nConstr_QP_, nVar_QP_);
    x_qp_ = make_shared<Vector>(nVar_QP_);
    y_qp_ = make_shared<Vector>(nConstr_QP_+nVar_QP_);
    int returnvalue;
    if (qptype != LP) {
        returnvalue = QPNew(&qp_, nVar_QP_, nConstr_QP_, nnz_g_QP,
                            nlp_info.nnz_h_lag);
        H_ = make_shared<SpHbMat>(nVar_QP_, nVar_QP_, true);
    } else {
        //if we are solving an LP, the number of nonzero in the Hessian is 0

        returnvalue = QPNew(&qp_, nVar_QP_, nConstr_QP_, nnz_g_QP, 0);
    }
    return returnvalue;
}

double* QOREInterface::get_optimal_solution() {
    rv_=QPGetDblVector(qp_,"primalsol",x_qp_->values());
    assert(rv_==QPSOLVER_OK);
    return x_qp_->values();
}

double QOREInterface::get_obj_value() {
//FIXME: QORE does not have existing emthod to return the qp obj, now it is calculated
// in Algorithm class for QOREInterface...
}


double* QOREInterface::get_multipliers() {
    rv_=QPGetDblVector(qp_,"dualsol",y_qp_->values());
    assert(rv_==QPSOLVER_OK);
    return y_qp_->values();
}

void QOREInterface::set_ub(int location, double value) {
    ub_->setValueAt(location, value);

}

void QOREInterface::set_g(int location, double value) {
    g_->setValueAt(location, value);
}

void QOREInterface::set_lb(int location, double value) {
    lb_->setValueAt(location, value);
}

void QOREInterface::set_A_structure(shared_ptr<const SpTripletMat> rhs,
                                    Identity2Info I_info) {

    A_->setStructure(rhs, I_info);
}

void
QOREInterface::set_A_values(shared_ptr<const SpTripletMat> rhs,
                            Identity2Info I_info) {

    A_->setMatVal(rhs, I_info);
}


}//SQP_HOTSTART

