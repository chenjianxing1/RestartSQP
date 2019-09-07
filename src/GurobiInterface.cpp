/* Copyright (C) 2019
 * All Rights Reserved.
 *
 * Authors: Xinyi Luo
 * Date:2019-09-06
 */
#include <sqphot/GurobiInterface.hpp>

namespace SQPhotstart {

GurobiInterface::GurobiInterface(Index_info nlp_info,
                                 QPType qptype,
                                 shared_ptr<const Options> options,
                                 Ipopt::SmartPtr<Ipopt::Journalist> jnlst):
    status_(UNSOLVED),
    qptype_(qptype),
    nConstr_QP_(nlp_info.nCon),
    nVar_QP_(nlp_info.nVar+2*nlp_info.nCon),
    options_(options),
    grb_mod_(NULL),
    A_(NULL),
    qobj_(0)
{
    grb_env_ = new GRBEnv();
    grb_mod_ = new GRBModel(*grb_env_);

    set_solver_options();
    grb_vars_ = grb_mod_->addVars(nVar_QP_);
    x_qp = make_shared<Vector>(nVar_QP_);
    y_qp = make_shared<Vector>(nConstr_QP_);
}

GurobiInterface::~GurobiInterface() {
    if (grb_mod_) delete grb_mod_;
    delete grb_env_;
}

void GurobiInterface::reset_model() {
    if (grb_mod_ != NULL) delete grb_mod_;
//    grb_vars_.clear();
    grb_mod_ = new GRBModel(*grb_env_);
    qobj_ = 0;
}

/**
 * @brief Solve a regular QP with given data and options.
 */
void GurobiInterface::optimizeQP(shared_ptr<Stats> stats)  {
    grb_mod_->setObjective(qobj_,GRB_MINIMIZE);
    grb_mod_->update();
    grb_mod_->optimize();
    for(int i =0; i<nVar_QP_; i++) {
        x_qp->setValueAt(i,grb_vars_[i].get(GRB_DoubleAttr_X));
    }

    x_qp->print("x_qp");
    for(int i=0; i<nConstr_QP_; i++) {
        y_qp->setValueAt(i,grb_mod_->getConstr(i).get(GRB_DoubleAttr_Pi));
    }

    y_qp->print("y_qp");

}



/**
 * @brief Solve a regular LP with given data and options
 *
 */

void GurobiInterface::optimizeLP(shared_ptr<Stats> stats) {
    grb_mod_->setObjective(qobj_,GRB_MINIMIZE);
    grb_mod_->optimize();
}


/**-------------------------------------------------------**/
/**                    Getters                            **/
/**-------------------------------------------------------**/
/**@name Getters*/
//@{
/**
 * @return the pointer to the optimal solution
 *
 */
double* GurobiInterface::get_optimal_solution()  {
    return x_qp->values();

}

/**
 *@brief get the objective value from the QP solvers
 *
 * @return the objective function value of the QP problem
 */
double GurobiInterface::get_obj_value()  {
    return grb_mod_->get(GRB_DoubleAttr_ObjVal);
}


/**
 * @brief get the pointer to the multipliers to the bounds constraints.
 */
double* GurobiInterface::get_multipliers_bounds()  {}

/**
 * @brief get the pointer to the multipliers to the regular constraints.
 */
double* GurobiInterface::get_multipliers_constr()  {}

/**
 * @brief copy the working set information
 * @param W_constr a pointer to an array of length (nCon_QP_) which will store the
 * working set for constraints
 * @param W_bounds a pointer to an array of length (nVar_QP_) which will store the
 * working set for bounds
 *
 */
void GurobiInterface::get_working_set(ActiveType* W_constr, ActiveType* W_bounds) {}

QPReturnType GurobiInterface::get_status() {}

//@}

/**-------------------------------------------------------**/
/**                    Setters                            **/
/**-------------------------------------------------------**/
/**@name Setters, by location and value*/
//@{
void GurobiInterface::set_lb(int location, double value) {

    if(value>-INF) {
        grb_vars_[location].set(GRB_DoubleAttr_LB, value);
    }

}

void GurobiInterface::set_ub(int location, double value) {
    if(value<INF) {
        grb_vars_[location].set(GRB_DoubleAttr_UB, value);
    }
}

void GurobiInterface::set_lbA(int location, double value) {
    assert(A_!=NULL);
    GRBLinExpr lterm, linlhs; //linear term
    A_->print_full("A");
    linlhs =0;
    for(int i=0; i<A_->EntryNum(); i++) {
        if(A_->RowIndex(i)==location+1) {

            lterm = A_->MatVal(i)*grb_vars_[A_->ColIndex(i)-1];
//                 printf("On Row % i : lterm = A_->MatVal(%i)*grb_vars_[%i]\n",
//                         location+1, i,
//                         A_->ColIndex(i)
//                 -1);
            //TODOL: check the index...
        }
        linlhs += lterm;
    }
    lterm = grb_vars_[A_->ColNum()+location];
//         printf("lterm = grb_vars_[%i]\n",A_->ColNum()+location);

    lterm += -1.0*grb_vars_[A_->ColNum()+nConstr_QP_+location];
//         printf("lterm = grb_vars_[%i]\n",A_->ColNum()+location+nConstr_QP_);
    linlhs += lterm;
    linlhs-=value;

    grb_mod_->addConstr(linlhs,GRB_GREATER_EQUAL,0,"lbA_"+to_string(location));

}

void GurobiInterface::set_ubA(int location, double value) {

    assert(A_!=NULL);
    GRBLinExpr lterm, linlhs; //linear term

    linlhs =0;
    for(int i=0; i<A_->EntryNum(); i++) {
        if(A_->RowIndex(i)==location+1) {
            lterm = A_->MatVal(i)*grb_vars_[A_->ColIndex(i)-1];
            //TODOL: check the index...
        }
        linlhs += lterm;
    }
    lterm = grb_vars_[A_->ColNum()+location];
    lterm += -1.0*grb_vars_[A_->ColNum()+nConstr_QP_+location];
    linlhs += lterm;

    linlhs-=value;
    grb_mod_->addConstr(linlhs,GRB_LESS_EQUAL,0,"ubA_"+to_string(location));

}

void GurobiInterface::set_g(int location, double value) {
    qobj_ +=value*grb_vars_[location];
}
//@}


/**@name Setters for matrix*/
//@{
void GurobiInterface::set_H_structure(shared_ptr<const SpTripletMat> rhs) {}

void GurobiInterface::set_H_values(shared_ptr<const SpTripletMat> rhs) {
    rhs->print_full("H");

    for(int i =0; i<rhs->EntryNum(); i++) {
        if(rhs->isSymmetric()) {
            if (rhs->ColIndex(i) == rhs->RowIndex(i)) {
                qobj_ += 0.5 * rhs->MatVal(i) * grb_vars_[rhs->ColIndex(i)
                         - 1] *
                         grb_vars_[rhs->RowIndex(i) - 1];
//                    printf("obj+=0.5*%10e * x[%i]*x[%i]\n",rhs->MatVal(i),rhs->ColIndex(i)
//                    -1,rhs->RowIndex(i)-1);
            }
            else {
                qobj_ += rhs->MatVal(i) * grb_vars_[rhs->ColIndex(i) - 1] *
                         grb_vars_[rhs->RowIndex(i) - 1];
//                    printf("obj+=%10e * x[%i]*x[%i]\n", rhs->MatVal(i), rhs->ColIndex(i)
//                                                                       - 1,
//                           rhs->RowIndex(i) - 1);
            }

        }
    }
}

void GurobiInterface::set_A_structure(shared_ptr<const SpTripletMat> rhs, Identity2Info
                                      I_info) {
}

void GurobiInterface::set_A_values(shared_ptr<const SpTripletMat> rhs, Identity2Info
                                   I_info) {


    A_ = rhs;
    I_info_ = I_info;
}
//@}

/**@name Gurobi model setup */
//@{

void GurobiInterface::set_solver_options() {}

//@}




/**@name Setters for dense vector, by vector value*/
//@{
void GurobiInterface::set_ub(shared_ptr<const Vector> rhs) {
}


void GurobiInterface::set_lb(shared_ptr<const Vector> rhs) {

}

void GurobiInterface::set_lbA(shared_ptr<const Vector> rhs) {
}

void GurobiInterface::set_ubA(shared_ptr<const Vector> rhs) {
}

void GurobiInterface::set_g(shared_ptr<const Vector> rhs) {}
//@}

#if DEBUG
/**@name Getters for private members*/
//@{
const shared_ptr<Vector>& GurobiInterface::getLb() const {}

const shared_ptr<Vector>& GurobiInterface::getUb() const {}

const shared_ptr<Vector>& GurobiInterface::getLbA() const {}

const shared_ptr<Vector>& GurobiInterface::getUbA() const  {}

const shared_ptr<Vector>& GurobiInterface::getG() const {}

const shared_ptr<const SpTripletMat> GurobiInterface::getH() const  {}

const shared_ptr<const SpTripletMat> GurobiInterface::getA() const  {}
//@}

#endif
/**-------------------------------------------------------**/
/**                  Data Writer                          **/
/**-------------------------------------------------------**/

void GurobiInterface::WriteQPDataToFile(Ipopt::SmartPtr<Ipopt::Journalist> jnlst,
                                        Ipopt::EJournalLevel level,
                                        Ipopt::EJournalCategory category) {}

};



