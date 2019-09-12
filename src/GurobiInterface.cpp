/* Copyright (C) 2019
 * All Rights Reserved.
 *
 * Authors: Xinyi Luo
 * Date:2019-09-06
 */
#include <sqphot/GurobiInterface.hpp>
#ifdef GUROBI

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
    qobj_(0),
    lterm_(0)
{
    grb_env_ = new GRBEnv();
    grb_mod_ = new GRBModel(*grb_env_);
    grb_constr.resize(2*nConstr_QP_);
    set_solver_options();
    grb_vars_ = grb_mod_->addVars(nVar_QP_);
    x_qp = make_shared<Vector>(nVar_QP_);
    y_qp = make_shared<Vector>(nConstr_QP_);
}

GurobiInterface::~GurobiInterface() {
    grb_constr.clear();
    delete grb_mod_;
    grb_mod_ = nullptr;
    delete grb_env_;
    grb_env_ = nullptr;
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
    grb_mod_->setObjective(qobj_+lterm_,GRB_MINIMIZE);
    grb_mod_->update();
    try {
        grb_mod_->optimize();
    }
    catch(GRBException exception) {
        THROW_EXCEPTION(GRB_SOLVER_FAILS,"Gurobi Fails due to internal errors");
    }
    for(int i =0; i<nVar_QP_; i++) {
        x_qp->setValueAt(i,grb_vars_[i].get(GRB_DoubleAttr_X));
    }

    if(grb_mod_->get(GRB_IntAttr_Status)!=GRB_OPTIMAL) {
        printf("qp is not optimal, the current status is %i", grb_mod_->get
               (GRB_IntAttr_Status));
    }
    for(int i=0; i<nConstr_QP_*2; i++) {
        if(grb_mod_->getConstr(i).get(GRB_DoubleAttr_Pi)>0)
            y_qp->setValueAt((int)i/2,grb_mod_->getConstr(i).get(GRB_DoubleAttr_Pi));
        else if(grb_mod_->getConstr(i).get(GRB_DoubleAttr_Pi)<0)
            y_qp->setValueAt((int)i/2,grb_mod_->getConstr(i).get(GRB_DoubleAttr_Pi));
    }



}



/**
 * @brief Solve a regular LP with given data and options
 *
 */

void GurobiInterface::optimizeLP(shared_ptr<Stats> stats) {
    grb_mod_->setObjective(lterm_,GRB_MINIMIZE);
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
double* GurobiInterface::get_multipliers_bounds()  {
    return nullptr;
}



/**
 * @brief get the pointer to the multipliers to the regular constraints.
 */
double* GurobiInterface::get_multipliers_constr()  {
    return y_qp->values();
}

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
//    A_->print_full("A");
    linlhs =0;
    for(int i=0; i<A_->EntryNum(); i++) {
        if(A_->RowIndex(i)==location+1) {
            lterm = A_->MatVal(i)*grb_vars_[A_->ColIndex(i)-1];
        }
        linlhs += lterm;
    }
    lterm = grb_vars_[A_->ColNum()+location];

    lterm += -1.0*grb_vars_[A_->ColNum()+nConstr_QP_+location];
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
    lterm_.remove(grb_vars_[location]);
    lterm_+=value*grb_vars_[location];
}
//@}


/**@name Setters for matrix*/
//@{
void GurobiInterface::set_H_structure(shared_ptr<const SpTripletMat> rhs) {}

void GurobiInterface::set_H_values(shared_ptr<const SpTripletMat> rhs) {
    qobj_=0;

    for(int i =0; i<rhs->EntryNum(); i++) {
        if(rhs->isSymmetric()) {
            if (rhs->ColIndex(i) == rhs->RowIndex(i)) {
                qobj_ += 0.5 * rhs->MatVal(i) * grb_vars_[rhs->ColIndex(i)
                         - 1] *
                         grb_vars_[rhs->RowIndex(i) - 1];
            }
            else {
                qobj_ += rhs->MatVal(i) * grb_vars_[rhs->ColIndex(i) - 1] *
                         grb_vars_[rhs->RowIndex(i) - 1];
            }

        }
    }
}

void GurobiInterface::set_A_structure(shared_ptr<const SpTripletMat> rhs, Identity2Info
                                      I_info) {
}

void GurobiInterface::reset_constraints() {
    for(int i =0; i<nConstr_QP_; i++) {
        grb_mod_->remove(grb_mod_->getConstrByName("lbA_"+to_string(i)));
        grb_mod_->remove(grb_mod_->getConstrByName("ubA_"+to_string(i)));
    }

}

void GurobiInterface::set_A_values(shared_ptr<const SpTripletMat> rhs, Identity2Info
                                   I_info) {

    A_ = rhs;
    I_info_ = I_info;
}
//@}

/**@name Gurobi model setup */
//@{

void GurobiInterface::set_solver_options() {
    grb_mod_->set(GRB_DoubleParam_TimeLimit,1000.0);
    grb_mod_->set(GRB_IntParam_OutputFlag,0);

}

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

/**-------------------------------------------------------**/
/**                  Data Writer                          **/
/**-------------------------------------------------------**/

void GurobiInterface::WriteQPDataToFile(Ipopt::SmartPtr<Ipopt::Journalist> jnlst,
                                        Ipopt::EJournalLevel level,
                                        Ipopt::EJournalCategory category) {}

};
#endif



