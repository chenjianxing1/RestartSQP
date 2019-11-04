/* Copyright (C) 2019
 * All Rights Reserved.
 *
 * Authors: Xinyi Luo
 * Data: 2019-09-09
 */

#include <sqphot/CplexInterface.hpp>

namespace SQPhotstart {
/** Constructor*/
CplexInterface::CplexInterface(NLPInfo nlp_info,
                               QPType qptype,
                               shared_ptr<const Options> options,
                               Ipopt::SmartPtr<Ipopt::Journalist> jnlst):
    status_(UNSOLVED),
    qptype_(qptype),
    nConstr_QP_(nlp_info.nCon),
    nVar_QP_(nlp_info.nVar+2*nlp_info.nCon),
    options_(options),
    A_(NULL),
    firstQPsolved_(false)
{
    cplex_env_ = make_shared<IloEnv>();
    cplex_model_ = make_shared<IloModel>(*cplex_env_);

    x_qp_ = make_shared<Vector>(nVar_QP_);
    y_qp_ = make_shared<Vector>(nConstr_QP_);

    cplex_vars_.resize(nVar_QP_);
    qobj_ = IloNumExpr(*cplex_env_);

    constraints_.resize(nConstr_QP_*2);
    lterm_.resize(nVar_QP_);
    for(int i =0; i<nVar_QP_; i++) {
        cplex_vars_.at(i) = IloNumVar(*cplex_env_);
        lterm_.at(i) = IloNumExpr(*cplex_env_);
    }
    for(int i = 0; i<nConstr_QP_*2;  i++)
        constraints_.at(i) = IloNumExpr(*cplex_env_);

    lb_ = make_shared<Vector>(nVar_QP_);
    ub_ = make_shared<Vector>(nVar_QP_);

    set_solver_options();
}


/** Destructor */
CplexInterface::~CplexInterface() {
    lterm_.clear();
//        delete qobj_; qobj_ = nullptr;
//        delete lterm_; lterm_ = nullptr;
}



void CplexInterface::optimizeQP(shared_ptr<Stats> stats)  {
    //setup the lower and upper bounds for the varaibles

    IloNumArray lb(*cplex_env_,nVar_QP_);
    IloNumArray ub(*cplex_env_,nVar_QP_);
    IloNumExpr obj(*cplex_env_);
    IloNumVarArray vars(*cplex_env_);
    IloNumArray x_start(*cplex_env_);
    IloNumArray y_start(*cplex_env_);
    IloCplex cplex(*cplex_env_);

    IloRangeArray constraints(*cplex_env_);



    for(int i = 0; i<nConstr_QP_; i++) {
        constraints.add(IloRange(*cplex_env_,0,constraints_[i]));
        constraints.add(IloRange(*cplex_env_,constraints_[i+nConstr_QP_],0));
    }
    cplex_model_->add(constraints);


    for(int i =0; i<nVar_QP_; i++) {
        lb[i] = lb_->values(i);
        ub[i] = ub_->values(i);
        cplex_vars_[i].setBounds(lb[i],ub[i]);
        vars.add(cplex_vars_[i]);

//        if(firstQPsolved_){
//            x_start.add(x_qp_->values(i));
//        }
//        else
        x_start.add(0);

    }
//    for(int i = 0; i<nConstr_QP_*2; i++){
//     y_start.add(0);
//    }
    //dual startint point? TODO
//    for(int i = 0; i<nConstr_QP_*2; i++) {
//        y_start.add(y_qp_->values(i));
//    }
    //setup the objectives

    cplex.setOut(cplex_env_->getNullStream());
    cplex.setParam(IloCplex::OptimalityTarget,2);
    obj = qobj_;
    for(int i = 0; i< nVar_QP_; i++) {
        obj += lterm_[i];
    }

    auto cplex_obj_ = IloMinimize(*cplex_env_,obj);
    cplex_model_->add(cplex_obj_);
    cplex.extract(*cplex_model_);

    cplex.setStart(x_start, 0, vars, 0, 0, 0);
    try {
        cplex.solve();
    }
    catch (IloException& e) {
        cerr << "Concert exception caught: " << e << endl;
    }
    if(cplex.getStatus()!=IloAlgorithm::Optimal) {
//        THROW_EXCEPTION(CPLEX_SOLVER_FAILS,"CPLEX does not solve the QP to optimal");
    }
    //get the primal solution

    for(int i = 0; i< nVar_QP_; i++) {
        if(cplex.isExtracted(vars[i]))
            x_qp_->setValueAt(i, cplex.getValue(vars[i]));
    }


    //get the dual solution
    try {
        cplex.getDuals(y_start,constraints);
    }
    catch (IloException& e) {
        cerr << "Concert exception caught: " << e << endl;
    }
    for(int i= 0; i<nConstr_QP_; i++) {
        if(y_start[i]>0)
            y_qp_->setValueAt(i,y_start[i]);
        else
            y_qp_->setValueAt(i,y_start[i*2+1]);
    }

//    y_qp_->print("y_qp");
    final_obj_ = cplex.getValue(obj);

    cplex_model_->remove(cplex_obj_);

}


void CplexInterface::optimizeLP(shared_ptr<Stats> stats)  {
    //setup the lower and upper bounds for the varaibles

    IloNumArray lb(*cplex_env_,nVar_QP_);
    IloNumArray ub(*cplex_env_,nVar_QP_);
    IloNumExpr obj(*cplex_env_);//the true objective of the QP
    IloCplex cplex(*cplex_env_);

    for(int i =0; i<nVar_QP_; i++) {
        lb[i] = lb_->values(i);
        ub[i] = ub_->values(i);
        cplex_vars_[i].setBounds(lb[i],ub[i]);
    }
    obj = 0;
    for(int i = 0; i< nVar_QP_; i++) {
        obj += lterm_.at(i);
    }
    auto cplex_obj_ = IloMinimize(*cplex_env_,obj);
    cplex_model_->add(cplex_obj_);
    cplex.extract(*cplex_model_);
    cplex.setOut(cplex_env_->getNullStream());
    cplex.solve();
    if(cplex.getStatus()!=IloAlgorithm::Optimal) {
        THROW_EXCEPTION(CPLEX_SOLVER_FAILS, "CPLEX does not solve the QP to optimal");
    }

    for(int i = 0; i< nVar_QP_; i++) {
        if(cplex.isExtracted(cplex_vars_[i]))
            x_qp_->setValueAt(i, cplex.getValue(cplex_vars_[i]));
    }
    x_qp_->print("x_qp");

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
double* CplexInterface::get_optimal_solution()  {
    return x_qp_->values();

}

/**
 *@brief get the objective value from the QP solvers
 *
 * @return the objective function value of the QP problem
 */
double CplexInterface::get_obj_value()  {
    return final_obj_;
}


/**
 * @brief get the pointer to the multipliers to the bounds constraints.
 */
double* CplexInterface::get_multipliers_bounds()  {
    return nullptr;
}



/**
 * @brief get the pointer to the multipliers to the regular constraints.
 */
double* CplexInterface::get_multipliers_constr()  {
    return y_qp_->values();
}

/**
 * @brief copy the working set information
 * @param W_constr a pointer to an array of length (nCon_QP_) which will store the
 * working set for constraints
 * @param W_bounds a pointer to an array of length (nVar_QP_) which will store the
 * working set for bounds
 *
 */
void CplexInterface::get_working_set(ActiveType* W_constr, ActiveType* W_bounds) {}

QPReturnType CplexInterface::get_status() {}

//@}

/**-------------------------------------------------------**/
/**                    Setters                            **/
/**-------------------------------------------------------**/
/**@name Setters, by location and value*/
//@{
void CplexInterface::set_lb(int location, double value) {
    lb_->setValueAt(location, value);
}

void CplexInterface::set_ub(int location, double value) {
    ub_->setValueAt(location, value);
}

void CplexInterface::set_lbA(int location, double value) {
    assert(A_!=NULL);

    constraints_.at(location) = cplex_vars_[A_->ColNum()+location];
    constraints_.at(location) += -1.0*cplex_vars_[A_->ColNum()+nConstr_QP_+location];
    constraints_.at(location) -= value;
    for(int i=0; i<A_->EntryNum(); i++) {
        if(A_->RowIndex(i)==location+1) {
            constraints_.at(location) += A_->MatVal(i)*cplex_vars_[A_->ColIndex(i)-1];
        }
    }
//    cplex_model_->add(constraints_.at(location)>=0);

}

void CplexInterface::set_ubA(int location, double value) {
    assert(A_ != NULL);
    int loc =  location +nConstr_QP_;
    constraints_.at(loc) = cplex_vars_[A_->ColNum()+location];
    constraints_.at(loc) += -1.0*cplex_vars_[A_->ColNum()+nConstr_QP_+location];
    constraints_.at(loc) -= value;
    for(int i=0; i<A_->EntryNum(); i++) {
        if(A_->RowIndex(i)==location+1) {
            constraints_.at(loc) += A_->MatVal(i)*cplex_vars_[A_->ColIndex(i)-1];
        }
    }
//    cplex_model_->add(constraints_.at(loc) <= 0);
}

void CplexInterface::set_g(int location, double value) {
    lterm_.at(location) = value* cplex_vars_[location];

}

//@}


/**@name Setters for matrix*/
//@{
void CplexInterface::set_H_structure(shared_ptr<const SpTripletMat> rhs) {}

void CplexInterface::set_H_values(shared_ptr<const SpTripletMat> rhs) {

    for(int i =0; i<rhs->EntryNum(); i++) {
        if(rhs->isSymmetric()) {
            if (rhs->ColIndex(i) == rhs->RowIndex(i)) {
                qobj_ += 0.5 * rhs->MatVal(i) * cplex_vars_[rhs->ColIndex(i)
                         - 1] *
                         cplex_vars_[rhs->RowIndex(i) - 1];
            }
            else {
                qobj_ += rhs->MatVal(i) * cplex_vars_[rhs->ColIndex(i) - 1] *
                         cplex_vars_[rhs->RowIndex(i) - 1];
            }

        }
    }
}

void CplexInterface::set_A_structure(shared_ptr<const SpTripletMat> rhs, IdentityInfo
                                     I_info) {
}

void CplexInterface::reset_constraints() {

    for(int i = 0; i<nConstr_QP_; i++) {
        cplex_model_->remove(constraints_[i] >= 0);
        cplex_model_->remove(constraints_[i+nConstr_QP_] <= 0);
    }

}

void CplexInterface::set_A_values(shared_ptr<const SpTripletMat> rhs, IdentityInfo
                                  I_info) {
    A_ = rhs;
    I_info_ = I_info;
}
//@}

/**@name Cplex model setup */
//@{

void CplexInterface::set_solver_options() {

}

//@}




/**@name Setters for dense vector, by vector value*/
//@{
void CplexInterface::set_ub(shared_ptr<const Vector> rhs) {}


void CplexInterface::set_lb(shared_ptr<const Vector> rhs) {}

void CplexInterface::set_lbA(shared_ptr<const Vector> rhs) {}

void CplexInterface::set_ubA(shared_ptr<const Vector> rhs) {}

void CplexInterface::set_g(shared_ptr<const Vector> rhs) {}


//@}

/**@name Getters for private members*/
//@{
const shared_ptr<Vector>& CplexInterface::getLb() const {}

const shared_ptr<Vector>& CplexInterface::getUb() const {}

const shared_ptr<Vector>& CplexInterface::getLbA() const {}

const shared_ptr<Vector>& CplexInterface::getUbA() const {}

const shared_ptr<Vector>& CplexInterface::getG() const {}

const shared_ptr<const SpTripletMat> CplexInterface::getH() const {}

const shared_ptr<const SpTripletMat> CplexInterface::getA() const {}
//@}

}

