/* Copyright (C) 2019
 * All Rights Reserved.
 *
 * Authors: Xinyi Luo
 * Date:2019-07
 */
#include <sqphot/QPhandler.hpp>


namespace SQPhotstart {
/**
 * Default constructor
 *
 */
QPhandler::QPhandler() {}

/**
 *Default destructor
 */
QPhandler::~QPhandler() {}

/**
 * Get the optimal solution from the QPhandler_interface
 *
 *This is only an interface for user to avoid call interface directly.
 * @param p_k       the pointer to an empty array with the length equal to the size
 *                  of the QP subproblem
 */
void QPhandler::GetOptimalSolution(double* p_k) {
    qp_interface_->get_optimal_solution(p_k);

}

/**
 *Get the multipliers from the QPhandler_interface
 *
 *This is only an interface for user to avoid call interface directly.
 *
 * @param y_k       the pointer to an empty array with the length equal to the size of
 * multipliers of the QP subproblem
 */
void QPhandler::GetMultipliers(double* y_k) {
    qp_interface_->get_multipliers(y_k);

}


/**
 * This function initializes all objects will be used in this class.
 *
 * @param nlp_info
 *          it contains the information about number of variables, number of
 *          constraints,  number of elements in the Hessian and that of Jacobian. The
 *          definition of Index_info is in Types.hpp
 * @param Constraint_type
 *          it specifies how the variables are bounded. Please check
 *           @ClassifyConstraintType in Algorithm.hpp
 *          for more details
 * @param qptype
 *           it specifies which type of QP is going to be solved. It can be either LP,
 *           or QP, or SOC
 */
void QPhandler::init(Index_info nlp_info, QPType qptype) {
    allocate(nlp_info, qptype);
    nlp_info_ = nlp_info;
    qptype_ = qptype;

}

/**
 * Setup the bounds for the QP subproblems according to the information from current
 * iterate. We have
 * 	c_l -c_k <=J_p+ u-v<=c_u-c_k
 * The bound is formulated as
 *   max(-delta, x_l-x_k)<=p<=min(delta,x_u-x_k)
 * and  u,v>=0
 * @param delta      trust region radius
 * @param x_k        current iterate point
 * @param c_k        current constraint value evaluated at x_k
 * @param x_l        the lower bounds for variables
 * @param x_u        the upper bounds for variables
 * @param c_l        the lower bounds for constraints
 * @param c_u        the upper bounds for constraints
 */
void QPhandler::setup_bounds(double delta, shared_ptr<const Vector> x_k,
                             shared_ptr<const Vector> x_l,
                             shared_ptr<const Vector> x_u,
                             shared_ptr<const Vector> c_k,
                             shared_ptr<const Vector> c_l,
                             shared_ptr<const Vector> c_u) {

    for (int i = 0; i < nlp_info_.nCon; i++) {
        qp_interface_->getLbA()->setValueAt(i, c_l->values()[i] - c_k->values()[i]);
        qp_interface_->getUbA()->setValueAt(i, c_u->values()[i] - c_k->values()[i]);
    }

    for (int i = 0; i < nlp_info_.nVar; i++) {
        qp_interface_->getLb()->setValueAt(i, std::max(
                                               x_l->values()[i] - x_k->values()[i], -delta));
        qp_interface_->getUb()->setValueAt(i, std::min(
                                               x_u->values()[i] - x_k->values()[i], delta));
    }
    /**
     * only set the upper bound for the last 2*nCon entries (those are slack variables).
     * The lower bounds are initialized as 0
     */
    for (int i = 0; i < nlp_info_.nCon * 2; i++)
        qp_interface_->getUb()->setValueAt(nlp_info_.nVar + i, INF);


}


/**
 * This function sets up the object vector g of the QP problem
 * The (2*nCon+nVar) vector g_^T in QP problem will be the same as
 * [grad_f^T, rho* e^T], where the unit vector is of length (2*nCon).
 * @param grad      Gradient vector from nlp class
 * @param rho       Penalty Parameter
 */
void QPhandler::setup_g(shared_ptr<const Vector> grad, double rho) {
    qp_interface_->getG()->assign(1, grad->Dim(), grad->values());
    qp_interface_->getG()->assign_n(grad->Dim() + 1, nlp_info_.nCon * 2, rho);
}


/**
 * @brief This function updates the constraint if there is any changes to
 * the iterates
 */
void QPhandler::update_constraints(double delta, shared_ptr<const Vector> x_l,
                                   shared_ptr<const Vector> x_u,
                                   shared_ptr<const Vector> c_k,
                                   shared_ptr<const Vector> c_l,
                                   shared_ptr<const Vector> c_u,
                                   shared_ptr<const Vector> x_k) {
    for (int i = 0; i < nlp_info_.nCon; i++) {
        qp_interface_->getLbA()->setValueAt(i, c_l->values()[i] - c_k->values()[i]);
        qp_interface_->getUbA()->setValueAt(i, c_u->values()[i] - c_k->values()[i]);
    }
    for (int i = 0; i < nlp_info_.nVar; i++) {
        qp_interface_->getLb()->setValueAt(i, std::max(
                                               x_l->values()[i] - x_k->values()[i], -delta));
        qp_interface_->getUb()->setValueAt(i, std::min(
                                               x_u->values()[i] - x_k->values()[i], delta));
    }
}

/**
 * @brief This function updates the vector g in the QP subproblem when there are any
 * change to the values of penalty parameter
 *
 * @param rho               penalty parameter
 * @param nVar              number of variables in NLP
 */

void QPhandler::update_penalty(double rho) {
    qp_interface_->getG()->assign_n(nlp_info_.nVar + 1, nlp_info_.nCon * 2, rho);
}


/**
 * @brief This function updates the vector g in the QP subproblem when there are any
 * change to the values of gradient in NLP
 *
 * @param grad              the gradient vector from NLP
 */
void QPhandler::update_grad(shared_ptr<const Vector> grad) {
    qp_interface_->getG()->assign(1, grad->Dim(), grad->values());
}


/**
 * @brief Set up the H for the first time in the QP problem.
 * It will be concatenated as [H_k 0]
 *          		         [0   0]
 * where H_k is the Lagragian hessian evaluated at x_k and lambda_k.
 *
 * This method should only be called for once.
 *
 * @param hessian the Lagragian hessian evaluated at x_k and lambda_k from nlp
 * readers.
 */
void QPhandler::setup_H(shared_ptr<const SpTripletMat> hessian) {
    qp_interface_->getH()->setStructure(hessian);
    qp_interface_->getH()->setMatVal(hessian);

}

/**
 * @brief This function sets up the matrix A in the QP subproblem
 * The matrix A in QP problem will be concatenate as [J I -I]
 * @param jacobian  the Matrix object for Jacobian from c(x)
 */
void QPhandler::setup_A(shared_ptr<const SpTripletMat> jacobian) {

    I_info_A.irow1 = 1;
    I_info_A.irow2 = 1;
    I_info_A.jcol1 = nlp_info_.nVar + 1;
    I_info_A.jcol2 = nlp_info_.nVar + nlp_info_.nCon + 1;
    I_info_A.size = nlp_info_.nCon;

    qp_interface_->getA()->setStructure(jacobian, I_info_A);
    qp_interface_->getA()->setMatVal(jacobian->MatVal(), I_info_A);

}

/**
 *@brief Solve the QP with objective and constraints defined by its class members
 */
void QPhandler::solveQP(shared_ptr<SQPhotstart::Stats> stats,
                        shared_ptr<Options> options) {
    qp_interface_->optimizeQP(stats, options);
}

double QPhandler::GetObjective() {
    return qp_interface_->get_obj_value();
}

/**
 * This function initializes all class members which
 * will be used in _qp_interface
 *
 * @param nlp_info
 *          it contains the information about number
 *          of variables, number of constraints, number of
 *          elements in the Hessian and that of
 *          Jacobian. The definition of Index_info is
 *          in Types.hpp
 * @param Variable_type
 *          it specifies how the variables are bounded.
 *          Please check @ClassifyConstraintType in Algorithm.hpp
 *          for more details
 * @param qptype
 *          it specifies which type of QP is going to
 *          be solved. It can be either LP, or QP, or SOC
 */
void
QPhandler::allocate(SQPhotstart::Index_info nlp_info, SQPhotstart::QPType qptype) {
    qp_interface_ = make_shared<qpOASESInterface>(nlp_info, qptype);
}

void QPhandler::update_H(shared_ptr<const SpTripletMat> Hessian) {
    qp_interface_->getH()->setMatVal(Hessian);
}

void QPhandler::update_A(shared_ptr<const SpTripletMat> Jacobian) {
    qp_interface_->getA()->setMatVal(Jacobian->MatVal(), I_info_A);
}

const shared_ptr<qpOASESInterface>& QPhandler::getQpInterface() const {
    return qp_interface_;
}


} // namespace SQPhotstart



