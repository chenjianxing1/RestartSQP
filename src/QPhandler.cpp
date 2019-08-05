/* Copyright (C) 2019
 * All Rights Reserved.
 *
 * Authors: Xinyi Luo
 * Date:2019-07
 */
#include <sqphot/QPhandler.hpp>


namespace SQPhotstart {
    QPhandler::QPhandler(Index_info nlp_info) :
            nlp_info_(nlp_info) {
        solverInterface_ = make_shared<qpOASESInterface>(nlp_info, QP);

    }

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
        solverInterface_->get_optimal_solution(p_k);

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
        solverInterface_->get_multipliers(y_k);

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
    void QPhandler::set_bounds(double delta, shared_ptr<const Vector> x_k,
                               shared_ptr<const Vector> x_l,
                               shared_ptr<const Vector> x_u,
                               shared_ptr<const Vector> c_k,
                               shared_ptr<const Vector> c_l,
                               shared_ptr<const Vector> c_u) {

        for (int i = 0; i < nlp_info_.nCon; i++) {
            solverInterface_->set_lbA(i, c_l->values()[i] - c_k->values()[i]);
            solverInterface_->set_ubA(i, c_u->values()[i] - c_k->values()[i]);
        }

        for (int i = 0; i < nlp_info_.nVar; i++) {
            solverInterface_->set_lb(i, std::max(
                    x_l->values()[i] - x_k->values()[i], -delta));
            solverInterface_->set_ub(i, std::min(
                    x_u->values()[i] - x_k->values()[i], delta));
        }
        /**
         * only set the upper bound for the last 2*nCon entries (those are slack variables).
         * The lower bounds are initialized as 0
         */
        for (int i = 0; i < nlp_info_.nCon * 2; i++)
            solverInterface_->set_ub(nlp_info_.nVar + i, INF);

    }


/**
 * This function sets up the object vector g of the QP problem
 * The (2*nCon+nVar) vector g_^T in QP problem will be the same as
 * [grad_f^T, rho* e^T], where the unit vector is of length (2*nCon).
 * @param grad      Gradient vector from nlp class
 * @param rho       Penalty Parameter
 */
    void QPhandler::set_g(shared_ptr<const Vector> grad, double rho) {
        for (int i = 0; i < nlp_info_.nVar + 2 * nlp_info_.nCon; i++)
            if (i < nlp_info_.nVar)
                solverInterface_->set_g(i, grad->values()[i]);

            else
                solverInterface_->set_g(i, rho);


    }


/**
 * @brief This function updates the constraint if there is any changes to
 * the iterates
 */
    void QPhandler::update_bounds(double delta, shared_ptr<const Vector> x_l,
                                  shared_ptr<const Vector> x_u,
                                  shared_ptr<const Vector> x_k,
                                  shared_ptr<const Vector> c_l,
                                  shared_ptr<const Vector> c_u,
                                  shared_ptr<const Vector> c_k) {
        for (int i = 0; i < nlp_info_.nCon; i++) {
            solverInterface_->set_lbA(i, c_l->values()[i] - c_k->values()[i]);
            solverInterface_->set_ubA(i, c_u->values()[i] - c_k->values()[i]);
        }

        for (int i = 0; i < nlp_info_.nVar; i++) {
            solverInterface_->set_lb(i, std::max(
                    x_l->values()[i] - x_k->values()[i], -delta));
            solverInterface_->set_ub(i, std::min(
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
        for (int i = nlp_info_.nVar; i < nlp_info_.nVar + nlp_info_.nCon * 2; i++)
            solverInterface_->set_g(i, rho);

    }


/**
 * @brief This function updates the vector g in the QP subproblem when there are any
 * change to the values of gradient in NLP
 *
 * @param grad              the gradient vector from NLP
 */
    void QPhandler::update_grad(shared_ptr<const Vector> grad) {
        for (int i = 0; i < nlp_info_.nVar; i++)
            solverInterface_->set_g(i, grad->values()[i]);
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
    void QPhandler::set_H(shared_ptr<const SpTripletMat> hessian) {
        solverInterface_->set_H_structure(hessian);
        solverInterface_->set_H_values(hessian);

    }

/**
 * @brief This function sets up the matrix A in the QP subproblem
 * The matrix A in QP problem will be concatenate as [J I -I]
 * @param jacobian  the Matrix object for Jacobian from c(x)
 */
    void QPhandler::set_A(shared_ptr<const SpTripletMat> jacobian) {

        I_info_A.irow1 = 1;
        I_info_A.irow2 = 1;
        I_info_A.jcol1 = nlp_info_.nVar + 1;
        I_info_A.jcol2 = nlp_info_.nVar + nlp_info_.nCon + 1;
        I_info_A.size = nlp_info_.nCon;

        solverInterface_->set_A_structure(jacobian, I_info_A);
        solverInterface_->set_A_values(jacobian, I_info_A);

    }

/**
 *@brief Solve the QP with objective and constraints defined by its class members
 */
    void QPhandler::solveQP(shared_ptr<SQPhotstart::Stats> stats,
                            shared_ptr<Options> options) {
        solverInterface_->optimizeQP(stats, options);
    }

    double QPhandler::GetObjective() {
        return solverInterface_->get_obj_value();
    }


    void QPhandler::update_H(shared_ptr<const SpTripletMat> Hessian) {
        solverInterface_->set_H_values(Hessian);
    }

    void QPhandler::update_A(shared_ptr<const SpTripletMat> Jacobian) {
        solverInterface_->set_A_values(Jacobian, I_info_A);
    }

    const shared_ptr<QPSolverInterface>& QPhandler::getQpInterface() const {
        return solverInterface_;
    }

    void QPhandler::update_delta(double delta, shared_ptr<const Vector> x_l,
                                 shared_ptr<const Vector> x_u,
                                 shared_ptr<const Vector> x_k) {

        for (int i = 0; i < nlp_info_.nVar; i++) {
            solverInterface_->set_lb(i, std::max(
                    x_l->values()[i] - x_k->values()[i], -delta));
            solverInterface_->set_ub(i, std::min(
                    x_u->values()[i] - x_k->values()[i], delta));
        }

    }


} // namespace SQPhotstart



