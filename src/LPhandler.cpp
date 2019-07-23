/* Copyright (C) 2019
 * All Rights Reserved.
 *
 * Authors: Xinyi Luo
 * Date:2019-07
 */

#include <sqphot/LPhandler.hpp>

namespace SQPhotstart {
    /**
     * Default constructor
     *
     */
    LPhandler::LPhandler() {
        //message for debugging purpose
        std::cout << "LP is created" << std::endl;
    }

    /**
     *Default destructor
     */
    LPhandler::~LPhandler() {
        std::cout << "LP is destroyed" << std::endl;
    }

    /**
     * Get the optimal solution from the LPhandler_interface
     *
     *This is only an interface for user to avoid call interface directly.
     * @param p_k       the pointer to an empty array with the length equal to the size of the QP subproblem
     */
    bool LPhandler::GetOptimalSolution(double *p_k) {
        lp_interface_->get_optimal_solution(p_k);

        //        print_("p_k", p_k,nlp_info_.nVar+2*nlp_info_.nCon);
        return true;
    }


    /**
     * This function initializes all objects will be used in this class.
     * (Probably more functionality can be added to this)
     *
     * @param nlp_info
     *          it contains the information about number of variables, number of constraints, number of
     *          elements in the Hessian and that of Jacobian. The definition of Index_info is in Types.hpp
     * @param Constraint_type
     *          it specifies how the variables are bounded. Please check @ClassifyConstraintType in Algorithm.hpp
     *          for more details
     * @param qptype
     *          it specifies which type of QP is going to be solved. It can be either LP, or QP, or SOC
     */
    bool LPhandler::init(Index_info nlp_info, QPType qptype) {
        assert(qptype == LP);
        allocate(nlp_info, qptype);
        nlp_info_ = nlp_info;
        qptype_ = qptype;

        return true;
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
    bool LPhandler::setup_bounds(double delta, shared_ptr<const Vector> x_k,
                                 shared_ptr<const Vector> x_l,
                                 shared_ptr<const Vector> x_u,
                                 shared_ptr<const Vector> c_k,
                                 shared_ptr<const Vector> c_l,
                                 shared_ptr<const Vector> c_u) {

        for (int i = 0; i < nlp_info_.nCon; i++) {
            lp_interface_->getLbA()->setValueAt(i, c_l->values()[i] - c_k->values()[i]);
            lp_interface_->getUbA()->setValueAt(i, c_u->values()[i] - c_k->values()[i]);
        }

        for (int i = 0; i < nlp_info_.nVar; i++) {
            lp_interface_->getLb()->setValueAt(i, std::max(
                    x_l->values()[i] - x_k->values()[i], -delta));
            lp_interface_->getUb()->setValueAt(i, std::min(
                    x_u->values()[i] - x_k->values()[i], delta));
        }
        /**
         * only set the upper bound for the last 2*nCon entries (those are slack variables).
         * The lower bounds are initialized as 0
         */
        for (int i = 0; i < nlp_info_.nCon * 2; i++)
            lp_interface_->getUb()->setValueAt(nlp_info_.nVar + i, INF);
        return true;
    }



    /**
     * This function updates the bounds on x if there is any changes to the values of trust-region radius
     *
     * @param delta      trust region radius
     * @param nVar               number of variables in NLP
     */
    bool
    LPhandler::update_constraints(double delta, shared_ptr<const Vector> x_l,
                                  shared_ptr<const Vector> x_u,
                                  shared_ptr<const Vector> c_k,
                                  shared_ptr<const Vector> c_l,
                                  shared_ptr<const Vector> c_u,
                                  shared_ptr<const Vector> x_k) {
        for (int i = 0; i < nlp_info_.nCon; i++) {
            lp_interface_->getLbA()->setValueAt(i, c_l->values()[i] - c_k->values()[i]);
            lp_interface_->getUbA()->setValueAt(i, c_u->values()[i] - c_k->values()[i]);
        }
        for (int i = 0; i < nlp_info_.nVar; i++) {
            lp_interface_->getLb()->setValueAt(i, std::max(
                    x_l->values()[i] - x_k->values()[i], -delta));
            lp_interface_->getUb()->setValueAt(i, std::min(
                    x_u->values()[i] - x_k->values()[i], delta));
        }
        return true;
    }

    /**
     * This function updates the vector g in the QP subproblem when there are any change to the values of penalty parameter
     *
     * @param rho               penalty parameter
     * @param nVar              number of variables in NLP
     */

    bool LPhandler::update_penalty(double rho) {
        lp_interface_->getG()->assign_n(nlp_info_.nVar + 1, nlp_info_.nCon * 2, rho);
        return true;
    }


    /**
     * This function updates the vector g in the QP subproblem when there are any change to the values of gradient in NLP
     *
     * @param grad              the gradient vector from NLP
     */
    bool LPhandler::update_grad(shared_ptr<const Vector> grad) {
        lp_interface_->getG()->assign(1, grad->Dim(), grad->values());
        return true;
    }


    /**
     * @brief This function sets up the matrix A in the LP subproblem
     * The matrix A in QP problem will be concatenate as [J I -I]
     * @param jacobian  the Matrix object for Jacobian from c(x)
     */
    bool LPhandler::setup_A(shared_ptr<const SpTripletMat> jacobian) {


        I_info_A.irow1 = 1;
        I_info_A.irow2 = 1;
        I_info_A.jcol1 = nlp_info_.nVar + 1;
        I_info_A.jcol2 = nlp_info_.nVar + nlp_info_.nCon + 1;
        I_info_A.size = nlp_info_.nCon;

        lp_interface_->getA()->setStructure(jacobian, I_info_A);
        lp_interface_->getA()->setMatVal(jacobian->MatVal(), I_info_A);

        return true;
    }


    /**
     * @brief solve the LP problem with bounds and objectives defined by the class
     * members
     */

    bool LPhandler::solveLP(shared_ptr<SQPhotstart::Stats> stats,
                            shared_ptr<Options> options) {
        lp_interface_->optimizeLP(stats, options);
        lp_obj_ = lp_interface_->get_obj_value();
        return true;
    }


    /**
     * @brief This function initializes all class members which
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
    bool
    LPhandler::allocate(SQPhotstart::Index_info nlp_info, SQPhotstart::QPType qptype) {
        lp_interface_ = make_shared<qpOASESInterface>(nlp_info, qptype);
        return true;
    }


    /**
     * @brief
     * @param Jacobian
     * @return
     */
    bool LPhandler::update_A(shared_ptr<const SpTripletMat> Jacobian) {
        lp_interface_->getA()->setMatVal(Jacobian->MatVal(), I_info_A);
        return true;
    }

    /**
     * @brief
     *
     * @param rhs
     * @return
     */
    bool LPhandler::copy_QP_info(shared_ptr<const QPhandler> rhs) {
        lp_interface_->getLbA()->copy_vector(rhs->getQpInterface()->getLbA());
        lp_interface_->getLb()->copy_vector(rhs->getQpInterface()->getLb());
        lp_interface_->getUb()->copy_vector(rhs->getQpInterface()->getUb());
        lp_interface_->getUbA()->copy_vector(rhs->getQpInterface()->getUbA());

        lp_interface_->getG()->assign(nlp_info_.nVar + 1, 2 * nlp_info_.nCon,
                                      rhs->getQpInterface()
                                              ->getG()->values() + nlp_info_.nVar);

        lp_interface_->getA()->copy(rhs->getQpInterface()->getA());

        return true;
    }

    /**
     * @brief This function sets up the object vector g of the LP problem
     * The (2*nCon+nVar) vector g_^T in QP problem will be the same as
     * [0 , rho* e^T], where the unit vector is of length (2*nCon).
     * @param rho       Penalty Parameter
     */
    bool LPhandler::setup_g(double rho) {
        lp_interface_->getG()->assign_n(nlp_info_.nVar+ 1, nlp_info_.nCon * 2, rho);
        return true;
    }


} // namespace SQPhotstart



