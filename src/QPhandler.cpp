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
    QPhandler::QPhandler() {
        //message for debugging purpose
        std::cout << "QP is created" << std::endl;
    }

    /**
     *Default destructor
     */
    QPhandler::~QPhandler() {
        std::cout << "QP is destroyed" << std::endl;
    }

    /**
     * Get the optimal solution from the QPhandler_interface
     *
     *This is only an interface for user to avoid call interface directly.
     * @param p_k       the pointer to an empty array with the length equal to the size
     * of the QP subproblem
     */
    bool QPhandler::GetOptimalSolution(double* p_k) {
        qp_interface_->get_optimal_solution(p_k);

        //        print_("p_k", p_k,nlp_info_.nVar+2*nlp_info_.nCon);
        return true;
    }

    /**
     *Get the multipliers from the QPhandler_interface
     *
     *This is only an interface for user to avoid call interface directly.
     *
     * @param y_k       the pointer to an empty array with the length equal to the size of multipliers of the QP subproblem
     */
    bool QPhandler::GetMultipliers(double* y_k) {
        qp_interface_->get_multipliers(y_k);
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
    bool QPhandler::init(Index_info nlp_info, QPType qptype) {
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
    bool QPhandler::setup_bounds(double delta, shared_ptr<const Vector> x_k,
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

        //        if(DEBUG){
        //            std::cout<<"the value of lb is"<<endl;
        //        qp_interface_->getLb()->print();
        //            divider();
        //                        std::cout<<"the value of ub is"<<endl;
        //        qp_interface_->getUb()->print();
        //            divider();
        //        }
        //
        return true;
    }


    /**
     * This function sets up the object vector g of the QP problem
     * The (2*nCon+nVar) vector g_^T in QP problem will be the same as
     * [grad_f^T, rho* e^T], where the unit vector is of length (2*nCon).
     * @param grad      Gradient vector from nlp class
     * @param rho       Penalty Parameter
     */
    bool QPhandler::setup_g(shared_ptr<const Vector> grad, double rho) {
        qp_interface_->getG()->assign(1, grad->Dim(), grad->values());
        qp_interface_->getG()->assign_n(grad->Dim() + 1, nlp_info_.nCon * 2, rho);
        //        std::cout<<"--------------------------"<<std::endl;
        //        qp_interface_->getG()->print();
        //        std::cout<<"--------------------------"<<std::endl;
        return true;
    }


    /**
     * @brief This function updates the constraint if there is any changes to
     * the iterates
     */
    bool
    QPhandler::update_constraints(double delta, shared_ptr<const Vector> x_l,
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
        return true;
    }

    /**
     * @brief This function updates the vector g in the QP subproblem when there are any
     * change to the values of penalty parameter
     *
     * @param rho               penalty parameter
     * @param nVar              number of variables in NLP
     */

    bool QPhandler::update_penalty(double rho) {
        qp_interface_->getG()->assign_n(nlp_info_.nVar + 1, nlp_info_.nCon * 2, rho);
        return true;
    }


    /**
     * @brief This function updates the vector g in the QP subproblem when there are any
     * change to the values of gradient in NLP
     *
     * @param grad              the gradient vector from NLP
     */
    bool QPhandler::update_grad(shared_ptr<const Vector> grad) {
        qp_interface_->getG()->assign(1, grad->Dim(), grad->values());
        return true;
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
    bool QPhandler::setup_H(shared_ptr<const SpTripletMat> hessian) {
        qp_interface_->getH()->setStructure(hessian);
        qp_interface_->getH()->setMatVal(hessian->MatVal());

        return true;
    }

    /**
     * @brief This function sets up the matrix A in the QP subproblem
     * The matrix A in QP problem will be concatenate as [J I -I]
     * @param jacobian  the Matrix object for Jacobian from c(x)
     */
    bool QPhandler::setup_A(shared_ptr<const SpTripletMat> jacobian) {

        I_info_A.irow1 = 1;
        I_info_A.irow2 = 1;
        I_info_A.jcol1 = nlp_info_.nVar + 1;
        I_info_A.jcol2 = nlp_info_.nVar + nlp_info_.nCon + 1;
        I_info_A.size = nlp_info_.nCon;

        qp_interface_->getA()->setStructure(jacobian, I_info_A);
        qp_interface_->getA()->setMatVal(jacobian->MatVal(), I_info_A);

        return true;
    }

    /**
     *@brief Solve the QP with objective and constraints defined by its class members
     */
    bool QPhandler::solveQP(shared_ptr<SQPhotstart::Stats> stats,
                            shared_ptr<Options> options) {
        qp_interface_->optimizeQP(stats, options);
        assert(qp_interface_->get_status() == OPTIMAL);
        return true;
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
    bool
    QPhandler::allocate(SQPhotstart::Index_info nlp_info, SQPhotstart::QPType qptype) {
        qp_interface_ = make_shared<qpOASESInterface>(nlp_info, qptype);
        return true;
    }

    bool QPhandler::update_H(shared_ptr<const SpTripletMat> Hessian) {
        qp_interface_->getH()->setMatVal(Hessian->MatVal(), I_info_H);
        return true;
    }

    bool QPhandler::update_A(shared_ptr<const SpTripletMat> Jacobian) {
        qp_interface_->getA()->setMatVal(Jacobian->MatVal(), I_info_A);
        return true;
    }

    const shared_ptr<qpOASESInterface>& QPhandler::getQpInterface() const {
        return qp_interface_;
    }


} // namespace SQPhotstart



