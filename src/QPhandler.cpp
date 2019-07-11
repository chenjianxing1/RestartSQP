#include <sqphot/QPhandler.hpp>

namespace SQPhotstart {
    /**
     * Default constructor
     *
     */
    QPhandler::QPhandler() {
        std::cout << "QP is created" << std::endl;
    }

    /** 
     *Default destructor
     */
    QPhandler::~QPhandler() {
        std::cout << "QP is destroyed" << std::endl;
    }

    /**
     * Get the optimal solution from the QPsolver_interface
     *
     *This is only an interface for user to avoid call interface directly.
     * @param p_k 	the pointer to an empty array with the length equal to the size of the QP subproblem
     */
    bool QPhandler::GetOptimalSolution(double* p_k) {
        qp_interface_->get_optimal_solution(p_k);
        return true;
    }

    /**
     *Get the multipliers from the QPsolver_interface
     *
     *This is only an interface for user to avoid call interface directly.
     *
     * @param y_k 	the pointer to an empty array with the length equal to the size of multipliers of the QP subproblem
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
     * 		it contains the information about number of variables, number of constraints, number of 
     * 		elements in the Hessian and that of Jacobian. The definition of Index_info is in Types.hpp
     * @param Constraint_type
     * 		it specifies how the variables are bounded. Please check @ClassifyConstraintType in Algorithm.hpp
     * 		for more details
     * @param qptype
     * 		it specifies which type of QP is going to be solved. It can be either LP, or QP, or SOC
     */
    bool QPhandler::init(Index_info nlp_info, QPType qptype) {
        allocate(nlp_info, qptype);
        nlp_info_ = nlp_info;
        qptype_ = qptype;

        return true;
    }

    /**
     * Setup the bounds for the QP subproblems according to the information from current iterate
     * The bound is formulated as
     *   max(-delta, x_l-x_k)<=p<=min(delta,x_u-x_k)
     * and  u,v>=0
     * @param delta 	 trust region radius	
     * @param x_k 	 current iterate point
     * @param c_k        current constraint value evaluated at x_k
     * @param x_l        the lower bounds for variables
     * @param x_u        the upper bounds for variables
     * @param c_l        the lower bounds for constraints
     * @param c_u        the upper bounds for constraints
     */
    bool QPhandler::setup_bounds(double delta, shared_ptr<const Vector> x_k, shared_ptr<const Vector> x_l,
                                 shared_ptr<const Vector> x_u) {
        for (int i = 0; i < nlp_info_.nVar; i++) {
            qp_interface_->getLb()->setValueAt(i, std::max(x_l->getValueAt(i) - x_k->getValueAt(i), -delta));
            qp_interface_->getUb()->setValueAt(i, std::min(x_u->getValueAt(i) - x_k->getValueAt(i), delta));
        }
        /**
         * only set the upper bound for the last 2*nCon entries (those are slack variables).
         * The lower bounds are initialized as 0
         */
        for (int i = 0; i < nlp_info_.nCon * 2; i++)
            qp_interface_->getUb()->setValueAt(nlp_info_.nVar + i, INF);

        return true;
    }


    /**
     * This function sets up the object vector g of the QP problem
     * The (2*nCon+nVar) vector g_^T in QP problem will be the same as
     * [grad_f^T, rho* e^T], where the unit vector is of length (2*nCon).
     * @param grad 	Gradient vector from nlp class
     * @param rho  	Penalty Parameter
     */
    bool QPhandler::setup_g(shared_ptr<const Vector> grad, double rho) {
        qp_interface_->getG()->assign(1, grad->Dim(), grad->values());
        qp_interface_->getG()->assign_n(grad->Dim() + 1, nlp_info_.nCon * 2, rho);
        return true;
    }


    /**
     * This function updates the bounds on x if there is any changes to the values of trust-region radius
     *
     * @param delta 	 trust region radius	
     * @param nVar 		 number of variables in NLP
     */
    bool QPhandler::update_bounds(double delta, shared_ptr<const Vector> x_k, shared_ptr<const Vector> x_l,
                                  shared_ptr<const Vector> x_u) {
        for (int i = 0; i < nlp_info_.nVar; i++) {
            qp_interface_->getLb()->setValueAt(i, std::max(x_l->getValueAt(i) - x_k->getValueAt(i), -delta));
            qp_interface_->getUb()->setValueAt(i, std::min(x_u->getValueAt(i) - x_k->getValueAt(i), delta));
        }
        return true;
    }

    /**
     * This function updates the vector g in the QP subproblem when there are any change to the values of penalty parameter 
     *
     * @param rho		penalty parameter
     * @param nVar 		number of variables in NLP
     */

    bool QPhandler::update_penalty(double rho) {
        qp_interface_->getG()->assign_n(nlp_info_.nVar + 1, nlp_info_.nCon * 2, rho);
        return true;
    }


    /**
     * This function updates the vector g in the QP subproblem when there are any change to the values of gradient in NLP
     *
     * @param grad		the gradient vector from NLP
     */
    bool QPhandler::update_grad(shared_ptr<const Vector> grad) {
        qp_interface_->getG()->assign(1, grad->Dim(), grad->values());
        return true;
    }


    /**
     *This function sets up the Hessian used for the QP problem
     *
     * @param hessian 	the Matrix object for Hessian from NLP
     */
    bool QPhandler::setup_H(shared_ptr<const SpMatrix> hessian) {
        Identity2Info I_info_H;
        qp_interface_->getH()->setStructure(hessian,I_info_H);
        qp_interface_->getH()->setMatVal(hessian->MatVal());
        return true;
    }

    /**
     * This function sets up the matrix A in the QP subproblem
     * The matrix A in QP problem will be concatenate as [J I -I]
     * @param jacobian 	the Matrix object for Jacobian from c(x)
     */
    bool QPhandler::setup_A(shared_ptr<const SpMatrix> jacobian) {

        Identity2Info I_info_A;
        I_info_A.irow1 = 1;
        I_info_A.irow2 = 1;
        I_info_A.jcol1 = nlp_info_.nVar+1;
        I_info_A.jcol2 = nlp_info_.nVar+nlp_info_.nCon+1;
        I_info_A.size = nlp_info_.nCon;
        
        qp_interface_->getA()->setStructure(jacobian,I_info_A);
        qp_interface_->getA()->setMatVal(jacobian->MatVal());
        return true;
    }

    /**
     *
     *
     * @param stats	the static used to record iterations numbers
     */
    bool QPhandler::solveQP(shared_ptr<SQPhotstart::Stats> stats, shared_ptr<Options> options) {
//        qp_interface_->optimizeQP(H_qpOASES_, g_, A_qpOASES_, lbA_, ubA_, lb_, ub_, stats, options);
        qp_obj_ = qp_interface_->get_obj_value();
        return true;
    }


    /**
     *This function initializes all class members which will be used in _qp_interface
     *
     * @param nlp_info 
     * 		it contains the information about number of variables, number of constraints, number of 
     * 		elements in the Hessian and that of Jacobian. The definition of Index_info is in Types.hpp
     * @param Variable_type
     * 		it specifies how the variables are bounded. Please check @ClassifyConstraintType in Algorithm.hpp
     * 		for more details
     * @param qptype
     * 		it specifies which type of QP is going to be solved. It can be either LP, or QP, or SOC
     */
    bool QPhandler::allocate(SQPhotstart::Index_info nlp_info, SQPhotstart::QPType qptype) {
        qp_interface_ = make_shared<qpOASESInterface>(nlp_info);
        return true;
    }

    bool QPhandler::update_H(shared_ptr<const SpMatrix> Hessian) {
        qp_interface_->getH()->setMatVal(Hessian->MatVal());
        return true;
    }

    bool QPhandler::update_A(shared_ptr<const SpMatrix> Jacobian) {
        qp_interface_->getA()->setMatVal(Jacobian->MatVal());
        return true;
    }


} // namespace SQPhotstart


