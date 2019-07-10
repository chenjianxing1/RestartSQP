#include <sqphot/QPhandler.hpp>
#include <sqphot/SQPTNLP.hpp>
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
     *
     * setup the bounds for the QP subproblems according to the information from current iterate
     * 
     * @param delta 	 trust region radius	
     * @param x_k 	 current iterate point
     * @param c_k        current constraint value evaluated at x_k
     * @param x_l        the lower bounds for variables
     * @param x_u        the upper bounds for variables
     * @param c_l        the lower bounds for constraints
     * @param c_u        the upper bounds for constraints
     */
    bool QPhandler::setup_bounds(double delta, shared_ptr<const Vector> x_k,
                                 shared_ptr<const Vector> c_k,
                                 shared_ptr<const Vector> x_l, shared_ptr<const Vector> x_u,
                                 shared_ptr<const Vector> c_l, shared_ptr<const Vector> c_u) {
        int nCon = nlp_info_.nCon;
        int nVar = nlp_info_.nVar;
        lb_->assign_n(1, nVar, -delta);
        ub_->assign_n(1, nVar, delta);
        ub_->assign_n(nVar + 1, ub_->Dim() - nVar, INF);


        if (lbA_->Dim() == nCon) {
            lbA_->copy_vector(c_l->vector());
            ubA_->copy_vector(c_u->vector());
            lbA_->subtract_vector(c_k->vector());
            ubA_->subtract_vector(c_k->vector());
        } else {
            lbA_->assign(1, nCon, c_l->vector());
            lbA_->subtract_subvector(1, nCon, c_k->vector());
            lbA_->assign(nCon + 1, nVar, x_l->vector());
            lbA_->subtract_subvector(nCon + 1, nVar, x_k->vector());
            ubA_->assign(1, nCon, c_l->vector());
            ubA_->subtract_subvector(1, nCon, c_k->vector());
            ubA_->assign(nCon + 1, nVar, x_u->vector());
            ubA_->subtract_subvector(nCon + 1, nVar, x_k->vector());
        }
        return true;
    }


    /**
     * This function sets up the object vector g of the QP problem 
     *
     * @param grad 	Gradient vector from nlp class
     * @param rho  	Penalty Parameter
     */
    bool QPhandler::setup_g(shared_ptr<const Vector> grad, double rho) {
        g_->assign(1, grad->Dim(), grad->vector());
        g_->assign_n(grad->Dim() + 1, g_->Dim() - grad->Dim(), rho);
        return true;
    }


    /**
     * This function updates the bounds on x if there is any changes to the values of trust-region radius
     *
     * @param delta 	 trust region radius	
     * @param nVar 		 number of variables in NLP
     */
    bool QPhandler::update_bounds(double delta) {
        lb_->assign_n(1, nlp_info_.nVar, -delta);
        ub_->assign_n(1, nlp_info_.nVar, delta);
        return true;
    }

    /**
     * This function updates the vector g in the QP subproblem when there are any change to the values of penalty parameter 
     *
     * @param rho		penalty parameter
     * @param nVar 		number of variables in NLP
     */

    bool QPhandler::update_penalty(double rho) {
        g_->assign_n(nlp_info_.nVar + 1, g_->Dim() - nlp_info_.nVar, rho);
        return true;
    }


    /**
     * This function updates the vector g in the QP subproblem when there are any change to the values of gradient in NLP
     *
     * @param grad		the gradient vector from NLP
     */
    bool QPhandler::update_grad(shared_ptr<const Vector> grad) {
        g_->assign(1, grad->Dim(), grad->vector());
        return true;
    }


    /**
     *This function sets up the Hessian used for the QP problem
     *
     * @param hessian 	the Matrix object for Hessian from NLP
     */
    bool QPhandler::setup_H(shared_ptr<const Matrix> hessian) {
        H_->copyMatrix(hessian);
        return true;
    }

    /**
     * This function sets up the matrix A in the QP subproblem
     *
     * @param jacobian 	the Matrix object for Jacobian from c(x)
     */
    bool QPhandler::setup_A(shared_ptr<const Matrix> jacobian) {
        // A_->copyMatrix(jacobian);
        // int nCon = jacobian->RowNum();
        // int nVar = jacobian->ColNum();
        //// FIXME: assign dimension in constructor
        ////  A_->assign_dimension(nCon, nVar + 2 * nCon);
        // if (A_->current_I_num() != A_->num_I()) {
        //     A_->add_I(1, nVar + 1, nCon, true);
        //     A_->add_I(1, nVar + nCon + 1, nCon, false);
        // }

        return true;
    }

    /**
     *
     *
     * @param stats	the static used to record iterations numbers
     */
    bool QPhandler::solveQP(shared_ptr<SQPhotstart::Stats> stats, shared_ptr<Options> options) {
        qp_interface_->optimizeQP(H_, g_, A_, lbA_, ubA_, lb_, ub_, stats, options);
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
        int numVar_QP;
        int numCon_QP;

         numVar_QP = nlp_info.nVar + 2 * nlp_info.nCon;
         numCon_QP = nlp_info.nCon;
        // A_ = make_shared<Matrix>(nlp_info.nnz_jac_g, 2);

        if (qptype == QP) {
            H_ = make_shared<Matrix>(nlp_info.nnz_h_lag, numVar_QP, numVar_QP);
        }

        qp_interface_ = std::make_shared<qpOASESInterface>(numVar_QP, numCon_QP);
        lbA_ = make_shared<Vector>(numCon_QP);
        ubA_ = make_shared<Vector>(numCon_QP);
        lb_ = make_shared<Vector>(numVar_QP);
        ub_ = make_shared<Vector>(numVar_QP);
        g_ = make_shared<Vector>(numVar_QP);
        return true;
    }

    bool QPhandler::update_H(shared_ptr<const Matrix> Hessian) {
        return false;
    }

    bool QPhandler::update_A(shared_ptr<const Matrix> Jacobian) {
        return false;
    }


} // namespace SQPhotstart


