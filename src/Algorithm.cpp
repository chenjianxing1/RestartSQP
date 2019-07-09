#include <sqphot/Algorithm.hpp>

namespace SQPhotstart {

    /**
     * Default Constructor
     */
    Algorithm::Algorithm() {
        cout << "ALG created" << endl;
    }

    /**
     * Default Destructor
     */
    Algorithm::~Algorithm() {
        delete[] cons_type_;
        cons_type_ = NULL;
        delete[] bound_cons_type_;
        bound_cons_type_ = NULL;
        cout << "ALG destroyed" << endl;
    }

    /**
     * This is the main function to optimize the NLP given as the input
     *
     * @param nlp: the nlp reader that read data of the function to be minimized;
     */
    bool Algorithm::Optimize(SmartPtr <Ipopt::TNLP> nlp) {
        allocate(nlp);
        Presolve();

        /** Main iteration */

        // while (stats->iter < options->iter_max) {
        /* setup the QPs and solve them */
        //setupQP();
        //based on the information given by NLP reader; otherwise, it will do nothing
        //myQP->solveQP(stats, options);//solve the SL1QP problems
        //get_search_direction(myQP);//
        //norm_p_k_ = p_k_->getInfNorm(); //calculate the infinity norm of the search direction

        /* Update the penalty parameter if necessary*/
        //penalty_update();
        //
        //     /* Calculate the the trial points, x_trial = x_k+p_k, */
        //     x_trial_->copy_vector(x_k_->vector());
        //     x_trial_->add_vector(p_k_->vector());
        //
        //     /** Calculate f_trial, c_trial and infea_measure_trial for the trial points x_trial*/
        //     nlp_->Eval_f(x_trial_, obj_value_trial_);
        //     nlp_->Eval_constraints(x_trial_, c_trial_);
        //
        //     infea_cal(true);
        //     ratio_test();
        //
        //     /* Calculate the second-order-correction steps*/
        //     second_order_correction();
        //
        //     /* Update the radius and the QP bounds if the radius has been changed*/
        //     radius_update();
        //
        //     stats->iter_addone();
        //
        //     /* output some information to the console, TODO: change it*/
        //     if (stats->iter % 10 == 0)log->print_header();
        //     if (options->printLevel > 0)
        //         log->print_main_iter(stats->iter, obj_value_, norm_p_k_, infea_measure_, delta_, rho_);
        //
        //     /**check if the current iterates is optimal and decide to exit the loop or not*/
        //     termination_check();
        //     if (exitflag_ != UNKNOWN) {
        //         break;
        //     }
        // }
        // /** print the final summary message to the console*/
        // log->print_final(stats->iter, stats->qp_iter, obj_value_, norm_p_k_, infea_measure_, exitflag_);
        return true;
    }

    /**
     *
     *This is the function that checks if the current point is optimal, and decides if to exit the loop or not
     *
     *@return
     * 	 if it decides the function is optimal, the class member _exitflag = OPTIMAL
     * 	 if it decides that there is an error during the function run or the function cannot be solved, it will assign _exitflag the
     * 	 	corresponding code according to the error type.
     */
    bool Algorithm::termination_check() {
        if (norm_p_k_ < options->tol) {
            exitflag_ = OPTIMAL;
        }
        return true;
    }


    /**
     * This function initializes the objects required by the SQP Algorithm, copies some parameters required by the algorithm,
     * obtains the function information for the first QP, and solve the first QP and LP.
     *
     */
    bool Algorithm::Presolve() {
        nlp_->Get_bounds_info(x_l_, x_u_, c_l_, c_u_);
        nlp_->Get_starting_point(x_k_, lambda_);
        nlp_->shift_starting_point(x_k_, x_l_, x_u_);
        nlp_->Eval_f(x_k_, obj_value_);
        nlp_->Eval_gradient(x_k_, grad_f_);
        nlp_->Eval_constraints(x_k_, c_k_);
	nlp_->Get_Structure_Hessian(x_k_, lambda_, hessian_);
        nlp_->Eval_Hessian(x_k_, lambda_, hessian_);
	nlp_->Get_Strucutre_Jacobian(x_k_,jacobian_);
        nlp_->Eval_Jacobian(x_k_, jacobian_);
        ClassifyConstraintType();

        infea_cal(false);
        /**set up the QP objects*/
        myQP->init(nlp_->nlp_info_, QP);
        myLP->init(nlp_->nlp_info_, LP);
        log->print_header();
	log->print_main_iter(stats->iter, obj_value_, 0.0, infea_measure_, delta_, rho_);
        return true;
    }

    /**
     * This function initializes all the shared pointer which will be used in the Algorithm::Optimize, and it copies all parameters
     * that might be changed during the run of the function Algorithm::Optimize
     *
     * @param nlp: the nlp reader that read data of the function to be minimized;
     */
    bool Algorithm::allocate(SmartPtr <Ipopt::TNLP> nlp) {
        nlp_ = make_shared<SQPTNLP>(nlp);
        nVar_ = nlp_->nlp_info_.nVar;
        nCon_ = nlp_->nlp_info_.nCon;
        cons_type_ = new ConstraintType[nCon_];
        bound_cons_type_ = new ConstraintType[nVar_];
        x_k_ = make_shared<Vector>(nVar_);
        x_trial_ = make_shared<Vector>(nVar_);
        p_k_ = make_shared<Vector>(nVar_);
	lambda_ = make_shared<Vector>(nCon_);
        c_k_ = make_shared<Vector>(nCon_);
        c_trial_ = make_shared<Vector>(nCon_);
        x_l_ = make_shared<Vector>(nVar_);
        x_u_ = make_shared<Vector>(nVar_);
        c_l_ = make_shared<Vector>(nCon_);
        c_u_ = make_shared<Vector>(nCon_);
        grad_f_ = make_shared<Vector>(nVar_);
        jacobian_ = make_shared<Matrix>(nlp_->nlp_info_.nnz_jac_g,nCon_,nVar_);
        hessian_ = make_shared<Matrix>(nlp_->nlp_info_.nnz_h_lag, nVar_,nVar_);

        options = make_shared<Options>();
        stats = make_shared<Stats>();
        log = make_shared<Log>();
        myQP = make_shared<QPhandler>();
        myLP = make_shared<QPhandler>();

        delta_ = options->delta;
        rho_ = options->rho;


        return true;
    }


    /**
     *
     * This function calculates the infeasibility measure for either tiral point or current iterate
     *
     *@param trial: true if the user are going to evaluate the infeasibility measure of the trial point _x_trial;
     *		        infea_fun_trial = norm(-max(c_trial-cu,0),1)+norm(-min(c_trial-cl,0),1);
     *		        infea_var_trial = norm(-min(x_trial-bl,0),1)+norm(-max(x_trial-bu,0),1);
     *		        infea_measure_trial = infea_fun_trial+ infea_var_trial;
     *
     * 	            false if the user are going to evaluate the infeasibility measure of the current iterates _x_k
     *	    		infea_fun = norm(-max(c-cu,0),1)+norm(-min(c-cl,0),1);
     *  	    		infea_var = norm(-min(x-bl,0),1)+norm(-max(x-bu,0),1);
     *  		    	infea_measure = infea_fun+infea_var;
     *
     */
    bool Algorithm::infea_cal(bool trial) {
        Number infea_con = 0.0;
        Number infea_var = 0.0;

        if (trial) {
            for (int i = 0; i < x_k_->Dim(); i++) {
                if (x_trial_->getEntryAt(i) < x_l_->getEntryAt(i))
                    infea_var += (x_l_->getEntryAt(i) - x_trial_->getEntryAt(i));
                else if (x_trial_->getEntryAt(i) > x_u_->getEntryAt(i))
                    infea_var += (x_trial_->getEntryAt(i) - x_u_->getEntryAt(i));
            }
            for (int i = 0; i < c_k_->Dim(); i++) {
                if (c_trial_->getEntryAt(i) < c_l_->getEntryAt(i))
                    infea_con += (c_l_->getEntryAt(i) - c_trial_->getEntryAt(i));
                else if (c_trial_->getEntryAt(i) > c_u_->getEntryAt(i))
                    infea_con += (c_trial_->getEntryAt(i) - c_u_->getEntryAt(i));
            }
            infea_measure_trial_ = infea_con + infea_var;
        } else {
            for (int i = 0; i < x_k_->Dim(); i++) {
                if (x_k_->getEntryAt(i) < x_l_->getEntryAt(i)) infea_var += (x_l_->getEntryAt(i) - x_k_->getEntryAt(i));
                else if (x_k_->getEntryAt(i) > x_u_->getEntryAt(i))
                    infea_var += (x_k_->getEntryAt(i) - x_u_->getEntryAt(i));
            }
            for (int i = 0; i < c_k_->Dim(); i++) {
                if (c_k_->getEntryAt(i) < c_l_->getEntryAt(i)) infea_con += (c_l_->getEntryAt(i) - c_k_->getEntryAt(i));
                else if (c_trial_->getEntryAt(i) > c_u_->getEntryAt(i))
                    infea_con += (c_k_->getEntryAt(i) - c_u_->getEntryAt(i));
            }
            infea_measure_ = infea_con + infea_var;
        }

        return true;
    }


    /**
     * This function extracts the search direction for NLP from the QP subproblem solved and copies it to the class member _p_k
     *
     * It will truncate the optimal solution of QP into two parts, the first half (with length equal to the number of variables)
     *to be the search direction.
     *
     * @param qphandler the QPhandler class object used for solving a QP subproblem with specified QP informations
     */
    bool Algorithm::get_search_direction(shared_ptr<SQPhotstart::QPhandler> qphandler) {
        double *tmp_p_k = new double[qphandler->A_->ColNum()]();
        qphandler->GetOptimalSolution(tmp_p_k);
        p_k_->copy_vector(tmp_p_k);
        delete[] tmp_p_k;
        return true;
    }

    /**
     * This function extracts the the Lagragian multipliers for constraints in NLP and copies it to the class member _lambda
     *
     *   Note that the QP subproblem will return a multiplier for the constraints and the bound in a single vector, so we only take
     *   the first #constraints number of elements as an approximation of multipliers for the nlp problem
     *
     * @param qphandler the QPsolver class object used for solving a QP subproblem with specified QP informations
     */
    bool Algorithm::get_multipliers(shared_ptr<QPhandler> qphandler) {
        double *tmp_lambda = new double[qphandler->A_->RowNum() + qphandler->A_->ColNum()]();
        qphandler->GetMultipliers(tmp_lambda);
        lambda_->copy_vector(tmp_lambda);
        delete[] tmp_lambda;
        return true;
    }

    /**
    * This function will set up the data for the QP subproblem if reset = true
     *
     * @param reset
     * 		if reset = true; the data in the QP subproblem will be reset.
     *
     */

    bool Algorithm::setupQP() {

        if (stats->iter == 0) {
            myQP->setup_bounds(delta_, x_k_, c_k_, x_l_, x_u_, c_l_, c_u_);
            myQP->setup_g(grad_f_, rho_);
            myQP->setup_H(hessian_);
            myQP->setup_A(jacobian_);
        } else {
            if (QPinfoFlag_.Update_A) {
                myQP->update_A(jacobian_);
                QPinfoFlag_.Update_A = false;
            }
            if (QPinfoFlag_.Update_H) {
                myQP->update_H(hessian_);
                QPinfoFlag_.Update_H = false;
            }
            if (QPinfoFlag_.Update_bounds) {
                myQP->update_bounds(delta_);
                QPinfoFlag_.Update_bounds = false;
            }

            if (QPinfoFlag_.Update_penalty) {
                myQP->update_penalty(rho_);
                QPinfoFlag_.Update_penalty = false;
            }
            if (QPinfoFlag_.Update_grad) {
                myQP->update_grad(grad_f_);
                QPinfoFlag_.Update_grad = false;
            }
        }
        return true;
    }

    /**
     *
     * This function performs the ratio test to determine if we should accept the trial point
     *
     * The ratio is calculated by (P_1(x_k;\rho)-P_1(x_trial;\rho))/(q_k(0;\rho)-q_k(p_k;rho),
     *
     * where
     * 	P_1(x,rho) = f(x) + rho* infeasibility_measure
     * is the l_1 merit function and
     * 	q_k(p; rho) = f_k+ g_k^Tp +1/2 p^T H_k p +rho* infeasibility_measure_model
     * is the quadratic model at x_k.
     *
     * Tne trial point  will be accepted if the ratio >= eta_s. If it is accepted, the function will also updates the gradient, Jacobian
     * 	information by reading from _nlp object. _update will be set to true, so the SOC direction will not be calculated. _reset_qp is
     * 	set to true, meaning that the data in the QP subproblem need to be reset.
     * Otherwise, there is no change except that _update will be set to false, the SOC direction will be calculated.
     *
     */
    bool Algorithm::ratio_test() {
        Number P1x = obj_value_ + rho_ * infea_measure_;
        Number P1_x_trial = obj_value_trial_ + rho_ * infea_measure_trial_;

        actual_reduction_ = P1x - P1_x_trial;
        pred_reduction_ = rho_ * infea_measure_ - myQP->qp_obj_;

        if (actual_reduction_ >= options->eta_s * pred_reduction_) {
            //succesfully update
            //copy information already calculated from the trial point
            infea_measure_ = infea_measure_trial_;
            obj_value_ = obj_value_trial_;
            x_k_->copy_vector(x_trial_->vector());
            c_k_->copy_vector(c_trial_->vector());
            //update function information by reading from nlp_ object
            get_multipliers(myQP);
            nlp_->Eval_gradient(x_k_, grad_f_);
            nlp_->Eval_Jacobian(x_k_, jacobian_);
            nlp_->Eval_Hessian(x_k_, lambda_, hessian_);
            QPinfoFlag_.Update_A = true;
            QPinfoFlag_.Update_H = true;
            QPinfoFlag_.Update_grad = true;
            isaccept_ = true;    //no need to calculate the SOC direction
        } else {
            isaccept_ = false;
        }
        return true;
    }

    /**
     *
     * This function update the trust-region radius when the ratio calculated by the ratio test is smaller than eta_c or bigger than
     * 	eta_e and the search_direction hits the trust-region bounds.
     *
     * If ratio<eta_c, the trust region radius will decrease by the parameter gamma_c, to be gamma_c*_delta
     * If ratio_test> eta_e and _delta = _norm_p_k, the trust-region radius will be increased by the parameter gamma_c.
     *
     * In either of these two cases, and if the trial point does not pass the ratio_test (the x_trial hasn't been accepted), only the
     * trust region parameter of the QP bounds will be update, and the _reset_qp will be set to false.
     *
     */

    bool Algorithm::radius_update() {
        if (actual_reduction_ < options->eta_c * pred_reduction_) {
            delta_ = options->gamma_c * delta_;
            //decrease the trust region radius. gamma_c is the parameter in options object
        } else {
            if (actual_reduction_ > options->
                    eta_e * pred_reduction_
                && options->tol > (delta_ - norm_p_k_)) {
                delta_ = std::min(options->gamma_e * delta_, options->delta_max);
            }
        }
        return true;
    }

/**
 *
 * This function checks how the constraints specified by the nlp readers are bounded
 * TODO: fix comments
 * If there is only upper bounds for constraints, c(x)<=c_u, then _Constraint_type = BOUNDED_ABOVE
 * If there is only lower bounds for constraints, c(x)>=c_l, then _Constraint_type = BOUNDED_BELOW
 * If there are both upper bounds and lower bounds, c_l<=c(x)<=c_u, then _Constraint_type = BOUNDED,
 * If there is no constraints on all of c_i(x), then _Constraint_type = UNBOUNDED;
 *
 * It is similar for variables:
 * If there is only upper bounds for variables, x<=x_u, then _Variable_type = BOUNDED_ABOVE
 * If there is only lower bounds for variables, x>=x_l, then _Variable_type = BOUNDED_BELOW
 * If there are both upper bounds and lower bounds, x_l<=x<=x_u, then _Variable_type = BOUNDED,
 * If there is no bounds on all of c_i(x), then _Variable_type = UNBOUNDED.
 *
 */
    bool Algorithm::ClassifyConstraintType() {
        for (int i = 0; i < nCon_; i++) {
            cons_type_[i] = classify_single_constraint(x_l_->getEntryAt(i), x_u_->getEntryAt(i));
        }
        for (int i = 0; i < nVar_; i++) {
            bound_cons_type_[i] = classify_single_constraint(x_l_->getEntryAt(i), x_u_->getEntryAt(i));
        }
        return true;
    }

}//END_NAMESPACE_SQPHOTSTART
