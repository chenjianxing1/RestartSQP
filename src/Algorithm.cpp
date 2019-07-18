/* Copyright (C) 2019
 * All Rights Reserved.
 *
 * Authors: Xinyi Luo
 * Date:2019-06
 */
#include <sqphot/Algorithm.hpp>
#include <sqphot/Utils.hpp>

namespace SQPhotstart {

    /**
     * Default Constructor
     */
    Algorithm::Algorithm() {
        setDefaultOption();
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
     * @brief This is the main function to optimize the NLP given as the input
     *
     * @param nlp: the nlp reader that read data of the function to be minimized;
     */
    bool Algorithm::Optimize(SmartPtr<Ipopt::TNLP> nlp) {
        allocate(nlp); //allocate memory to class members
        initilization();

        /** Main iteration */

        while (stats->iter < options->iter_max) {
            //setup the QPs problems
//            std::cout<<QPinfoFlag_.Update_A;
//            std::cout<<QPinfoFlag_.Update_bounds;
//            std::cout<<QPinfoFlag_.Update_H;
//            std::cout<<QPinfoFlag_.Update_grad;

            setupQP();
            //based on the information given by NLP reader otherwise, it will do nothing
            myQP->solveQP(stats, options);//solve the QP subproblem and update the stats

            //get the search direction from the solution of the QPsubproblem
            get_search_direction(myQP);

            //calculate the infinity norm of the search direction
            norm_p_k_ = p_k_->getInfNorm();

            //Update the penalty parameter if necessary
            penalty_update();

            // Calculate the the trial points, x_trial = x_k+p_k
            x_trial_->copy_vector(x_k_->values());
            x_trial_->add_vector(p_k_->values());

            // Calculate f_trial, c_trial and infea_measure_trial for the trial points
            // x_trial
            nlp_->Eval_f(x_trial_, obj_value_trial_);
            nlp_->Eval_constraints(x_trial_, c_trial_);

            infea_cal(true);

//                        std::cout<<"infea_measure is"<<infea_measure_<<std::endl;
//                        std::cout<<"infea_measure_trial is"<<infea_measure_trial_<<std::endl;

            ratio_test();

            // Calculate the second-order-correction steps
            second_order_correction();

            // Update the radius and the QP bounds if the radius has been changed

            stats->iter_addone();

            /* output some information to the console, TODO: change it*/
            if (stats->iter % 10 == 0)log->print_header();
            if (options->printLevel > 0)
                log->print_main_iter(stats->iter, obj_value_, norm_p_k_, infea_measure_,
                                     delta_, rho_);

            radius_update();

            //check if the current iterates is optimal and decide to
            //exit the loop or not
            termination_check();

            if (exitflag_ != UNKNOWN) {
                break;
            }
        }
        // print the final summary message to the console
        log->print_final(stats->iter, stats->qp_iter, obj_value_, norm_p_k_,
                         infea_measure_, exitflag_);
        return true;
    }

    /**
     * @brief This is the function that checks if the current point is optimal, and
     * decides if to exit the loop or not
     * *@return if it decides the function is optimal, the class member _exitflag =
     * OPTIMAL
     * if it decides that there is an error during the function run or the
     *  function cannot be solved, it will assign _exitflag the	corresponding
     *  code according to the error type.
     */
    bool Algorithm::termination_check() {
        if (norm_p_k_ < options->tol) {
            exitflag_ = OPTIMAL;
        }
        return true;
    }


    /**
     *  @brief This function initializes the objects required by the SQP Algorithm,
     *  copies some parameters required by the algorithm, obtains the function
     *  information for the first QP.
     *
     */
    bool Algorithm::initilization() {
        nlp_->Get_bounds_info(x_l_, x_u_, c_l_, c_u_);
        nlp_->Get_starting_point(x_k_, lambda_);
        //shift starting point to satisfy the bound constraint
        nlp_->shift_starting_point(x_k_, x_l_, x_u_);
        nlp_->Eval_f(x_k_, obj_value_);
        nlp_->Eval_gradient(x_k_, grad_f_);
        nlp_->Eval_constraints(x_k_, c_k_);
        nlp_->Get_Structure_Hessian(x_k_, lambda_, hessian_);
        nlp_->Eval_Hessian(x_k_, lambda_, hessian_);
        nlp_->Get_Strucutre_Jacobian(x_k_, jacobian_);
        nlp_->Eval_Jacobian(x_k_, jacobian_);
        ClassifyConstraintType();

        infea_cal(false); //calculate the infeasibility measure for x_k
        // initializes QP objects*/
        myQP->init(nlp_->nlp_info_, QP);
        myLP->init(nlp_->nlp_info_, LP);
        log->print_header();
        log->print_main_iter(stats->iter, obj_value_, 0.0, infea_measure_, delta_, rho_);
        return true;
    }

    /**
     * @brief alloocate memory for class members.
     * This function initializes all the shared pointer which will be used in the 
     * Algorithm::Optimize, and it copies all parameters that might be changed during 
     * the run of the function Algorithm::Optimize.
     *
     * @param nlp: the nlp reader that read data of the function to be minimized;
     */
    bool Algorithm::allocate(SmartPtr<Ipopt::TNLP> nlp) {
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
        jacobian_ = make_shared<SpTripletMat>(nlp_->nlp_info_.nnz_jac_g, nCon_, nVar_);
        hessian_ = make_shared<SpTripletMat>(nlp_->nlp_info_.nnz_h_lag, nVar_, nVar_);

        options = make_shared<Options>();
        stats = make_shared<Stats>();
        log = make_shared<Log>();
        myQP = make_shared<QPhandler>();
        myLP = make_shared<LPhandler>();

        delta_ = options->delta;
        rho_ = options->rho;


        return true;
    }


    /**
     * @brief This function calculates the infeasibility measure for either trial point
     * or current iterate x_k
     *
     *@param trial: true if the user are going to evaluate the infeasibility measure of
     * the trial point _x_trial;
     *	infea_measure_trial = norm(-max(c_trial-cu,0),1)+norm(-min(c_trial-cl,0),1)
     *
     * 	            false if the user are going to evaluate the infeasibility measure
     * 	            of the current iterates _x_k
     *	infea_measure = norm(-max(c-cu,0),1)+norm(-min(c-cl,0),1);
     *
     */
    bool Algorithm::infea_cal(bool trial) {

        if (trial) {
            infea_measure_trial_ = 0;
            for (int i = 0; i < c_k_->Dim(); i++) {
                if (c_trial_->values()[i] < c_l_->values()[i])
                    infea_measure_trial_ += (c_l_->values()[i] - c_trial_->values()[i]);
                else if (c_trial_->values()[i] > c_u_->values()[i])
                    infea_measure_trial_ += (c_trial_->values()[i] - c_u_->values()[i]);
            }
        } else {
            infea_measure_ = 0;
            for (int i = 0; i < c_k_->Dim(); i++) {
                if (c_k_->values()[i] < c_l_->values()[i])
                    infea_measure_ += (c_l_->values()[i] - c_k_->values()[i]);
                else if (c_trial_->values()[i] > c_u_->values()[i])
                    infea_measure_ += (c_k_->values()[i] - c_u_->values()[i]);
            }
        }
        return true;
    }


    /**
     * @brief This function extracts the search direction for NLP from the QP subproblem
     * solved, and copies it to the class member _p_k
     *
     * It will truncate the optimal solution of QP into two parts, the first half (with
     * length equal to the number of variables) to be the search direction.
     *
     * @param qphandler the QPhandler class object used for solving a QP subproblem with
     * specified QP information
     */
    bool Algorithm::get_search_direction(shared_ptr<SQPhotstart::QPhandler> qphandler) {

        double *tmp_p_k = new double[nVar_ + 2 * nCon_];
        qphandler->GetOptimalSolution(tmp_p_k);
        p_k_->copy_vector(tmp_p_k);
        if (options->penalty_update)
            infea_measure_model_ = oneNorm(tmp_p_k + nVar_, 2 * nCon_);
        //FIXME:calculate somewhere else?
        delete[] tmp_p_k;
        return true;
    }


    /**
     * @brief This function extracts the the Lagragian multipliers for constraints
     * in NLP and copies it to the class member lambda_
     *
     *   Note that the QP subproblem will return a multiplier for the constraints
     *   and the bound in a single vector, so we only take the first #constraints
     *   number of elements as an approximation of multipliers for the nlp problem
     *
     * @param qphandler the QPsolver class object used for solving a QP subproblem
     * with specified QP information
     */

    bool Algorithm::get_multipliers(shared_ptr<QPhandler> qphandler) {
        double *tmp_lambda = new double[nVar_ + 3 * nCon_];
        qphandler->GetMultipliers(tmp_lambda);
        lambda_->copy_vector(tmp_lambda);

        //        lambda_->print();
        delete[] tmp_lambda;
        return true;
    }

    /**
     * @brief This function will set up the data for the QP subproblem
     *
     * It will initialize all the data at once at the beginning of the Algorithm. After
     * that, the data in the QP problem will be updated according to the class
     * member QPinfoFlag_
     */

    bool Algorithm::setupQP() {

        if (stats->iter == 0) {
            myQP->setup_bounds(delta_, x_k_, x_l_, x_u_, c_k_, c_l_, c_u_);
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
                myQP->update_bounds(delta_, x_k_, x_l_, x_u_, c_k_, c_l_, c_u_);
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
     * @brief This function performs the ratio test to determine if we should accept
     * the trial point
     *
     * The ratio is calculated by
     * (P_1(x_k;\rho)-P_1( x_trial;\rho))/(q_k(0;\rho)-q_k(p_k;rho), where
     * P_1(x,rho) = f(x) + rho* infeasibility_measure is the l_1 merit function and
     * q_k(p; rho) = f_k+ g_k^Tp +1/2 p^T H_k p+rho* infeasibility_measure_model is the
     * quadratic model at x_k.
     * The trial point  will be accepted if the ratio >= eta_s.
     * If it is accepted, the function will also updates the gradient, Jacobian
     * information by reading from nlp_ object. The corresponding flags of class member
     * QPinfoFlag_ will set to be true.
     */
    bool Algorithm::ratio_test() {
        using namespace std;
        Number P1x = obj_value_ + rho_ * infea_measure_;
        Number P1_x_trial = obj_value_trial_ + rho_ * infea_measure_trial_;
        actual_reduction_ = P1x - P1_x_trial;
//                cout<<"actual reduction is"<<actual_reduction_<<endl;

        pred_reduction_ = rho_ * infea_measure_ - myQP->get_obj();

//                cout<<"pred reduction is"<<pred_reduction_<<endl;
        if (actual_reduction_ >= options->eta_s * pred_reduction_) {
            //succesfully update
            //copy information already calculated from the trial point
            infea_measure_ = infea_measure_trial_;
            obj_value_ = obj_value_trial_;
            x_k_->copy_vector(x_trial_->values());
            c_k_->copy_vector(c_trial_->values());
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
     * @brief Update the trust region radius.
     * 
     * This function update the trust-region radius when the ratio calculated by the 
     * ratio test is smaller than eta_c or bigger than eta_e and the search_direction 
     * hits the trust-region bounds. 
     * If ratio<eta_c, the trust region radius will decrease by the parameter 
     * gamma_c, to be gamma_c* delta_
     * If ratio_test> eta_e and delta_= norm_p_k_ the trust-region radius will be 
     * increased by the parameter gamma_c.
     * 
     * If trust region radius has changed, the corresponding flags will be set to be
     * true;
     */

    bool Algorithm::radius_update() {
        if (actual_reduction_ < options->eta_c * pred_reduction_) {
            delta_ = options->gamma_c * delta_;
            QPinfoFlag_.Update_bounds = true;
            //decrease the trust region radius. gamma_c is the parameter in options object
        } else {
            if (actual_reduction_ > options->
                    eta_e * pred_reduction_
                && options->tol > (delta_ - norm_p_k_)) {
                delta_ = std::min(options->gamma_e * delta_, options->delta_max);
                QPinfoFlag_.Update_bounds = true;
            }
        }
        return true;
    }

    /**
     *
     * @brief This function checks how each constraint specified by the nlp readers are
     * bounded.
     * If there is only upper bounds for a constraint, c_i(x)<=c^i_u, then 
     * cons_type_[i]= BOUNDED_ABOVE
     * If there is only lower bounds for a constraint, c_i(x)>=c^i_l, then 
     * cons_type_[i]= BOUNDED_BELOW
     * If there are both upper bounds and lower bounds, c^i_l<=c_i(x)<=c^i_u, and 
     * c^i_l<c^i_u then cons_type_[i]= BOUNDED,
     * If there is no constraints on all 
     * of c_i(x), then cons_type_[i]= UNBOUNDED;
     * 
     * The same rules are also applied to the bound-constraints.
     */


    bool Algorithm::ClassifyConstraintType() {
        for (int i = 0; i < nCon_; i++) {
            cons_type_[i] = classify_single_constraint(x_l_->values()[i],
                                                       x_u_->values()[i]);
        }
        for (int i = 0; i < nVar_; i++) {
            bound_cons_type_[i] = classify_single_constraint(x_l_->values()[i],
                                                             x_u_->values()[i]);
        }
        return true;
    }

    /**
     * @brief update the penalty parameter for the algorithm.
     *
     */
    bool Algorithm::penalty_update() {
        if (options->penalty_update) {
            if (infea_measure_model_ > options->penalty_update_tol) {
                myLP->copy_QP_info(myQP);
                myLP->solveLP(stats, options);
                shared_ptr<Vector> sol_tmp = make_shared<Vector>(nVar_ + 2 * nCon_);

                double rho_trial = rho_;
                //calculate the infea_measure of the LP
                get_full_search_direction(myLP, sol_tmp);
                double infea_measure_infty = oneNorm(sol_tmp->values() + nVar_,
                                                     2 * nCon_);
                if (infea_measure_infty <= options->penalty_update_tol) {
                    //try to increase the penalty parameter to a number such that the
                    // infeasibility measure of QP model with such penalty parameter
                    // becomes zero
                    while (infea_measure_model_ > options->penalty_update_tol) {
                        rho_trial = options->increase_parm * rho_trial; //increase rho
                        stats->penalty_change_trial_addone();
                        myQP->update_penalty(rho_trial);
                        myQP->solveQP(stats, options);
                        get_full_search_direction(myQP, sol_tmp);
                        infea_measure_model_ = oneNorm(sol_tmp->values() + nVar_,
                                                       2 * nCon_);
                    }

                } else {

                    while ((infea_measure_ - infea_measure_model_ <
                            options->eps1 * (infea_measure_ - infea_measure_infty) &&
                            (stats->penalty_change_trial < options->penalty_iter_max))) {
                        rho_trial = options->increase_parm * rho_trial; //increase rho
                        stats->penalty_change_trial_addone();
                        myQP->update_penalty(rho_trial);
                        myQP->solveQP(stats, options);
                        get_full_search_direction(myQP, sol_tmp);
                        sol_tmp->print();
                        infea_measure_model_ = oneNorm(sol_tmp->values() + nVar_,
                                                       2 * nCon_);

                    }
                }
                //if any change occurs
                if(rho_trial>rho_) {
                    if(rho_trial*infea_measure_-myQP->get_obj()>=
                    options->eps2*rho_trial*(infea_measure_-infea_measure_model_)){
                        stats->penalty_change_Succ_addone();
                        options->eps1+=(1-options->eps1)*0.1;
                        p_k_->copy_vector(sol_tmp);
                        rho_ = rho_trial;
                    }
                    else{
                        stats->penalty_change_Fail_addone();
                    }
                }
            }
        }
        return true;
    }

    bool Algorithm::setDefaultOption() {
        roptions = new Ipopt::RegisteredOptions();
        roptions->SetRegisteringCategory("rust-region");
        roptions->AddNumberOption("eta_c", "trust-region parameter for the ratio test.",
                                  0.25,
                                  "If ratio<=eta_c, then the trust-region radius for the next "
                                  "iteration will be decreased for the next iteration.");
        roptions->AddNumberOption("eta_s", "trust-region parameter for the ratio test.",
                                  1.0e-8,
                                  "The trial point will be accepted if ratio>= eta_s. ");
        roptions->AddNumberOption("eta_e", "trust-region parameter for the ratio test.",
                                  0.75,
                                  "If ratio>=eta_e and the search direction hits the  "
                                  "trust-region boundary, the trust-region radius will "
                                  " be increased for the next iteration.");
        roptions->AddNumberOption("gamma_c", "radius update parameter",
                                  0.5,
                                  "If the trust-region radius is going to be decreased,"
                                  " then it will be set as gamma_c*delta, where delta "
                                  "is current trust-region radius.");
        roptions->AddNumberOption("gamma_e", "radius update parameter",
                                  2.0,
                                  "If the trust-region radius is going to be "
                                  "increased, then it will be set as gamma_e*delta,"
                                  "where delta is current trust-region radius.");
        roptions->AddNumberOption("delta_0", "initial trust-region radius value", 1.0);
        roptions->AddNumberOption("delta_max", "the maximum value of trust-region radius"
                                               " allowed for the radius update", 1.0e8);

        roptions->SetRegisteringCategory("Penalty Update");
        roptions->AddNumberOption("eps1", "penalty update parameter", 0.3, "");
        roptions->AddNumberOption("eps2", "penalty update parameter", 1.0e-6, "");
        roptions->AddNumberOption("print_level_penalty_update", "print level for penalty "
                                                                "update",
                                  0);

        roptions->SetRegisteringCategory("Optimality Test");

        roptions->AddNumberOption("opt_tol", "", 1.0e-5);
        roptions->AddNumberOption("opt_compl_tol", "", 1.0e-6);
        roptions->AddNumberOption("opt_dual_fea_tol", " ", 1.0e-6);
        roptions->AddNumberOption("opt_prim_fea_tol", " ", 1.0e-5);
        roptions->AddNumberOption("opt_second_tol", " ", 1.0e-8);

        roptions->SetRegisteringCategory("General");
        roptions->AddNumberOption("iter_max", "maximum number of iteration for the "
                                              "algorithm", 10);
        roptions->AddNumberOption("print_level", "print level for main algorithm", 2);
        roptions->AddStringOption2(
                "second_order_correction",
                "Tells the algorithm to calculate the second-order correction step "
                "during the main iteration"
                "yes",
                "no", "not calculate the soc steps",
                "yes", "will calculate the soc steps",
                "");

        roptions->SetRegisteringCategory("QPsolver");
        roptions->AddNumberOption("iter_max_qp", "maximum number of iteration for the "
                                                 "QP solver in solving each QP", 10);
        roptions->AddNumberOption("print_level_qp", "print level for QP solver", 0);

        return true;
    }

    bool Algorithm::get_full_search_direction(
            shared_ptr<SQPhotstart::QPhandler> qphandler,
            shared_ptr<SQPhotstart::Vector> search_direction) {
        qphandler->GetOptimalSolution(search_direction->values());
        return true;
    }

    /**
     * @brief
     * @param lphandler
     * @param search_direction
     * @return
     */
    bool Algorithm::get_full_search_direction(shared_ptr<LPhandler> lphandler,
                                              shared_ptr<Vector> search_direction) {

        lphandler->GetOptimalSolution(search_direction->values());
        return false;
    }

}//END_NAMESPACE_SQPHOTSTART
