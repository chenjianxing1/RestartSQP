#ifndef SQPHOTSTART_ALG_HPP_
#define SQPHOTSTART_ALG_HPP_

#include "IpTNLP.hpp"
#include "Stats.hpp"
#include "Types.hpp"
#include "Options.hpp"
#include "QPhandler.hpp"
#include "Utils.hpp"
#include "Log.hpp"
#include "SQPTNLP.hpp"
#include "Vector.hpp"
#include "Matrix.hpp"


namespace SQPhotstart {
    /**
     *
     * This is the class with method solve a NLP problem by using SQP(SL1QP)
     *
     *It can solve a problem in the following format
     *
     *	minimize 	f(x)
     *	subject   c_l<=c(x)<=c_u,
     *		  x_l<= x  <=u
     *
     *
     *To use this method, call @Optimize and input NLP class object.
     *
     */
    class Algorithm {
    public:
        /** Default Constructor*/
        Algorithm();

        /** Default Destructor*/
        virtual ~Algorithm();

        /**
         * This is the main method to optimize the NLP given as the input
         *
         * @param nlp: the nlp reader that read data of the function to be minimized;
         */
        virtual bool Optimize(SmartPtr <Ipopt::TNLP> nlp);




        /* Private methods*/
    private:

        /** Copy Constructor */
        Algorithm(const Algorithm &);

        /** Overloaded Equals Operator */
        void operator=(const Algorithm &);

        /**
             *
             *This is the function that checks if the current point is optimal, and decides if to exit the loop or not
             *
             *@return
             * 	 if it decides the function is optimal, the class member _exitflag = OPTIMAL
             * 	 if it decides that there is an error during the function run or the function cannot be solved, it will assign _exitflag the
             * 	 	corresponding code according to the error type.
             */
        virtual bool termination_check();

        /**
         * This function initializes the objects required by the SQP Algorithm, copies some parameters required by the algorithm,
         * obtains the function information for the first QP, and solve the first QP and LP.
         *
         */
        virtual bool Presolve();

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
         *			infea_fun = norm(-max(c-cu,0),1)+norm(-min(c-cl,0),1);
         *  			infea_var = norm(-min(x-bl,0),1)+norm(-max(x-bu,0),1);
         *  			infea_measure = infea_fun+infea_var;
         *
         */

        virtual bool infea_cal(bool trial);

        virtual bool second_order_correction() { return false; };

        /** Perform the ratio test to determine if we should accept the trial point*/
        virtual bool ratio_test();

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
        virtual bool radius_update();

        /**
         * Update the penalty parameter
         */
        virtual bool penalty_update() {
            return false;
        };

        /**
         * This function extracts the search direction for NLP from the QP subproblem solved and copies it to the class member _p_k
         *
         * It will truncate the optimal solution of QP into two parts, the first half (with length equal to the number of variables)
         *to be the search direction.
         *
         * @param qphandler the qphandler class object used for solving a QP subproblem with specified QP informations
         */
        bool get_search_direction(shared_ptr<QPhandler> qphandler);


        /**
         * This function will set up the data for the QP subproblem if reset = true
         *
         * @param reset
         * 		if reset = true; the data in the QP subproblem will be reset.
         *
         */

        bool setupQP();

        /**
         * This function extracts the the Lagragian multipliers for constraints in NLP and copies it to the class
         * member _lambda
         *
         *   Note that the QP subproblem will return a multiplier for the constraints and the bound in a single vector,
         *   so we only take the first #constraints number of elements as an approximation of multipliers for the nlp
         *   problem
         *
         * @param qphandler the QPsolver class object used for solving a QP subproblem with specified QP informations
         */

        bool get_multipliers(shared_ptr<QPhandler> qphandler);

        /**
         * This function initializes all the shared pointer which will be used in the Algorithm::Optimize, and it copies all parameters
         * that might be changed during the run of the function Algorithm::Optimize
         *
         * @param nlp: the nlp reader that read data of the function to be minimized;
         */
        bool allocate(SmartPtr <Ipopt::TNLP> nlp);

        /** Check how the constraints are bounded
         *  If there is only upper bounds for constraints, c(x)<=c_u, then _Constraint_type = BOUNDED_ABOVE
         *  If there is only lower bounds for constraints, c(x)>=c_l, then _Constraint_type = BOUNDED_BELOW
         *  If there are both upper bounds and lower bounds, c_l<=c(x)<=c_u, then _Constraint_type = BOUNDED,
         *  If there is no constraints on all of c_i(x), then _Constraint_type = UNBOUNDED;
         *
         *  It is similar for variables:
         *  If there is only upper bounds for variables, x<=x_u, then _Variable_type = BOUNDED_ABOVE
         *  If there is only lower bounds for variables, x>=x_l, then _Variable_type = BOUNDED_BELOW
         *  If there are both upper bounds and lower bounds, x_l<=x<=x_u, then _Variable_type = BOUNDED,
         *  If there is no bounds on all of c_i(x), then _Variable_type = UNBOUNDED;
         */
        bool ClassifyConstraintType();

        /* public class members */
    private:
        Index nVar_; /* number of variables*/
        Index nCon_; /* number of constraints*/
        shared_ptr<SQPTNLP> nlp_;

        shared_ptr<Vector> lambda_;/*multiplier for constraints evaluated at x_k*/
        shared_ptr<Vector> grad_f_;/*gradient evaluated at x_k*/

        Number infea_measure_;/* the measure of infeasibility evaluated at x_k*/
        Number obj_value_;/*the objective corresponding to the x_k*/
        shared_ptr<Vector> x_k_; /* current iterate point*/
        shared_ptr<Vector> c_k_; /* the constraints' value evaluated at x_k_*/

        Number infea_measure_trial_;/* the measure of infeasibility evaluated at x_trial*/
        Number obj_value_trial_;/*the objective corresponding to the x_trial*/
        shared_ptr<Vector> x_trial_;/* the trial point from the search direction, x_trial = x_k+p_k*/
        shared_ptr<Vector> c_trial_;/* the constraints' value evaluated at x_trial_*/

        shared_ptr<Vector> c_l_; /* the lower bounds for constraints*/
        shared_ptr<Vector> c_u_; /* the upper constraints vector*/
        shared_ptr<Vector> x_l_; /* the lower bounds for variables*/
        shared_ptr<Vector> x_u_; /* the upper bounds for variables*/
        shared_ptr<Vector> p_k_; /* search direction at x_k*/
        shared_ptr<Matrix> hessian_;/* the Matrix object for hessain of f(x)+sum_{i=1}^m lambda_i c_i(x) //TODO: check the sign*/
        shared_ptr<Matrix> jacobian_;/*the Matrix object for Jacobian from c(x)*/
        ConstraintType* cons_type_; /* the constraints type, it can be either bounded, bounded above,bounded below, or unbounded*/
        ConstraintType* bound_cons_type_;/* the variables type, it can be either bounded, bounded above,bounded below, or unbounded*/
        shared_ptr<Options> options;/* the default options used for now. TODO: modify it*/
        shared_ptr<Stats> stats;
        shared_ptr<QPhandler> myQP;
        shared_ptr<QPhandler> myLP;
        shared_ptr<Log> log;

        Number norm_p_k_;/*the infinity norm of p_k*/
        Number delta_;/*trust-region radius*/
        Number rho_; /*penalty parameter*/
        Number pred_reduction_;/* the predicted reduction evaluated at x_k and p_k*/
        Number actual_reduction_; /* the actual_reduction evaluated at x_k and p_k*/
        Exitflag exitflag_ = UNKNOWN;
        bool isaccept_; // is the new point accepted?
        UpdateFlags QPinfoFlag_; /*indicates which QP problem bounds should be updated*/
    };//END_OF_ALG_CLASS



}//END_NAMESPACE_SQPHOTSTART




#endif

