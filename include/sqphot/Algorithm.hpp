/* Copyright (C) 2019
 * All Rights Reserved.
 *
 * Authors: Xinyi Luo
 * Date:    2019-07
 */
#ifndef SQPHOTSTART_ALG_HPP_
#define SQPHOTSTART_ALG_HPP_

#include <sqphot/SQPDebug.hpp>
#include <coin/IpTNLP.hpp>
#include <coin/IpRegOptions.hpp>
#include <coin/IpOptionsList.hpp>
#include <sqphot/Stats.hpp>
#include <sqphot/Types.hpp>
#include <sqphot/Options.hpp>
#include <sqphot/QPhandler.hpp>
#include <sqphot/LPhandler.hpp>
#include <sqphot/Utils.hpp>
#include <sqphot/Log.hpp>
#include <sqphot/SQPTNLP.hpp>
#include <sqphot/Vector.hpp>
#include <sqphot/Matrix.hpp>


namespace SQPhotstart {
/**
 *
 * @brief This is the class with method solve a NLP problem by using SQP(SL1QP)
 *
 *It can solve a problem in the following format
 *
 *	minimize 	f(x)
 *	subject   c_l<=c(x)<=c_u,
 *	      	  x_l<= x  <=x_u
 *
 *
 *To use this method, call @Optimize and input NLP class object.
 *
 */
class Algorithm {
public:
    /** @name constructor/destructor*/
    //@{
    /** Default Constructor*/
    Algorithm();

    /** Default Destructor*/
    ~Algorithm();
    //@}

    /** temporarily use Ipopt options*/
    //@{
    //

    SmartPtr<RegisteredOptions> getRoptions() {
        return roptions;
    }

    SmartPtr<OptionsList> getRoptions2() {
        return roptions2_;
    }

    SmartPtr<Journalist> getJnlst() {
        return jnlst_;
    }

    //@}


    /**
     * @brief This is the main method to optimize the NLP given as the input
     *
     * @param nlp: the nlp reader that read data of the function to be minimized;
     */
    virtual void Optimize(SmartPtr<Ipopt::TNLP> nlp);

    /** @name Set the corresponding option to the user-defined value */
    //@{
    template<typename T>
    bool setOptions(const std::string& name, T value) {
        return false;
    };
//@}
//
    /** @name Getters*/
    //@{
    inline Exitflag getExitFlag() {
        return exitflag_;
    }

    int* getActiveSetConstraints() const {
        return NULL;
    };

    int* getActiveSetBounds() const {
        return NULL;
    };

    //@}
    /* Private methods*/
private:

    /** Copy Constructor */
    Algorithm(const Algorithm&);

    /** Overloaded Equals Operator */
    void operator=(const Algorithm&);

    /**
     * @brief set the default option values
     */
    void setDefaultOption();

    /**
     * @brief This is the function that checks if the current point is optimal, and
     * decides if to exit the loop or not
     *
     * If it decides the function is optimal, the class member _exitflag =
     * OPTIMAL
     * if it decides that there is an error during the function run or the
     *  function cannot be solved, it will assign _exitflag the	corresponding
     *  code according to the error type.
     */
    void termination_check();

    /**
     *  @brief This function initializes the objects required by the SQP Algorithm,
     *  copies some parameters required by the algorithm, obtains the function
     *  information for the first QP.
     *
     */
    void initialization(SmartPtr<Ipopt::TNLP> nlp);


    /**
     * @brief This function calculates the infeasibility measure for  current
     * iterate x_k
     *
     *	infea_measure = norm(-max(c-cu,0),1)+norm(-min(c-cl,0),1);
     *
     */
    void cal_infea();

    /**
     * @brief This function calculates the infeasibility measure for  current
     * iterate x_k
     *
     *   infea_measure_trial = norm(-max(c_trial-cu,0),1)+norm(-min(c_trial-cl,0),1)
     *
     *
     */
    void cal_infea_trial();

    /**
     * @brief Calculate the trial point based on current search direction,
     * x_trial = x_k+p_k,
     *       and get the funcion value, constraints value, and infeasibility measure
     *       at the trial point
     */
    void get_trial_point_info();;

    /**
     * @brief calculate the second order correction step and decide if accept the
     * calculated step or not.
     * It will be calculated only if the second_order_correction in options is set
     * to be true.
     *
     */

    void second_order_correction();

    /**
     *
     * @brief This function performs the ratio test to determine if we should accept
     * the trial point
     *
     * The ratio is calculated by
     * (P_1(x_k;\rho)-P_1( x_trial;\rho))/(q_k(0;\rho)-q_k(p_k;rho), where
     * P_1(x,rho) = f(x) + rho* infeasibility_measure is the l_1 merit function and
     * q_k(p; rho) = f_k+ g_k^Tp +1/2 p^T H_k p+rho* infeasibility_measure_model is
     * the quadratic model at x_k.
     * The trial point  will be accepted if the ratio >= eta_s.
     * If it is accepted, the function will also updates the gradient, Jacobian
     * information by reading from nlp_ object. The corresponding flags of class
     * member QPinfoFlag_ will set to be true.
     */
    void ratio_test();

    /**
     * @brief Update the trust region radius.
     *
     * This function update the trust-region radius when the ratio calculated by the
     * ratio test is smaller than eta_c or bigger than eta_e and the search_direction
     * hits the trust-region bounds.
     * If ratio<eta_c, the trust region radius will decrease by the parameter
     * gamma_c, to be gamma_c* delta_
     * If ratio > eta_e and delta_= norm_p_k_ the trust-region radius will be
     * increased by the parameter gamma_c.
     *
     * If trust region radius has changed, the corresponding flags will be set to be
     * true;
     */
    void update_radius();

    /**
     * @brief Update the penalty parameter
     */
    void update_penalty_parameter();


    /**@name Get the search direction from the LP/QP handler*/

    //@{
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

    void get_search_direction();


    /**
     * @brief get the full search direction(including the slack variables) from the
     * QPhandler
     * @param qphandler a QP handler objbect which has already solved a QP
     * @param search_direction a Vector object which will be used to store the full
     * direction
     */
    void get_full_direction_QP(shared_ptr<SQPhotstart::Vector> search_direction);

    /**
     * @brief get the full search direction(including the slack variables) from the
     * LPhandler
     * @param lphandler an LP handler objbect which has already solved a LP
     * @param search_direction a Vector object which will be used to store the full
     * direction
     */
    void get_full_direction_LP(shared_ptr<Vector> search_direction);

    //@}

    /**
    * @brief This function will set up the data for the QP subproblem
    *
    * It will initialize all the data at once at the beginning of the Algorithm. After
    * that, the data in the QP problem will be updated according to the class
    * member QPinfoFlag_
    */


    void setupQP();

    /**
     * @brief This function will setup the data for the LP subproblem
     */
    void setupLP();
    /**
     * @brief This function extracts the Lagragian multipliers for constraints
     * in NLP and copies it to the class member multiplier_cons_.
     *
     *   Note that the QP subproblem will return a multiplier for the constraints
     *   and the bound in a single vector, so we only take the first #constraints
     *   number of elements as an approximation of multipliers for the nlp problem
     *
     * @param qphandler the QPhandler class object used for solving a QP subproblem
     * with specified QP information
     */

    void get_multipliers();

    /**
     * @brief alloocate memory for class members.
     * This function initializes all the shared pointer which will be used in the
     * Algorithm::Optimize, and it copies all parameters that might be changed during
     * the run of the function Algorithm::Optimize.
     *
     * @param nlp: the nlp reader that read data of the function to be minimized;
     */
    void allocate_memory(SmartPtr<Ipopt::TNLP> nlp);

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
    void classify_constraints_types();

    void handle_error(const char* error = NULL);


    /** @name Optimality Test */
//@{
    void IdentifyActiveSet();

    bool Check_KKTConditions(double infea_measure = INF,
                             bool isConstraintChanged = false,
                             bool isPointChanged = false) ;

    /**
     * @brief Test the Second-order optimality conditions
     */
//    bool Check_SecondOrder() ;

    /**
     * @brief Check the Feasibility conditions;
     */

    bool Check_Feasibility(double infea_measure = INF) ;

    /**
     * @brief Check the sign of the multipliers
     */
    bool Check_Dual_Feasibility() ;


    /**
     * @brief Check the complementarity condition
     *
     */

    bool Check_Complementarity() ;

    /**
     * @brief Check the Stationarity condition
     */
    bool Check_Stationarity() ;

//@}
    /* public class members */
private:
    Index nVar_; /**< number of variables*/
    Index nCon_; /**< number of constraints*/
    shared_ptr<SQPTNLP> nlp_;
    shared_ptr<Vector> multiplier_cons_;/**< multiplier for constraints*/
    shared_ptr<Vector> multiplier_vars_;/**< multipliers for variables*/
    shared_ptr<Vector> grad_f_;/**< gradient evaluated at x_k*/
    Number infea_measure_;/**< the measure of infeasibility evaluated at x_k*/
    Number infea_measure_model_; /**< the one norm of all slack variables in the QP */
    Number obj_value_;/**<the objective corresponding to the x_k*/
    shared_ptr<Vector> x_k_; /**< current iterate point*/
    shared_ptr<Vector> c_k_; /**< the constraints' value evaluated at x_k_*/

    Number infea_measure_trial_;/**< the measure of infeasibility evaluated at
                                     * x_trial*/
    Number obj_value_trial_;/**<the objective corresponding to the x_trial*/
    Number qp_obj_;/**< the objective value of current qp*/
    shared_ptr<Vector> x_trial_;/**< the trial point from the search direction,
                                     * x_trial = x_k+p_k*/
    shared_ptr<Vector> c_trial_;/* the constraints' value evaluated at x_trial_*/

    shared_ptr<Vector> c_l_; /* the lower bounds for constraints*/
    shared_ptr<Vector> c_u_; /* the upper constraints vector*/
    shared_ptr<Vector> x_l_; /* the lower bounds for variables*/
    shared_ptr<Vector> x_u_; /* the upper bounds for variables*/
    shared_ptr<Vector> p_k_; /* search direction at x_k*/
    shared_ptr<SpTripletMat> hessian_;/**< the SparseMatrix object for hessain of
                                       * f(x)-sum_{i=1}^m lambda_i c_i(x)*/
    shared_ptr<SpTripletMat> jacobian_;/** <the SparseMatrix object for Jacobian from
                                        * c(x)*/
    ConstraintType* cons_type_; /**<the constraints type, it can be either bounded,
                                     * bounded above,bounded below, or unbounded*/
    ConstraintType* bound_cons_type_;/**< the variables type, it can be either
                                          * bounded, bounded above,bounded below, or
                                          * unbounded*/
    shared_ptr<Options> options;/**< the default options used for now. */
    Ipopt::SmartPtr<Ipopt::RegisteredOptions> roptions; /**FIXME put here for testing.
                                                             it will replace options in
                                                             the future*/
    Ipopt::SmartPtr<Ipopt::OptionsList> roptions2_;
    Ipopt::SmartPtr<Ipopt::Journalist> jnlst_;
    shared_ptr<Stats> stats_;
    shared_ptr<QPhandler> myQP_;
    shared_ptr<LPhandler> myLP_;
    shared_ptr<Log> log_;
    Number norm_p_k_;/**< the infinity norm of p_k*/
    Number delta_;/**< trust-region radius*/
    Number rho_; /**< penalty parameter*/
    Number pred_reduction_;/**< the predicted reduction evaluated at x_k and p_k*/
    Number actual_reduction_; /**< the actual_reduction evaluated at x_k and p_k*/
    int* Active_Set_bounds_;
    int* Active_Set_constraints_;

    Exitflag exitflag_ = UNKNOWN;
    OptimalityStatus opt_status_;
    bool isaccept_; // is the new point accepted?
    UpdateFlags QPinfoFlag_; /**<indicates which QP problem bounds should be updated*/
};//END_OF_ALG_CLASS



}//END_NAMESPACE_SQPHOTSTART




#endif

