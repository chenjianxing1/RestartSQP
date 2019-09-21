/* Copyright (C) 2019
 * All Rights Reserved.
 *
 * Authors: Xinyi Luo
 * Data: 2019-09-09
 */
#ifndef SQPHOTSTART_CPLEX_INTERFACE_HPP
#define SQPHOTSTART_CPLEX_INTERFACE_HPP


#ifdef USE_CPLEX
#include <ilcplex/ilocplex.h>
#endif

#include <sqphot/QPsolverInterface.hpp>

namespace SQPhotstart {
DECLARE_STD_EXCEPTION(CPLEX_SOLVER_FAILS);

class CplexInterface : public QPSolverInterface {
    /**-------------------------------------------------------**/
    /**                  Public Methods                       **/
    /**-------------------------------------------------------**/
public:
    /**@name Getters for private members*/
    //@{
    const shared_ptr<Vector>& getLb() const override;

    const shared_ptr<Vector>& getUb() const override;

    const shared_ptr<Vector>& getLbA() const override;

    const shared_ptr<Vector>& getUbA() const override;

    const shared_ptr<Vector>& getG() const override;

    const shared_ptr<const SpTripletMat> getH() const override;

    const shared_ptr<const SpTripletMat> getA() const override;
    //@}


    CplexInterface(Index_info nlp_info,
                   QPType qptype,
                   shared_ptr<const Options> options,
                   Ipopt::SmartPtr<Ipopt::Journalist> jnlst);

    /** Default destructor*/
    ~CplexInterface() override;

    /**
     * @brief Solve a regular QP with given data and options.
     *
     */
    void
    optimizeQP(shared_ptr<Stats> stats) override;

    /**
     * @brief Solve a regular LP with given data and options
     */

    void
    optimizeLP(shared_ptr<Stats> stats) override;


    /**-------------------------------------------------------**/
    /**                    Getters                            **/
    /**-------------------------------------------------------**/
    /**@name Getters*/
    //@{
    /**
     * @return the pointer to the optimal solution
     *
     */
    double* get_optimal_solution() override;

    /**
     *@brief get the objective value from the QP solvers
     *
     * @return the objective function value of the QP problem
     */
    double get_obj_value() override;


    /**
     * @brief get the pointer to the multipliers to the bounds constraints.
     */
    double* get_multipliers_bounds() override;

    /**
     * @brief get the pointer to the multipliers to the regular constraints.
     */
    double* get_multipliers_constr() override;

    /**
     * @brief copy the working set information
     * @param W_constr a pointer to an array of length (nCon_QP_) which will store the
     * working set for constraints
     * @param W_bounds a pointer to an array of length (nVar_QP_) which will store the
     * working set for bounds
     *
     */
    void get_working_set(ActiveType* W_constr, ActiveType* W_bounds)override;

    Exitflag get_status() override;

    //@}

    /**-------------------------------------------------------**/
    /**                    Setters                            **/
    /**-------------------------------------------------------**/
    /**@name Setters, by location and value*/
    //@{
    void set_lb(int location, double value) override;

    void set_ub(int location, double value) override;

    void set_lbA(int location, double value) override;

    void set_ubA(int location, double value) override;

    void set_g(int location, double value) override;
    //@}

    /**@name Setters for dense vector, by vector value*/
    //@{
    void set_ub(shared_ptr<const Vector> rhs) override;

    void set_lb(shared_ptr<const Vector> rhs) override;

    void set_lbA(shared_ptr<const Vector> rhs) override;

    void set_ubA(shared_ptr<const Vector> rhs) override;

    void set_g(shared_ptr<const Vector> rhs) override;
    //@}

    /**@name Setters for matrix*/
    //@{
    void set_H_structure(shared_ptr<const SpTripletMat> rhs) override;

    void set_H_values(shared_ptr<const SpTripletMat> rhs) override;

    void set_A_structure(shared_ptr<const SpTripletMat> rhs, Identity2Info
                         I_info) override;

    void set_A_values(shared_ptr<const SpTripletMat> rhs, Identity2Info
                      I_info) override;
    //@}

    void reset_constraints() override;

    /**-------------------------------------------------------**/
    /**                  Data Writer                          **/
    /**-------------------------------------------------------**/

    void WriteQPDataToFile(Ipopt::SmartPtr<Ipopt::Journalist> jnlst,
                           Ipopt::EJournalLevel level,
                           Ipopt::EJournalCategory category) override;

    /**-------------------------------------------------------**/
    /**                  Private Methods                      **/
    /**-------------------------------------------------------**/
private:
    /** Default constructor*/
    CplexInterface();

    /** Copy Constructor */
    CplexInterface(const CplexInterface &);

    /** Overloaded Equals Operator */
    void operator=(const CplexInterface &);

    void set_solver_options();
    /**-------------------------------------------------------**/
    /**                  Private Members                      **/
    /**-------------------------------------------------------**/

private:
    Identity2Info I_info_;
    Ipopt::SmartPtr<Ipopt::Journalist> jnlst_;
    Exitflag status_;
    QPType qptype_;
    bool firstQPsolved_;
    double final_obj_;
    int nConstr_QP_;
    int nVar_QP_;
    shared_ptr<Vector> lb_;
    shared_ptr<Vector> ub_;
    shared_ptr<Vector> x_qp_;
    shared_ptr<Vector> y_qp_;
    shared_ptr<const Options> options_;
    shared_ptr<const SpTripletMat> A_;
#ifdef USE_CPLEX
    shared_ptr<IloEnv> cplex_env_;
    shared_ptr<IloModel> cplex_model_;
    IloNumExpr qobj_;/**< quadratic part of the objecitve*/
    vector<IloNumExpr> constraints_;
    vector<IloNumExpr> lterm_;/**< linear part of the objective */
    vector<IloNumVar> cplex_vars_;
#endif


};


}//SQPHOTSTART
#endif


