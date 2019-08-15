/* Copyright (C) 2019
 * All Rights Reserved.
 *
 * Authors: Xinyi Luo
 * Date:    2019-08
 */
#ifndef __SQPHOTSTART_QOREINTERFACE_HPP__
#define __SQPHOTSTART_QOREINTERFACE_HPP__

extern "C"{
#include "qpsolver.h"
}
#include <sqphot/QPsolverInterface.hpp>


namespace SQPhotstart {
class QOREInterface :
    public QPSolverInterface {

    ///////////////////////////////////////////////////////////
    //                      PUBLIC METHODS                   //
    ///////////////////////////////////////////////////////////
public:

#if DEBUG
#if GET_QPOASES_MEMBERS

    const shared_ptr<Vector>& getG() const override {};

    const shared_ptr<Vector>& getLb() const override {};

    const shared_ptr<Vector>& getLbA() const override {};

    const shared_ptr<Vector>& getUb() const override {};

    const shared_ptr<Vector>& getUbA() const override {};


#endif
#endif

    /**Constructor*/
    QOREInterface(Index_info nlp_info, QPType qptype);


    /** Default destructor*/
    ~QOREInterface();


    /**
     * @brief Solve a regular QP with given data and options.
     *
     * overload this method to optimize a QP with the data specified, update the
     * stats by adding the iteration number used to solve this QP to stats.qp_iter
     */

    void optimizeQP(shared_ptr<Stats> stats, shared_ptr<Options> options) override;


    /**
     * @brief Solve a regular LP with given data and options
     * overload this method to optimize a LP with the data specified, update the
     * stats by adding the iteration number used to solve this QP to stats.qp_iter
     */

    void optimizeLP(shared_ptr<Stats> stats, shared_ptr<Options> options) override; 

    /**@name Getters */
    //@{

    /**
     * @brief copy the optimal solution of the QP to the input pointer
     *
     * @param x_optimal a pointer to an empty array with allocated memory euqal to
     * sizeof(double)*number_variables
     */
    void get_optimal_solution(double* x_optimal) override {}


    /**
     *@brief get the objective value from the QP solvers
     *
     * @return the objective function value of the QP problem
     */
    double get_obj_value() override {}


    /**
     * @brief copy the multipliers of the QP to the input pointer
     *
     * @param y_k   a pointer to an array with allocated memory
     */
    void get_multipliers(double* y_optimal) override {}

    QPReturnType get_status() override {};

    //@}
    //
    
    /** @name Setters */
    //@{
    void set_H_structure(shared_ptr<const SpTripletMat> rhs) override {};


    void set_H_values(shared_ptr<const SpTripletMat> rhs) override {};


    void set_g(int location, double value) override {};


    void set_g(shared_ptr<const Vector> rhs) override {};


    void set_lb(int location, double value) override {};


    void set_lb(shared_ptr<const Vector> rhs) override {};


    void set_ub(int location, double value) override {};


    void set_ub(shared_ptr<const Vector> rhs) override {};

    void set_A_structure(shared_ptr<const SpTripletMat> rhs,
                         Identity2Info I_info) override {};


    void set_A_values(shared_ptr<const SpTripletMat> rhs, Identity2Info
                      I_info) override {};




    /** Just overload from base class, does not use them here though..*/
    //@{
    void set_lbA(int location, double value) override {};


    void set_lbA(shared_ptr<const Vector> rhs) override {};




    void set_ubA(int location, double value) override {};


    void set_ubA(shared_ptr<const Vector> rhs) override {};
//@}

    //@}

    ///////////////////////////////////////////////////////////
    //                      PRIVATE METHODS                  //
    ///////////////////////////////////////////////////////////
private:
    /** Default Constructor*/
    QOREInterface();

    /** Copy Constructor */
    QOREInterface(const QOREInterface&);


    /** Overloaded Equals Operator */
    void operator=(const QOREInterface&);
    
    void setQP_options(shared_ptr<Options> options);

    /**
     * @brief Allocate memory for the class members
     * @param nlp_index_info  the struct that stores simple nlp dimension info
     * @param qptype is the problem to be solved QP or LP?
     */
    void allocate(Index_info nlp_info, QPType qptype);



    ///////////////////////////////////////////////////////////
    //                      PRIVATE MEMBERS                  //
    ///////////////////////////////////////////////////////////
private:
    QPReturnType status_;
    QoreProblem* qp_;
    bool firstQPsolved_ = false;
    int nConstr_QP_;
    int nVar_QP_;
    shared_ptr<SpHbMat> A_;
    shared_ptr<SpHbMat> H_;
    shared_ptr<Vector> g_;
    shared_ptr<Vector> lb_;
    shared_ptr<Vector> ub_;
};



}


#endif
