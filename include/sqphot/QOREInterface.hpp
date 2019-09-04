/* Copyright (C) 2019
 * All Rights Reserved.
 *
 * Authors: Xinyi Luo
 * Date:    2019-08
 */
#ifndef __SQPHOTSTART_QOREINTERFACE_HPP__
#define __SQPHOTSTART_QOREINTERFACE_HPP__

extern "C" {
#include "qpsolver.h"
}

#include <sqphot/QPsolverInterface.hpp>



DECLARE_STD_EXCEPTION(INVALID_RETURN_TYPE);
namespace SQPhotstart {
class QOREInterface :
    public QPSolverInterface {

    ///////////////////////////////////////////////////////////
    //                      PUBLIC METHODS                   //
    ///////////////////////////////////////////////////////////
public:

#if DEBUG
#if GET_QP_INTERFACE_MEMBERS or COMPARE_QP_SOLVER

    const shared_ptr<Vector>& getG() const override {
        return g_;
    };

    const shared_ptr<Vector>& getLb() const override {
        return lb_;
    };

    const shared_ptr<Vector>& getLbA() const override {
        THROW_EXCEPTION(INVALID_RETURN_TYPE,INVALID_RETURN_TYPE_MSG);
    }

    const shared_ptr<Vector>& getUb() const override {
        return ub_;
    };

    const shared_ptr<Vector>& getUbA() const override {
        THROW_EXCEPTION(INVALID_RETURN_TYPE,INVALID_RETURN_TYPE_MSG);
    }


    const shared_ptr<const SpTripletMat> getH()const override {
        H_triplet_->convert2Triplet(H_);
        return H_triplet_;
    };

    const shared_ptr<const SpTripletMat> getA() const override {
        A_triplet_->convert2Triplet(A_);
        return A_triplet_;
    };

#endif
#endif

    /**Constructor*/

    QOREInterface(Index_info nlp_info,
                  QPType qptype,
                  shared_ptr<const Options> options,
                  Ipopt::SmartPtr<Ipopt::Journalist> jnlst);


    /** Default destructor*/
    ~QOREInterface();


    /**
     * @brief Solve a regular QP with given data and options.
     *
     * overload this method to optimize a QP with the data specified, update the
     * stats by adding the iteration number used to solve this QP to stats.qp_iter
     */

    void optimizeQP(shared_ptr<Stats> stats) override;


    /**
     * @brief Solve a regular LP with given data and options
     * overload this method to optimize a LP with the data specified, update the
     * stats by adding the iteration number used to solve this QP to stats.qp_iter
     */

    void optimizeLP(shared_ptr<Stats> stats) override;

    /**@name Getters */
    //@{

    /**
     * @brief copy the optimal solution of the QP to the input pointer
     *
     * @param x_optimal a pointer to an empty array with allocated memory euqal to
     * sizeof(double)*number_variables
     */
    double* get_optimal_solution() override;


    /**
     *@brief get the objective value from the QP solvers
     *
     * @return the objective function value of the QP problem
     */
    double get_obj_value() override;


    /**
     * @brief copy the multipliers of the QP to the input pointer
     *
     * @param y_k   a pointer to an array with allocated memory
     */
    double* get_multipliers() override;

    QPReturnType get_status() override;;

    //@}
    //

    /** @name Setters */
    //@{
    void set_H_structure(shared_ptr<const SpTripletMat> rhs) override;;


    void set_H_values(shared_ptr<const SpTripletMat> rhs) override;;


    void set_g(int location, double value) override;;


    void set_g(shared_ptr<const Vector> rhs) override {};


    void set_lb(int location, double value) override;;


    void set_lb(shared_ptr<const Vector> rhs) override {
    };


    void set_ub(int location, double value) override;;

    void set_ub(shared_ptr<const Vector> rhs) override {};

    void set_A_structure(shared_ptr<const SpTripletMat> rhs,
                         Identity2Info I_info) override;;


    void set_A_values(shared_ptr<const SpTripletMat> rhs, Identity2Info
                      I_info) override;;




    /** Just overload from base class, does not use them here though..*/
    //@{
    void set_lbA(int location, double value) override {};


    void set_lbA(shared_ptr<const Vector> rhs) override {};


    void set_ubA(int location, double value) override {};


    void set_ubA(shared_ptr<const Vector> rhs) override {};
    //@}

    //@}
    void WriteQPDataToFile(Ipopt::SmartPtr<Ipopt::Journalist> jnlst,
                           Ipopt::EJournalLevel level,
                           Ipopt::EJournalCategory category) override ;

    void GetWorkingSet(ActiveType* W_constr, ActiveType* W_bounds) override ;
    ///////////////////////////////////////////////////////////
    //                      PRIVATE METHODS                  //
    ///////////////////////////////////////////////////////////
private:
    /** Default Constructor*/
    QOREInterface();

    /** Copy Constructor */
    QOREInterface(const QOREInterface &);


    /** Overloaded Equals Operator */
    void operator=(const QOREInterface &);


    void setQP_options(shared_ptr<Options> options);

    /**
     * @brief Allocate memory for the class members
     * @param nlp_index_info  the struct that stores simple nlp dimension info
     * @param qptype is the problem to be solved QP or LP?
     */
    void allocate_memory(Index_info nlp_info, QPType qptype);




    ///////////////////////////////////////////////////////////
    //                      PRIVATE MEMBERS                  //
    ///////////////////////////////////////////////////////////
private:
    int status_;
    QoreProblem* solver_;
    bool firstQPsolved_ = false;
    int nConstr_QP_;
    int nVar_QP_;
    shared_ptr<SpHbMat> A_;
    shared_ptr<SpHbMat> H_;
    shared_ptr<Vector> g_;
    shared_ptr<Vector> lb_;
    shared_ptr<Vector> ub_;
    shared_ptr<Vector> x_qp_;
    shared_ptr<Vector> y_qp_;
    Ipopt::SmartPtr<Ipopt::Journalist> jnlst_;
    int rv_;//temporarily placed here, for recording the return value from the solver
    int* working_set_;
#if DEBUG
#if GET_QP_INTERFACE_MEMBERS or COMPARE_QP_SOLVER
    shared_ptr<SpTripletMat> H_triplet_;
    shared_ptr<SpTripletMat> A_triplet_;
#endif
#endif
};


}


#endif
