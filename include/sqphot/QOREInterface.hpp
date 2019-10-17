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


    const shared_ptr<Vector>& getG() const override {
        return g_;
    };

    const shared_ptr<Vector>& getLb() const override {
        return lb_;
    };

    const shared_ptr<Vector>& getUb() const override {
        return ub_;
    };

    const shared_ptr<Vector>& getLbA() const override {
        THROW_EXCEPTION(INVALID_RETURN_TYPE,INVALID_RETURN_TYPE_MSG);
    }

    const shared_ptr<Vector>& getUbA() const override {
        THROW_EXCEPTION(INVALID_RETURN_TYPE,INVALID_RETURN_TYPE_MSG);
    }

    shared_ptr<const SpHbMat> getH()const override {
        return H_;
    };

    shared_ptr<const SpHbMat> getA() const override {
        return A_;
    };


    /**Constructor*/

    QOREInterface(NLPInfo nlp_info,
                  QPType qptype,
                  shared_ptr<const Options> options,
                  Ipopt::SmartPtr<Ipopt::Journalist> jnlst);


    /** Default destructor*/
    ~QOREInterface();


    /**
     * @brief Solve a regular QP with given data and options.
     */

    void optimizeQP(shared_ptr<Stats> stats) override;


    /**
     * @brief Solve a regular LP with given data and options
     */

    void optimizeLP(shared_ptr<Stats> stats) override;

    /**@name Getters */
    //@{
    /**
     *@brief get the objective value from the QP solvers
     *
     * @return the objective function value of the QP problem
     */
    double get_obj_value() override;

    /**
     * @brief copy the optimal solution of the QP to the input pointer
     *
     * @param x_optimal a pointer to an empty array with allocated memory euqal to
     * sizeof(double)*number_variables
     */
    inline double* get_optimal_solution() override {
        return x_qp_->values();
    };

    /**
     * @brief get the pointer to the multipliers to the bounds constraints.
     */
    inline double* get_multipliers_bounds()override {
        return y_qp_->values();
    };

    /**
     * @brief get the pointer to the multipliers to the regular constraints.
     */
    inline double* get_multipliers_constr()override {
        return y_qp_->values()+ nVar_QP_;
    };

    Exitflag get_status() override;

    void get_working_set(ActiveType* W_constr, ActiveType* W_bounds) override;

    //@}
    //

    /** @name Setters */
    //@{
    void set_g(int location, double value) override {
        value = value < INF ? value : INF;
        g_->setValueAt(location, value);
    };

    void set_lb(int location, double value) override {
        value = value > -INF ? value : -INF;
        lb_->setValueAt(location, value);
    };

    void set_ub(int location, double value) override {
        value = value < INF ? value : INF;
        ub_->setValueAt(location, value);
    };

    void set_A_structure(shared_ptr<const SpTripletMat> rhs,
                         Identity2Info I_info) override {
        A_->setStructure(rhs, I_info);
    };

    void set_A_values(shared_ptr<const SpTripletMat> rhs, Identity2Info
                      I_info) override {
        A_->setMatVal(rhs, I_info);
    };

    void set_H_structure(shared_ptr<const SpTripletMat> rhs) override {
        H_->setStructure(rhs);
    };

    void set_H_values(shared_ptr<const SpTripletMat> rhs) override {
        H_->setMatVal(rhs);
    };

    //@}
    void WriteQPDataToFile(Ipopt::EJournalLevel level,
                           Ipopt::EJournalCategory category,
                           const string filename) override ;



    /** Just overload from base class, does not use them here though..*/
    //@{
    void set_g(shared_ptr<const Vector> rhs) override {};
    void set_lb(shared_ptr<const Vector> rhs) override {};
    void set_ub(shared_ptr<const Vector> rhs) override {};
    void set_lbA(int location, double value) override {};
    void set_lbA(shared_ptr<const Vector> rhs) override {};
    void set_ubA(int location, double value) override {};
    void set_ubA(shared_ptr<const Vector> rhs) override {};
    void reset_constraints() override {};
    //@}

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

    /**
     * @brief set options of QP solver based on the user-defined values
     */
    void set_solver_options(shared_ptr<const Options> options);

    /**
     * @brief Handle errors based on current status
     */
    void handle_error(QPType qptype);
    /**
     * @brief Allocate memory for the class members
     * @param nlp_index_info  the struct that stores simple nlp dimension info
     * @param qptype is the problem to be solved QP or LP?
     */
    void allocate_memory(NLPInfo nlp_info, QPType qptype);


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

    int qpiter_[1];
    Ipopt::SmartPtr<Ipopt::Journalist> jnlst_;
    int rv_;//temporarily placed here, for recording the return value from the solver
    int* working_set_;
    //TODO: for debugging use only
    shared_ptr<SpTripletMat> H_triplet_;
    shared_ptr<SpTripletMat> A_triplet_;
};


}


#endif
