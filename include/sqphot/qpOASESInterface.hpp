/* Copyright (C) 2019
 * All Rights Reserved.
 *
 * Authors: Xinyi Luo
 * Date:    2019-07
 */

#ifndef __QPOASES_INTERFACE_HPP__
#define __QPOASES_INTERFACE_HPP__

#include <sqphot/QPsolverInterface.hpp>


namespace SQPhotstart {
/**
 * @brief This is a derived class of QPsolverInterface.
 * It uses qpOASES as the QP solver  which features the hotstart option. It is used
 * as the default QP solver for SQPhostart.
 */
enum QPMatrixType {
    UNDEFINED,
    FIXED,
    VARIED
};


class qpOASESInterface :
    public QPSolverInterface {

    ///////////////////////////////////////////////////////////
    //                      PUBLIC METHODS                   //
    ///////////////////////////////////////////////////////////
public:

    /**
     * @brief Constructor which also initializes the qpOASES SQProblem objects
     * @param nlp_index_info the struct that stores simple nlp dimension info
     * @param qptype  is the problem to be solved QP or LP?
     */
    qpOASESInterface(Index_info nlp_index_info, QPType qptype,
                     shared_ptr<const Options> options);    //number of constraints in the QP problem

    /** Defualt Destructor */
    ~qpOASESInterface() override;


#if DEBUG
#if GET_QP_INTERFACE_MEMBERS or COMPARE_QP_SOLVER
    const shared_ptr<Vector>& getLb() const override {
        return lb_;
    };

    const shared_ptr<Vector>& getUb() const override {
        return ub_;
    };


    const shared_ptr<Vector>& getLbA() const override {
        return lbA_;
    };

    const shared_ptr<Vector>& getUbA() const override {
        return ubA_;
    };

    const shared_ptr<Vector>& getG() const override {
        return g_;
    };


    const shared_ptr<const SpTripletMat> getH() const override {
        H_triplet_->convert2Triplet(H_);
        return H_triplet_;
    };

    const shared_ptr<const SpTripletMat> getA() const override {
        A_triplet_->convert2Triplet(A_);
        return A_triplet_;
    };

#endif
#endif




    void optimizeQP(shared_ptr<Stats> stats) override;


    /**
     * @brief optimize the LP problem whose objective and constraints are defined
     * in the class members.
     */

    void optimizeLP(shared_ptr<Stats> stats) override;


    /**
    * @brief copy the optimal solution of the QP to the input pointer
    *
    * @param p_k a pointer to an empty array with allocated memory equals to
    * sizeof(double)*number_variables
    *
    */

    double* get_optimal_solution() override;


    /**
     * @brief copy the multipliers of the QP to the input pointer
     *
     * @param y_k   a pointer to an array with allocated memory equals to
     * sizeof(double)*(num_variable+num_constraint)
     */
    double* get_multipliers() override;


    /**
     *@brief get the objective value from the QP solvers
     *
     * @return the objective function value of the QP problem
     */

    double get_obj_value() override;


    /**
     * @brief get the final return status of the QP problem
     */

    QPReturnType get_status() override;

    /** @name Setters */
    //@{
    void set_lb(int location, double value) override;


    void set_lb(shared_ptr<const Vector> rhs) override;


    void set_ub(int location, double value) override;


    void set_ub(shared_ptr<const Vector> rhs) override;


    void set_lbA(int location, double value) override;


    void set_lbA(shared_ptr<const Vector> rhs) override;


    void set_ubA(int location, double value) override;


    void set_ubA(shared_ptr<const Vector> rhs) override;


    void set_g(int location, double value) override;


    void set_g(shared_ptr<const Vector> rhs) override;


    void set_H_structure(shared_ptr<const SpTripletMat> rhs) override;


    void set_H_values(shared_ptr<const SpTripletMat> rhs) override;


    void set_A_structure(shared_ptr<const SpTripletMat> rhs, Identity2Info I_info)
    override;


    void
    set_A_values(shared_ptr<const SpTripletMat> rhs, Identity2Info I_info) override;


    void GetWorkingSet(ActiveType* W_constr, ActiveType* W_bounds) override ;


    void WriteQPDataToFile(
        Ipopt::SmartPtr<Ipopt::Journalist> jnlst,
        Ipopt::EJournalLevel level,
        Ipopt::EJournalCategory category) override ;
    //@}


    ///////////////////////////////////////////////////////////
    //                      PRIVATE METHODS                  //
    ///////////////////////////////////////////////////////////
private:

    /** default constructor*/
    qpOASESInterface();



    void handler_error(QPType qptype, shared_ptr<Stats> stats);


    /**
     * @brief get the final return status of the QP problem
     */

    void reset_flags();


    /**
     * @brief obtain an exit status from QP solver and change the class member status_
     */
    void obtain_status();


    /**
     * @brief Allocate memory for the class members
     * @param nlp_index_info  the struct that stores simple nlp dimension info
     * @param qptype is the problem to be solved QP or LP?
     */
    void allocate(Index_info nlp_index_info, QPType qptype);


    void setQP_options();



    void get_Matrix_change_status();


    /** Copy Constructor */
    qpOASESInterface(const qpOASESInterface &);


    /** Overloaded Equals Operator */
    void operator=(const qpOASESInterface &);



    ///////////////////////////////////////////////////////////
    //                      PRIVATE MEMBERS                  //
    ///////////////////////////////////////////////////////////
private:


    Ipopt::SmartPtr<Ipopt::Journalist> jnlst_;
    QPMatrixType new_QP_matrix_status_ = UNDEFINED;
    QPMatrixType old_QP_matrix_status_ = UNDEFINED;
    QPReturnType status_;
    UpdateFlags data_change_flags_;
    bool firstQPsolved_ = false; /**< if the first QP has been solved? */
    int nConstr_QP_;  /**< number of constraints for QP*/
    int nVar_QP_;  /**< number of variables for QP*/
    shared_ptr<const Options> options_;
    shared_ptr<qpOASES::SymSparseMat> H_qpOASES_;/**< the Matrix object that qpOASES
                                                       * taken in, it only contains the
                                                       * pointers to array stored in
                                                       * the class members of H_*/

    shared_ptr<qpOASES::SQProblem> solver_;/**< the qpOASES object used for solving a qp*/
    shared_ptr<qpOASES::SparseMatrix> A_qpOASES_;/**< the Matrix object that qpOASES
                                                       * taken in, it only contains the
                                                       * pointers to array stored in
                                                       * the class members of A_*/
    shared_ptr<Vector> g_;   /**< the grad used for QPsubproblem*/
    shared_ptr<Vector> lbA_; /**< lower bounds of Ax */
    shared_ptr<Vector> lb_;  /**< lower bounds of x */
    shared_ptr<Vector> ubA_; /**< upper bounds of Ax */
    shared_ptr<Vector> ub_;  /**< upper bounds of x */
    shared_ptr<Vector> x_qp_; /** the qp solution */
    shared_ptr<Vector> y_qp_; /** the multiplier corresponding to the optimal solution */
    shared_ptr<SpHbMat> H_;/**< the Matrix object stores the QP data H in
                                          * Harwell-Boeing Sparse Matrix format*/
    shared_ptr<SpHbMat> A_;/**< the Matrix object stores the QP data A in
                                          * Harwell-Boeing Sparse Matrix format*/
#if DEBUG
#if GET_QP_INTERFACE_MEMBERS or COMPARE_QP_SOLVER
    shared_ptr<SpTripletMat> H_triplet_;
    shared_ptr<SpTripletMat> A_triplet_;
#endif
#endif

};
}
#endif
