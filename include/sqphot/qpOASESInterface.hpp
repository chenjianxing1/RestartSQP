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
     * @param nlp_info the struct that stores simple nlp dimension info
     * @param qptype  is the problem to be solved QP or LP?
     */
    qpOASESInterface(NLPInfo nlp_info, QPType qptype,
                     shared_ptr<const Options> options,
                     Ipopt::SmartPtr<Ipopt::Journalist> jnlst);    //number of constraints in the QP problem


    qpOASESInterface(shared_ptr<SpHbMat> H,
                     shared_ptr<SpHbMat> A,
                     shared_ptr<Vector> g,
                     shared_ptr<Vector> lb,
                     shared_ptr<Vector> ub,
                     shared_ptr<Vector> lbA,
                     shared_ptr<Vector> ubA,
                     shared_ptr<Options> options = nullptr);

    /** Defualt Destructor */
    ~qpOASESInterface() override;




    void optimizeQP(shared_ptr<Stats> stats = nullptr) override;



    /**
     * @brief optimize the LP problem whose objective and constraints are defined
     * in the class members.
     */

    void optimizeLP(shared_ptr<Stats> stats = nullptr) override;

    /** @name Getters*/
//@{
    /**
    * @brief copy the optimal solution of the QP to the input pointer
    */

    double* get_optimal_solution() override;

    /**
     * @brief get the pointer to the multipliers to the bounds constraints.
     */
    double* get_multipliers_bounds()override;

    /**
     * @brief get the pointer to the multipliers to the regular constraints.
     */
    double* get_multipliers_constr()override;


    void get_working_set(ActiveType* W_constr, ActiveType* W_bounds) override;

    /**
     *@brief get the objective value from the QP solvers
     *
     * @return the objective function value of the QP problem
     */

    double get_obj_value() override;


    /**
     * @brief get the final return status of the QP problem
     */

    Exitflag get_status() override;
//@}

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


    void set_H(shared_ptr<const SpTripletMat> rhs) override;

    void set_A(shared_ptr<const SpTripletMat> rhs, IdentityInfo I_info) override;

    //@}

    void WriteQPDataToFile(Ipopt::EJournalLevel level,
                           Ipopt::EJournalCategory category,
                           const string filename) override ;


    void reset_constraints() override;

    //@{
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


    shared_ptr<const SpHbMat> getH() const override {
        return H_;
    };

    shared_ptr<const SpHbMat> getA() const override {
        return A_;
    };
    //@}

    bool test_optimality(ActiveType* W_c = NULL, ActiveType* W_b = NULL) override ;


    ///////////////////////////////////////////////////////////
    //                      PRIVATE METHODS                  //
    ///////////////////////////////////////////////////////////
private:

    /** default constructor*/
    qpOASESInterface();



    void handle_error(QPType qptype, shared_ptr<Stats> stats = nullptr);


    void reset_flags();


    /**
     * @brief Allocate memory for the class members
     * @param nlp_info  the struct that stores simple nlp dimension info
     * @param qptype is the problem to be solved QP or LP?
     */
    void allocate_memory(NLPInfo nlp_info, QPType qptype);


    void set_solver_options();



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
    UpdateFlags data_change_flags_;
    bool firstQPsolved_ = false; /**< if the first QP has been solved? */
    int nConstr_QP_;  /**< number of constraints for QP*/
    int nVar_QP_;  /**< number of variables for QP*/
//    OptimalityStatus qpOptimalStatus_;
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

};
}
#endif

