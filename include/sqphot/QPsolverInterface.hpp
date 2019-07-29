/* Copyright (C) 2019
 * All Rights Reserved.
 *
 * Authors: Xinyi Luo      2019-07
 */
#ifndef SQPHOTSTART_QPSOLVER_INTERFACE_HPP
#define SQPHOTSTART_QPSOLVER_INTERFACE_HPP

#include <memory>
#include <vector>
#include <coin/IpException.hpp>
#include <sqphot/Vector.hpp>
#include <sqphot/Matrix.hpp>
#include <qpOASES.hpp>
#include <sqphot/Stats.hpp>
#include <sqphot/Options.hpp>
#include <sqphot/Types.hpp>
//#include <sqphot/SQPDebug.hpp>


using namespace std;
namespace SQPhotstart {

DECLARE_STD_EXCEPTION(QP_NOT_OPTIMAL);
/**
 * @brief Base class for all standard QP solvers that use standard triplet matrix
 * form and dense vectors.
 *
 * It can optimize QP problem in the following format
 *
 *  minimize 1/2 x^T H x + g^T x
 *  subject  lb_A<=Ax<=ub_A
 *              lb<=x<=ub
 */
class QPSolverInterface {

public:
    /** Default constructor*/
    QPSolverInterface() {}


    /** Default destructor*/
    virtual ~QPSolverInterface() {}

    /**
     * @brief Solve a regular QP with given data and options.
     *
     * overload this method to optimize a QP with the data specified, update the
     * stats by adding the iteration number used to solve this QP to stats.qp_iter
     */
    virtual bool
    optimizeQP(shared_ptr<Stats> stats, shared_ptr<Options> options) {
        return false;
    }

    /**
     * @brief Solve a regular LP with given data and options
     * overload this method to optimize a LP with the data specified, update the
     * stats by adding the iteration number used to solve this QP to stats.qp_iter
     */

    virtual bool
    optimizeLP(shared_ptr<Stats> stats, shared_ptr<Options> options) {
        return false;
    }

    /**
     * @brief copy the optimal solution of the QP to the input pointer
     *
     * @param x_optimal a pointer to an empty array with allocated memory euqal to
     * sizeof(double)*number_variables
     */
    virtual bool get_optimal_solution(double* x_optimal) = 0;

    /**
     *@brief get the objective value from the QP solvers
     *
     * @return the objective function value of the QP problem
     */
    virtual double get_obj_value() = 0;


    /**
     * @brief copy the multipliers of the QP to the input pointer
     *
     * @param y_k   a pointer to an array with allocated memory
     */
    virtual bool get_multipliers(double* y_optimal) = 0;


    /**
     * Return private class members info
     */

    //@{
    virtual shared_ptr<Vector>& getLb() = 0;

    virtual shared_ptr<Vector>& getUb() = 0;

    virtual shared_ptr<Vector>& getLbA() = 0;

    virtual shared_ptr<Vector>& getUbA() = 0;

    virtual shared_ptr<Vector>& getG() = 0;
    //@}

private:

    /** Copy Constructor */
    QPSolverInterface(const QPSolverInterface&);

    /** Overloaded Equals Operator */
    void operator=(const QPSolverInterface&);
};


DECLARE_STD_EXCEPTION(QP_INTERNAL_ERROR);

/**
 * @brief This is a derived class of QPsolverInterface.
 * It uses qpOASES as the QP solver  which features the hotstart option. It is used
 * as the default QP solver for SQPhostart.
 */
class qpOASESInterface : public QPSolverInterface {
public:

    virtual ~qpOASESInterface();

    /**
     * @brief Constructor which also initializes the qpOASES SQProblem objects
     * @param nlp_index_info the struct that stores simple nlp dimension info
     * @param qptype  is the problem to be solved QP or LP or SOC?
     */
    qpOASESInterface(Index_info nlp_index_info,
                     QPType qptype);    //number of constraints in the QP problem


    bool optimizeQP(shared_ptr<Stats> stats, shared_ptr<Options> options) override;

    /**
     * @brief optimize the LP problem whose objective and constraints are defined
     * in the class members.
     */

    bool optimizeLP(shared_ptr<Stats> stats, shared_ptr<Options> options) override;

    /**
     * @brief copy the optimal solution of the QP to the input pointer
     *
     * @param x_optimal a pointer to an empty array with allocated memory equals to
     * sizeof(double)*number_variables
     *
     */

    bool get_optimal_solution(double* p_k) override;

    /**
     * @brief copy the multipliers of the QP to the input pointer
     *
     * @param y_k   a pointer to an array with allocated memory equals to
     * sizeof(double)*(num_variable+num_constraint)
     */
    bool get_multipliers(double* y_k) override;

    /**
     *@brief get the objective value from the QP solvers
     *
     * @return the objective function value of the QP problem
     */

    double get_obj_value() override;

    /**
     * @brief get the final return status of the QP problem
     */

    inline QPReturnType get_status();

    //@{
    shared_ptr<Vector>& getLb() override;

    shared_ptr<Vector>& getUb() override;

    shared_ptr<Vector>& getLbA() override;

    shared_ptr<Vector>& getUbA() override;

    shared_ptr<Vector>& getG() override;

    shared_ptr<qpOASESSparseMat>& getH();

    shared_ptr<qpOASESSparseMat>& getA();

    //@}

public:
    shared_ptr<qpOASES::SQProblem> qp_;// the qpOASES object used for solving a qp


private:
    shared_ptr<qpOASES::SymSparseMat> H_qpOASES_;/**< the Matrix object that qpOASES
                                                       * taken in, it only contains the
                                                       * pointers to array stored in
                                                       * the class members of H_*/

    shared_ptr<qpOASES::SparseMatrix> A_qpOASES_;/**< the Matrix object that qpOASES
                                                       * taken in, it only contains the
                                                       * pointers to array stored in
                                                       * the class members of A_*/
    shared_ptr<Vector> lb_;  /**< lower bounds of x */
    shared_ptr<Vector> ub_;  /**< upper bounds of x */
    shared_ptr<Vector> lbA_; /**< lower bounds of Ax */
    shared_ptr<Vector> ubA_; /**< upper bounds of Ax */
    shared_ptr<Vector> g_;   /**< the grad used for QPsubproblem*/
    shared_ptr<qpOASESSparseMat> H_;/**< the Matrix object stores the QP data H in
                                          * Harwell-Boeing Sparse Matrix format*/
    shared_ptr<qpOASESSparseMat> A_;/**< the Matrix object stores the QP data A in
                                          * Harwell-Boeing Sparse Matrix format*/
    bool firstQPsolved = false; /**< if the first QP has been solved? */

private:
    /** default constructor*/
    qpOASESInterface();

    /**
     * @brief Allocate memory for the class members
     * @param nlp_index_info  the struct that stores simple nlp dimension info
     * @param qptype is the problem to be solved QP or LP or SOC?
     */
    bool allocate(Index_info nlp_index_info, QPType qptype);

    /** Copy Constructor */
    qpOASESInterface(const qpOASESInterface&);

    /** Overloaded Equals Operator */
    void operator=(const qpOASESInterface&);

};
}//SQPHOTSTART
#endif


