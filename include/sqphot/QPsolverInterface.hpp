/* Copyright (C) 2019
 * All Rights Reserved.
 *
 * Authors: Xinyi Luo      2019-07
 */
#ifndef SQPHOTSTART_QPSOLVER_INTERFACE_HPP
#define SQPHOTSTART_QPSOLVER_INTERFACE_HPP

#include <memory>
#include <vector>
#include <IpException.hpp>
#include <sqphot/Vector.hpp>
#include <sqphot/Utils.hpp>
#include <sqphot/Matrix.hpp>
#include <qpOASES.hpp>
#include <sqphot/Stats.hpp>
#include <sqphot/Options.hpp>
#include <sqphot/Types.hpp>
//#include <sqphot/SQPDebug.hpp>


using namespace std;
namespace SQPhotstart {

DECLARE_STD_EXCEPTION(QP_NOT_OPTIMAL);

DECLARE_STD_EXCEPTION(LP_NOT_OPTIMAL);

DECLARE_STD_EXCEPTION(QP_INTERNAL_ERROR);

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
#if DEBUG
#if GET_QPOASES_MEMBERS
    virtual const shared_ptr<Vector>& getLb() const = 0;

    virtual const shared_ptr<Vector>& getUb() const = 0;

    virtual const shared_ptr<Vector>& getLbA() const =0;

    virtual const shared_ptr<Vector>& getUbA() const =0;

    virtual const shared_ptr<Vector>& getG() const = 0;



#endif
#endif
    QPSolverInterface() {}

    /** Default destructor*/
    virtual ~QPSolverInterface() {}

    /**
     * @brief Solve a regular QP with given data and options.
     *
     * overload this method to optimize a QP with the data specified, update the
     * stats by adding the iteration number used to solve this QP to stats.qp_iter
     */
    virtual void
    optimizeQP(shared_ptr<Stats> stats, shared_ptr<Options> options) = 0;

    /**
     * @brief Solve a regular LP with given data and options
     * overload this method to optimize a LP with the data specified, update the
     * stats by adding the iteration number used to solve this QP to stats.qp_iter
     */

    virtual void
    optimizeLP(shared_ptr<Stats> stats, shared_ptr<Options> options) = 0;

    /**
     * @brief copy the optimal solution of the QP to the input pointer
     *
     * @param x_optimal a pointer to an empty array with allocated memory euqal to
     * sizeof(double)*number_variables
     */
    virtual double* get_optimal_solution() = 0;

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
    virtual double* get_multipliers() = 0;


    /**
     * Return private class members info
     */

    virtual void set_lb(int location, double value) = 0;

    virtual void set_lb(shared_ptr<const Vector> rhs) = 0;

    virtual void set_ub(int location, double value) = 0;

    virtual void set_ub(shared_ptr<const Vector> rhs) = 0;

    virtual void set_lbA(int location, double value) = 0;

    virtual void set_lbA(shared_ptr<const Vector> rhs) = 0;

    virtual void set_ubA(int location, double value) = 0;

    virtual void set_ubA(shared_ptr<const Vector> rhs) = 0;

    virtual void set_g(int location, double value) = 0;

    virtual void set_g(shared_ptr<const Vector> rhs) = 0;

    virtual void set_H_structure(shared_ptr<const SpTripletMat> rhs) = 0;

    virtual void set_H_values(shared_ptr<const SpTripletMat> rhs) = 0;

    virtual void set_A_structure(shared_ptr<const SpTripletMat> rhs, Identity2Info
                                 I_info) = 0;

    virtual void set_A_values(shared_ptr<const SpTripletMat> rhs, Identity2Info
                              I_info) = 0;

    virtual QPReturnType get_status() = 0;

    //@}

private:
    /** Copy Constructor */
    QPSolverInterface(const QPSolverInterface &);

    /** Overloaded Equals Operator */
    void operator=(const QPSolverInterface &);
};


}//SQPHOTSTART
#endif


