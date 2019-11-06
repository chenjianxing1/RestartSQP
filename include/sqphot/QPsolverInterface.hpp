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
#include <sqphot/SpHbMat.hpp>

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

DECLARE_STD_EXCEPTION(INVALID_WORKING_SET);
/**
 * @brief Virtual base class for all standard QP solvers that use standard triplet matrix
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

    /**@name Getters for private members*/
    virtual const shared_ptr<Vector>& getLb() const = 0;

    virtual const shared_ptr<Vector>& getUb() const = 0;

    virtual const shared_ptr<Vector>& getLbA() const =0;

    virtual const shared_ptr<Vector>& getUbA() const =0;

    virtual const shared_ptr<Vector>& getG() const = 0;

    virtual shared_ptr<const SpHbMat> getH() const =0;

    virtual shared_ptr<const SpHbMat> getA() const =0;
    //@}

    /** Default constructor*/
    QPSolverInterface() = default;

    /** Default destructor*/
    virtual ~QPSolverInterface() = default;

    /**
     * @brief Solve a regular QP with given data and options.
     *
     * overload this method to optimize a QP with the data specified, update the
     * stats by adding the iteration number used to solve this QP to stats.qp_iter
     */
    virtual void
    optimizeQP(shared_ptr<Stats> stats) = 0;

    /**
     * @brief Solve a regular LP with given data and options
     *
     * overload this method to optimize a LP with the data specified, update the
     * stats by adding the iteration number used to solve this LP to stats.qp_iter
     */

    virtual void
    optimizeLP(shared_ptr<Stats> stats) = 0;


    /**-------------------------------------------------------**/
    /**                    Getters                            **/
    /**-------------------------------------------------------**/
    /**@name Getters*/
    //@{
    /**
     * @return the pointer to the optimal solution
     *
     */
    virtual double* get_optimal_solution() = 0;

    /**
     *@brief get the objective value from the QP solvers
     *
     * @return the objective function value of the QP problem
     */
    virtual double get_obj_value() = 0;


    /**
     * @brief get the pointer to the multipliers to the bounds constraints.
     */
    virtual double* get_multipliers_bounds() = 0;

    /**
     * @brief get the pointer to the multipliers to the regular constraints.
     */
    virtual double* get_multipliers_constr() = 0;

    /**
     * @brief copy the working set information
     * @param W_constr a pointer to an array of length (nCon_QP_) which will store the
     * working set for constraints
     * @param W_bounds a pointer to an array of length (nVar_QP_) which will store the
     * working set for bounds
     *
     * overload this method by assign each entry of W_constr and W_bounds to be
     * ACTIVE_ABOVE, ACTIVE_BELOW, INACTIVE, or ACTIVE_BOTH_SIDE based on the working
     * set information from the QP solver
     */
    virtual void get_working_set(ActiveType* W_constr, ActiveType* W_bounds)=0;

    virtual Exitflag get_status() = 0;

    virtual bool test_optimality(ActiveType* W_c = NULL, ActiveType* W_b = NULL) = 0;


    OptimalityStatus get_optimality_status() {
        return qpOptimalStatus_;
    }

    //@}

    /**-------------------------------------------------------**/
    /**                    Setters                            **/
    /**-------------------------------------------------------**/
    /**@name Setters, by location and value*/
    //@{
    virtual void set_lb(int location, double value) = 0;

    virtual void set_ub(int location, double value) = 0;

    virtual void set_lbA(int location, double value) = 0;

    virtual void set_ubA(int location, double value) = 0;

    virtual void set_g(int location, double value) = 0;
    //@}

    /**@name Setters for dense vector, by vector value*/
    //@{
    virtual void set_ub(shared_ptr<const Vector> rhs) = 0;

    virtual void set_lb(shared_ptr<const Vector> rhs) = 0;

    virtual void set_lbA(shared_ptr<const Vector> rhs) = 0;

    virtual void set_ubA(shared_ptr<const Vector> rhs) = 0;

    virtual void set_g(shared_ptr<const Vector> rhs) = 0;
    //@}

    /**@name Setters for matrix*/
    //@{

    virtual void set_H(shared_ptr<const SpTripletMat> rhs) = 0;

    virtual void set_A(shared_ptr<const SpTripletMat> rhs, IdentityInfo I_info) = 0;
    //@}

    virtual void reset_constraints() =0;

    /**-------------------------------------------------------**/
    /**                  Data Writer                          **/
    /**-------------------------------------------------------**/

    virtual void WriteQPDataToFile(Ipopt::EJournalLevel level,
                                   Ipopt::EJournalCategory category,
                                   const string filename) = 0;

protected:
    OptimalityStatus qpOptimalStatus_;

private:
    /** Copy Constructor */
    QPSolverInterface(const QPSolverInterface &);

    /** Overloaded Equals Operator */
    void operator=(const QPSolverInterface &);
};


}//SQPHOTSTART
#endif


