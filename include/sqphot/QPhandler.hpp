
#ifndef SQPHOTSTART_QP_HPP_
#define SQPHOTSTART_QP_HPP_

#include <sqphot/Utils.hpp>
#include <qpOASES.hpp>
#include <sqphot/Stats.hpp>
#include <sqphot/Options.hpp>
#include <sqphot/QPsolverInterface.hpp>
#include <sqphot/SQPTNLP.hpp>

namespace SQPhotstart {
/**
 *
 * This is a class for setting up and solving the SQP
 * subproblems for Algorithm::Optimize
 *
 * It contains the methods to setup the QP problem in
 * the following format or similar ones,
 *
 * 	minimize  1/2 p_k^T H_k p^k + g_k^T p_k+rho*(u+v)
 * 	subject to c_l<=c_k+J_k p+u-v<=c_u,
 * 		  x_l<=x_k+p_k<=x_u,
 *		  -delta<=p_i<=delta,
 *		  u,v>=0
 *
 *
 * and transform them into sparse matrix (triplet) and
 * dense vectors' format which can be taken as input
 * of standard QPhandler.
 *
 * It also contains a method which interfaces to the
 * QP solvers that users choose. The interface will pre
 * -process the data required by individual solver.
 */

class QPhandler {
public:
    /** Default constructor */
    QPhandler();

    /** Default destructor */
    virtual ~QPhandler();


    /**
     * @brief Get the optimal solution from the QPsolverinterface
     *
     * This is only an interface for user to avoid call interface directly.
     * @param p_k 	the pointer to an empty array with the length equal to the size
     * of the QP subproblem
     */
    virtual void GetOptimalSolution(double* p_k);


    virtual /**
         *@brief Get the multipliers from the QPhandler_interface
         *
         * This is only an interface for user to avoid call
         * interface directly.
         *
         * @param y_k 	the pointer to an empty array with
         * the length equal to the size of multipliers of the QP
         * subproblem
         */
    void GetMultipliers(double* y_k);

    /**
     * @brief Get the objective value of the QP
     * @param qp_obj the reference to a double variable which will hold the
     * objective value of the qp problem
     */

    double GetObjective();

    /**
     * @brief This function initializes all objects will be used in this class.
     *
     * @param nlp_info
     * 		it contains the information about number of variables, number of
     * 		constraints, number of elements in the Hessian and that of Jacobian.
     * 		The definition of Index_info is in Types.hpp
     * @param Constraint_type
     * 		it specifies how the variables are bounded. Please check
     * 		@ClassifyConstraintType in Algorithm.hpp
     * 		for more details
     * @param qptype
     * 		it specifies which type of QP is going to be solved. It can be either
     * 		QP, or SOC
     */
    virtual void init(Index_info nlp_info, QPType qptype);

    /**
     *
     * @brief setup the bounds for the QP subproblems
     * according to the information from current iterate
     *
     * @param delta      trust region radius
     * @param x_k 	     current iterate point
     * @param c_k        current constraint value evaluated at x_k
     * @param x_l        the lower bounds for variables
     * @param x_u        the upper bounds for variables
     * @param c_l        the lower bounds for constraints
     * @param c_u        the upper bounds for constraints
     */
    virtual void setup_bounds(double delta,
                              shared_ptr<const Vector> x_k,
                              shared_ptr<const Vector> x_l,
                              shared_ptr<const Vector> x_u,
                              shared_ptr<const Vector> c_k,
                              shared_ptr<const Vector> c_l,
                              shared_ptr<const Vector> c_u);

    /**
     * @brief This function sets up the object vector g
     * of the QP problem
     *
     * @param grad 	Gradient vector from nlp class
     * @param rho  	Penalty Parameter
     */
    void setup_g(shared_ptr<const Vector> grad, double rho);

    /**
     * Set up the H for the first time in the QP
     * problem.
     * It will be concatenated as [H_k 0]
     * 			                  [0   0]
     * where H_k is the Lagragian hessian evaluated at x_k
     * and lambda_k.
     *
     * This method should only be called for once.
     *
     * @param hessian the Lagragian hessian evaluated
     * at x_k and lambda_k from nlp readers.
     * @return
     */

    virtual void setup_H(shared_ptr<const SpTripletMat> hessian);

    /** @brief setup the matrix A for the QP subproblems according to the
     * information from current iterate*/
    virtual void setup_A(shared_ptr<const SpTripletMat> jacobian);


    /**
     * @brief solve the QP subproblem according to the bounds setup before,
     * assuming the first QP subproblem has been solved.
     * */

    void
    solveQP(shared_ptr<SQPhotstart::Stats> stats, shared_ptr<Options> options);


    /**
    * @brief This function updates the bounds on x if there is any changes to the
    * values of trust-region or the iterate
    */
    virtual void update_constraints(double delta,
                                    shared_ptr<const Vector> x_l,
                                    shared_ptr<const Vector> x_u,
                                    shared_ptr<const Vector> c_k,
                                    shared_ptr<const Vector> c_l,
                                    shared_ptr<const Vector> c_u,
                                    shared_ptr<const Vector> x_k);


    /**
     * @brief This function updates the vector g in the QP subproblem when there
     * are any change to the values of penalty parameter
     *
     * @param rho		penalty parameter
     */
    virtual void update_penalty(double rho);

    /**
     * @brief This function updates the vector g in the QP subproblem when there
     * are any change to the values of gradient in NLP
     *
     * @param grad		the gradient vector from NLP
     */
    virtual void update_grad(shared_ptr<const Vector> grad);

    /*  @brief Update the SparseMatrix H of the QP
     *  problems when there is any change to the
     *  true function Hessian
     *
     *  */
    virtual void update_H(shared_ptr<const SpTripletMat> Hessian);

    /**
     * @brief Update the Matrix H of the QP problems
     * when there is any change to the Jacobian to the constraints.
     */
    virtual void update_A(shared_ptr<const SpTripletMat> Jacobian);


    const shared_ptr<qpOASESInterface>& getQpInterface() const;


    inline QPReturnType GetStatus() {
        return (qp_interface_->status());
    }
private:
    /**
     * @brief allocate memory to class members except QP objects
     * */
    virtual void
    allocate(SQPhotstart::Index_info nlp_info, SQPhotstart::QPType qptype);

    //@{
    /** Copy Constructor */
    QPhandler(const QPhandler&);

    /** Overloaded Equals Operator */
    void operator=(const QPhandler&);
    //@}

    /**public class member*/

    /** QP problem will be in the following form
     * min 1/2x^T H x+ g^Tx
     * s.t. lbA <= A x <= ubA,
     *      lb  <=   x <= ub.
     *
     */

private:
    //bounds that can be represented as vectors
    Identity2Info I_info_A;
    Identity2Info I_info_H;
    Index_info nlp_info_;
    QPType qptype_;
    shared_ptr<qpOASESInterface> qp_interface_; /**<an interface to the standard QP
                                                         solver specified by the user*/
};


} // namespace SQPhotstart
#endif //SQPHOTSTART_QP_HPP_ 
