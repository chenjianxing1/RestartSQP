#ifndef SQPHOTSTART_LPHANDLER_HPP_
#define SQPHOTSTART_LPHANDLER_HPP_

#include <sqphot/QPhandler.hpp>
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

    class LPhandler : public QPhandler {
    public:
        /** Default constructor */
        LPhandler();

        /** Default destructor */
        virtual ~LPhandler();


        /**
         * @brief Get the optimal solution from the QPsolverinterface
         *
         * This is only an interface for user to avoid call interface directly.
         * @param p_k 	the pointer to an empty array with the length equal to the size
         * of the QP subproblem
         */
        virtual bool GetOptimalSolution(double* p_k);


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
         */
        bool init(Index_info nlp_info, QPType qptype);

        /**
         *
         * @brief setup the bounds for the QP subproblems according to the information
         * from current iterate x_k
         *
         * @param delta      trust region radius
         * @param x_k 	     current iterate point
         * @param c_k        current constraint value evaluated at x_k
         * @param x_l        the lower bounds for variables
         * @param x_u        the upper bounds for variables
         * @param c_l        the lower bounds for constraints
         * @param c_u        the upper bounds for constraints
         */
        bool setup_bounds(double delta,
                          shared_ptr<const Vector> x_k,
                          shared_ptr<const Vector> x_l,
                          shared_ptr<const Vector> x_u,
                          shared_ptr<const Vector> c_k,
                          shared_ptr<const Vector> c_l,
                          shared_ptr<const Vector> c_u) override;


        /**
         * @brief This function sets up the object vector g
         * of the QP problem
         *
         * @param grad 	Gradient vector from nlp class
         * @param rho  	Penalty Parameter
         */
        bool setup_g(double rho);

        /** @brief setup the matrix A for the QP subproblems according to the
         * information from current iterate*/
        bool setup_A(shared_ptr<const SpTripletMat> jacobian) override;


        /**
         * @brief solve the QP subproblem according to the bounds setup before,
         * assuming the first QP subproblem has been solved.
         * */

        bool solveLP(shared_ptr<SQPhotstart::Stats> stats, shared_ptr<Options> options);

        /**
         * @brief This function copies the information from another LPhandler object
         * By default, it will copy the vectors of bounds matrix A information
         *
	     * If the class member qptype_ ==LP, then it will not copy the Hessian, and the
         * first half of the g)object
         *
         * @param a LPhandler object own a QPsolverInterface specific to QP.
         */
        virtual bool copy_QP_info(shared_ptr<const QPhandler> rhs);


        /**
         * @brief This function updates the bounds on x if there is any changes to the
	     * values of trust-region
         * radius
         *
         * @param delta 	 trust region radius
         * @param nVar 		 number of variables in NLP
         */
        virtual bool
        update_constraints(double delta, shared_ptr<const Vector> x_l,
                           shared_ptr<const Vector> x_u, shared_ptr<const Vector> c_k,
                           shared_ptr<const Vector> c_l, shared_ptr<const Vector> c_u,
                           shared_ptr<const Vector> x_k); //the trust region radius




        /**
         * @brief This function updates the vector g in the
         * QP subproblem when there are any change to the values
         * of penalty parameter
         *
         * @param rho		penalty parameter
         * @param nVar 		number of variables in NLP
         */
        bool update_penalty(double rho) override;

        /**
         * @brief This function updates the vector g in the
         * QP subproblem when there are any change to the values
         * of gradient in NLP
         *
         * @param grad		the gradient vector from NLP
         */
        bool update_grad(shared_ptr<const Vector> grad) override;


        /**
         * @brief Update the Matrix H of the QP problems
         * when there is any change to the Jacobian to the constraints.
         */
        bool update_A(shared_ptr<const SpTripletMat> Jacobian) override;

    private:
        /**
         * @brief allocate memory to class members except QP objects
         * */
        virtual bool
        allocate(SQPhotstart::Index_info nlp_info, SQPhotstart::QPType qptype);

        /**free all the memory*/
        virtual bool freeMemory() { return false; };

        /** Copy Constructor */
        LPhandler(const LPhandler&);

        /** Overloaded Equals Operator */
        void operator=(const LPhandler&);
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
        shared_ptr<qpOASESInterface> lp_interface_; //an interface to the standard LP

        // solver specified by the user
        double lp_obj_;        // the optimal objectives from LPhandler
    };


} // namespace SQPhotstart
#endif //SQPHOTSTART_QP_HPP_ 
