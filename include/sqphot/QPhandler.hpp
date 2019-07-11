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
     * This is a class for setting up and solving the SQP subproblems for Algorithm::Optimize
     *
     * It contains the methods to setup the QP problem in the following format or similar ones,
     *
     * 	minimize  1/2 p_k^T H_k p^k + g_k^T p_k+rho*(u+v+w+z)
     * 	subject to c_l<=c_k+J_k p+u-v<=c_u,
     * 		  x_l<=x_k+p +w-z<=x_u,
     *		  -delta<=p_i<=delta,
     *		  w,v,w,z>=0
     *
     *
     * and transform them into sparse matrix (triplet) and dense vectors' format which can be
     * taken as input of standard qpSolver.
     *
     * It also contains a method which interfaces to the QP solvers that users choose. The
     * interface will pre-process the data required by individual solver.
     */

    class QPhandler {
    public:
        /** Default constructor */
        QPhandler();

        /** Default destructor */
        virtual ~QPhandler();


        /**
         * @name Get the optimal solution from the QPsolver_interface
         *
         *This is only an interface for user to avoid call interface directly.
         * @param p_k 	the pointer to an empty array with the length equal to the size of the QP subproblem
         */
        virtual bool GetOptimalSolution(double* p_k);


        /**
         *@name Get the multipliers from the QPsolver_interface
         *
         *This is only an interface for user to avoid call interface directly.
         *
         * @param y_k 	the pointer to an empty array with the length equal to the size of multipliers of the QP subproblem
         */
        virtual bool GetMultipliers(double* y_k);

        /**
         * @name This function initializes all objects will be used in this class.
         *
         * @param nlp_info
         * 		it contains the information about number of variables, number of constraints, number of
         * 		elements in the Hessian and that of Jacobian. The definition of Index_info is in Types.hpp
         * @param Constraint_type
         * 		it specifies how the variables are bounded. Please check @ClassifyConstraintType in Algorithm.hpp
         * 		for more details
         * @param qptype
         * 		it specifies which type of QP is going to be solved. It can be either LP, or QP, or SOC
         */
        virtual bool init(Index_info nlp_info, QPType qptype);

        /**
         *
         * @name setup the bounds for the QP subproblems according to the information from current iterate
         *
         * @param delta 	 trust region radius
         * @param x_k 	 current iterate point
         * @param c_k        current constraint value evaluated at x_k
         * @param x_l        the lower bounds for variables
         * @param x_u        the upper bounds for variables
         * @param c_l        the lower bounds for constraints
         * @param c_u        the upper bounds for constraints
         */
        virtual bool setup_bounds(double delta, shared_ptr<const Vector> x_k, shared_ptr<const Vector> x_l,
                                  shared_ptr<const Vector> x_u);

        /**
         * @name This function sets up the object vector g of the QP problem
         *
         * @param grad 	Gradient vector from nlp class
         * @param rho  	Penalty Parameter
         */
        virtual bool setup_g(shared_ptr<const Vector> grad_f, double rho);


        /**
         *
         * @param hessian
         * @return
         */

        virtual bool setup_H(shared_ptr<const SpMatrix> hessian);


        /** @name setup the matrix A for the QP subproblems according to the information from current iterate*/
        virtual bool setup_A(shared_ptr<const SpMatrix> jacobian);


        /**
         * @name solve the QP subproblem according to the bounds setup before, assuming the first QP
         * subproblem has been solved.
         * */

        virtual bool solveQP(shared_ptr<SQPhotstart::Stats> stats, shared_ptr<Options> options);

        /**
         * @name This function copies the information from another QPsolver objects.
         * By default, it will copy the vectors of bounds matrix A information
         *
         * If qptype ==LP, then it will not copy the Hessian, and the first half of the g)
         * objects in qpOASES
         */
        virtual bool
        copy_QP_info(shared_ptr<const QPhandler> rhs,// the QPhandler object that are going to be copied from
                     QPType qptype) { return false; };                // is the object rhs an LP?

        /**
         * @name This function updates the bounds on x if there is any changes to the values
	 * of trust-region radius
         *
         * @param delta 	 trust region radius
         * @param nVar 		 number of variables in NLP
         */
        virtual bool update_bounds(double delta, shared_ptr<const Vector> x_k, shared_ptr<const Vector> x_l,
                                   shared_ptr<const Vector> x_u); //the trust region radius




        /**
         * @name This function updates the vector g in the QP subproblem when there are any change to the values of penalty parameter
         *
         * @param rho		penalty parameter
         * @param nVar 		number of variables in NLP
         */
        virtual bool update_penalty(double rho);

        /**
         * @name This function updates the vector g in the QP subproblem when there
         * are any change to the values of gradient in NLP
         *
         * @param grad		the gradient vector from NLP
         */
        virtual bool update_grad(shared_ptr<const Vector> grad);

        /*  @name Update the SparseMatrix H of the QP problems when there is any change to the
         *  true function Hessian
         *
         *  */
        virtual bool update_H(shared_ptr<const SpMatrix> Hessian);

        /**
         * @name Update the Matrix H of the QP problems when there is any change to the
         *  Jacobian to the constraints.
         */
        virtual bool update_A(shared_ptr<const SpMatrix> Jacobian);

        inline double get_obj() { return qp_obj_; }

    private:
        /** allocate memory to class members except QP objects*/
        virtual bool allocate(SQPhotstart::Index_info nlp_info, SQPhotstart::QPType qptype);

        /**free all the memory*/
        virtual bool freeMemory() { return false; };

        /** Copy Constructor */
        QPhandler(const QPhandler &);

        /** Overloaded Equals Operator */
        void operator=(const QPhandler &);
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
        Index_info nlp_info_;
        QPType qptype_;
        shared_ptr<qpOASESInterface> qp_interface_; //an interface to the standard QP solver specified by the user
        double qp_obj_;        // the optimal objectives from QPhandler
    };


} // namespace SQPhotstart
#endif //SQPHOTSTART_QP_HPP_ 
