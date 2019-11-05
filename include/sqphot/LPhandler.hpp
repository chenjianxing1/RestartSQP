#ifndef SQPHOTSTART_LPHANDLER_HPP_
#define SQPHOTSTART_LPHANDLER_HPP_

#include <qpOASES.hpp>
#include <sqphot/Options.hpp>
#include <sqphot/QPhandler.hpp>
#include <sqphot/SQPTNLP.hpp>
#include <sqphot/Stats.hpp>
#include <sqphot/Utils.hpp>

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

    ///////////////////////////////////////////////////////////
    //                      PUBLIC METHODS                   //
    ///////////////////////////////////////////////////////////

public:

    LPhandler(NLPInfo nlp_info, shared_ptr<const Options> options,
              Ipopt::SmartPtr<Ipopt::Journalist> jnlst);

    /** Default destructor */
    ~LPhandler() override;

    /**
     * @brief solve the QP subproblem according to the bounds setup before,
     * assuming the first QP subproblem has been solved.
     * */

    void solveLP(shared_ptr<SQPhotstart::Stats> stats);

    /**@name Setters */
    //@{
//    void set_bounds(double delta, shared_ptr<const Vector> x_l,
//                    shared_ptr<const Vector> x_u, shared_ptr<const Vector> x_k,
//                    shared_ptr<const Vector> c_l, shared_ptr<const Vector> c_u,
//                    shared_ptr<const Vector> c_k);


    void set_g(double rho);;

//    void set_A(shared_ptr<const SpTripletMat> jacobian) override;
    //@}

    /**
     * @brief Get the optimal solution from the QPsolverinterface
     *
     * This is only an interface for user to avoid call interface directly.
     * @param p_k 	the pointer to an empty array with the length equal to the size
     * of the QP subproblem
     */

    ///////////////////////////////////////////////////////////
    //                      PRIVATE METHODS                  //
    ///////////////////////////////////////////////////////////

private:
    /**
     * @brief allocate memory to class members except QP objects
     * */

    /** Default constructor */
    LPhandler();


    /** Copy Constructor */
    LPhandler(const LPhandler &);

    /** Overloaded Equals Operator */
    void operator=(const LPhandler &);

    ///////////////////////////////////////////////////////////
    //                      PRIVATE MEMBERS                  //
    ///////////////////////////////////////////////////////////

    /** LP problem will be in the following form
     * min g^Tx
     * s.t. lbA <= A x <= ubA,
     *      lb  <=   x <= ub.
     *
     */
private:

    Solver LPsolverChoice_;
    //bounds that can be represented as vectors
    IdentityInfo I_info_A_;
    const NLPInfo nlp_info_;
    int nConstr_LP_;
    int nVar_LP_;
    Ipopt::SmartPtr<Ipopt::Journalist> jnlst_;
//    shared_ptr<QPSolverInterface> solverInterface_; /**<an interface to the standard
//                                                              QP solver specified by the user*/
    //bounds that can be represented as vectors
    bool isinitialised = false;
};


} // namespace SQPhotstart
#endif //SQPHOTSTART_QP_HPP_ 
