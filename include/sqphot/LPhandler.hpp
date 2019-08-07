#ifndef SQPHOTSTART_LPHANDLER_HPP_
#define SQPHOTSTART_LPHANDLER_HPP_

#include <sqphot/QPhandler.hpp>
#include <sqphot/Utils.hpp>
#include <qpOASES.hpp>
#include <sqphot/Stats.hpp>
#include <sqphot/Options.hpp>
#include <sqphot/qpOASESInterface.hpp>

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
    LPhandler(Index_info nlp_info);

    /** Default destructor */
    ~LPhandler() override;

    void set_bounds(double delta,
                    shared_ptr<const Vector> x_k,
                    shared_ptr<const Vector> x_l,
                    shared_ptr<const Vector> x_u,
                    shared_ptr<const Vector> c_k,
                    shared_ptr<const Vector> c_l,
                    shared_ptr<const Vector> c_u) override;


    void set_g(double rho);;


    void set_A(shared_ptr<const SpTripletMat> jacobian) override;

    /**
     * @brief Get the optimal solution from the QPsolverinterface
     *
     * This is only an interface for user to avoid call interface directly.
     * @param p_k 	the pointer to an empty array with the length equal to the size
     * of the QP subproblem
     */
    void GetOptimalSolution(double* p_k) override;


    void update_delta(double delta,
                      shared_ptr<const Vector> x_l,
                      shared_ptr<const Vector> x_u,
                      shared_ptr<const Vector> x_k) override {
//        QPhandler::update_delta(delta,x_l,x_u,x_k);
    }

    const shared_ptr<QPSolverInterface>& getSolverInterface() const;;

    void update_bounds(double delta, shared_ptr<const Vector> x_l,
                       shared_ptr<const Vector> x_u, shared_ptr<const Vector> x_k,
                       shared_ptr<const Vector> c_l, shared_ptr<const Vector> c_u,
                       shared_ptr<const Vector> c_k) override {
//        QPhandler::update_bounds(delta,x_l,x_u,x_k,c_l,c_u,c_k);
    };

    void update_penalty(double rho) override {
//        QPhandler::update_penalty(rho);
    };


    void update_A(shared_ptr<const SpTripletMat> Jacobian) override {
//        QPhandler::update_A(Jacobian);
    };

    /**
     * @brief solve the QP subproblem according to the bounds setup before,
     * assuming the first QP subproblem has been solved.
     * */

    void solveLP(shared_ptr<SQPhotstart::Stats> stats, shared_ptr<Options> options);


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
    Index_info nlp_info_;
    shared_ptr<QPSolverInterface> solverInterface_; //an interface to the standard LP
    bool isAinitialised = false;//TODO: delete it later
};


} // namespace SQPhotstart
#endif //SQPHOTSTART_QP_HPP_ 
