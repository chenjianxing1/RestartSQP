#ifndef SQPHOTSTART_QPHANDLER_HPP_
#define SQPHOTSTART_QPHANDLER_HPP_

#include "sqphot/SqpNlpBase.hpp"
#include "sqphot/Statistics.hpp"
#include "sqphot/Utils.hpp"

#include "sqphot/QPsolverInterface.hpp"

// TODO: We should get ride of this
#include "sqphot/SQPDebug.hpp"

#include "IpOptionsList.hpp"

#ifdef DEBUG
#ifdef COMPARE_QP_SOLVER
#include "sqphot/QOREInterface.hpp"
#include "sqphot/qpOASESInterface.hpp"
#endif
#endif

namespace SQPhotstart {

/** Class for storing which quantities in a QP need to be updated. */
class QpUpdateTracker
{
public:
  /** Constructor.  Initializes everything to false. */
  QpUpdateTracker()
   : update_gradient_(false)
   , update_penalty_parameter_(false)
   , update_bounds_(false)
   , update_trust_region_radius_(false)
   , update_jacobian_(false)
   , update_hessian_(false)
  {
  }

  /** Accessor methods */
  //@{
  bool need_gradient_update() const
  {
    return update_gradient_;
  }
  bool need_penalty_parameter_update() const
  {
    return update_penalty_parameter_;
  }
  bool need_bounds_update() const
  {
    return update_bounds_;
  }
  bool need_trust_region_radius_update() const
  {
    return update_trust_region_radius_;
  }
  bool need_jacobian_update() const
  {
    return update_jacobian_;
  }
  bool need_hessian_update() const
  {
    return update_hessian_;
  }
  //@}

  /** Set methods. */
  //@{
  void trigger_gradient_update()
  {
    update_gradient_ = true;
  }
  void trigger_penalty_parameter_update()
  {
    update_penalty_parameter_ = true;
  }
  void trigger_bounds_update()
  {
    update_bounds_ = true;
  }
  void trigger_trust_region_radius_update()
  {
    update_trust_region_radius_ = true;
  }
  void trigger_jacobian_update()
  {
    update_jacobian_ = true;
  }
  void trigger_hessian_update()
  {
    update_hessian_ = true;
  }
  //@}

  /** Mark all quantities as needed to be updated. */
  void trigger_all_updates()
  {
    update_gradient_ = true;
    update_penalty_parameter_ = true;
    update_bounds_ = true;
    update_trust_region_radius_ = true;
    update_jacobian_ = true;
    update_hessian_ = true;
  }

  /** Reset all updates to false. */
  void reset()
  {
    update_gradient_ = false;
    update_penalty_parameter_ = false;
    update_bounds_ = false;
    update_trust_region_radius_ = false;
    update_jacobian_ = false;
    update_hessian_ = false;
  }

  /** Check if any quantities need to be updated. */
  bool need_update()
  {
    return (update_gradient_ || update_penalty_parameter_ || update_bounds_ ||
            update_trust_region_radius_ || update_jacobian_ || update_hessian_);
  }

  ~QpUpdateTracker()
  {
  }

private:
  /** Hide default methods. */
  //@{
  /** Copy Constructor */
  QpUpdateTracker(const QpUpdateTracker&);
  /** Overloaded Equals Operator */
  void operator=(const QpUpdateTracker&);
  //@}

  /** Flag indicating whether the full gradient needs to be updated. */
  bool update_gradient_;
  /** Flag indicating whether only the penalty part of the gradient needs to be
   * updated. */
  bool update_penalty_parameter_;
  /** Flag indicating whether the bounds need to be updated. */
  bool update_bounds_;
  /** Flag indicating whether only the trust region radius needs to be updated.
   */
  bool update_trust_region_radius_;
  /** Flag indicating whether the Jacobian needs to be updated. */
  bool update_jacobian_;
  /** Flag indicating whether the Hessian needs to be updated. */
  bool update_hessian_;
};

/** Forward Declaration */

/**
 *
 * This is a class for setting up and solving the SQP
 * subproblems for Algorithm::Optimize
 *
 *  It contains the methods to setup the QP problem in
 *  the following format.
 *
 *  not NEW_FORMULATION:
 *
 * 	minimize  1/2 p_k^T H_k p^k + g_k^T p_k+rho*e^T *(u+v)
 * 	subject to c_l<=c_k+J_k p+u-v<=c_u,
 * 		       max(-delta,x_l) <= x_k+p_k <= min(delta,x_u)
 *		       u,v>=0
 *
 *  Order of variables: x = (p, u, v)
 *
 *  NEW_FORMULATION:
 *
 * 	minimize  1/2 p_k^T H_k p^k + g_k^T p_k+rho*[e^T *(u+v) + e^T(w+t)]
 * 	subject to c_l<=c_k+J_k p+u-v<=c_u,
 * 		       x_l<=x_k+p_k+w-t<=x_u,
 *		       -delta<=p_i<=delta,   (these are bounds)
 *		       u,v,w,t>=0
 *
 *  Order of variables: x = (p, u, v, w, t)
 *
 * and transform them into sparse matrix (triplet) and
 * dense vectors' format which can be taken as input
 * of standard QPhandler.
 *
 * It also contains a method which interfaces to the
 * QP solvers that users choose. The interface will pre
 * -process the data required by individual solver.
 */

class QPhandler
{

  ///////////////////////////////////////////////////////////
  //                      PUBLIC METHODS                   //
  ///////////////////////////////////////////////////////////
public:
  QPhandler(std::shared_ptr<const SqpNlpSizeInfo>, QPType qptype,
            Ipopt::SmartPtr<Ipopt::Journalist> jnlst,
            Ipopt::SmartPtr<const Ipopt::OptionsList> options);

  /** Default destructor */
  virtual ~QPhandler();

  /**
   * @brief Solve the QP subproblem according to the bounds setup before.
   * */
  QpSolverExitStatus solve_qp(std::shared_ptr<Statistics> stats);

  /** Solve the LP that has been set up before. */
  QpSolverExitStatus solve_lp(std::shared_ptr<Statistics> stats)
  {
    return solverInterface_->optimize_lp(stats);
  }
  /** @name Getters */
  //@{
  /**
   * @brief Get the optimal solution from the QPsolverinterface
   *
   * This is only an interface for user to avoid call interface directly.
   * @param p_k 	the pointer to an empty array with the length equal to
   * the
   * size
   * of the QP subproblem
   */
  std::shared_ptr<const Vector> get_primal_solution() const
  {
    return solverInterface_->get_primal_solution();
  }

  /**
   * @brief Get the infeasibility measure of the quadratic model
   * @return The one norm of the last (2*nConstr) varaibles of the QP solution
   */
  double get_model_infeasibility() const;
  /**
   *@brief Get the multipliers corresponding to the bound variables
   */
  std::shared_ptr<const Vector> get_bounds_multipliers() const
  {
    return solverInterface_->get_bounds_multipliers();
  }

  /**
   * @brief Get the multipliers corresponding to the constraints
   */
  std::shared_ptr<const Vector> get_constraints_multipliers() const
  {
    return solverInterface_->get_constraints_multipliers();
  }

  /**
   * @brief Get the objective value of the QP
   */
  double get_objective();

  /**
   * @brief manually calculate the active set from the class member
   * solverInterface
   * @param A_c  pointer to the active set corresponding to the constraints
   * @param A_b  pointer to the active set corresponding to the bound
   * constraints
   * @param x    solution for QP problem(optional)
   * @param Ax   constraint evaluation at current QP solution x(optional)
   */
  void get_active_set(ActivityStatus* A_c, ActivityStatus* A_b,
                      std::shared_ptr<Vector> x = nullptr,
                      std::shared_ptr<Vector> Ax = nullptr);

  /**
   * @brief Get the return status of QPsolver
   */
  Exitflag get_status();

  //@}

  /** @name Setters*/
  //@{
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
  void set_bounds(double delta, std::shared_ptr<const Vector> x_l,
                  std::shared_ptr<const Vector> x_u,
                  std::shared_ptr<const Vector> x_k,
                  std::shared_ptr<const Vector> c_l,
                  std::shared_ptr<const Vector> c_u,
                  std::shared_ptr<const Vector> c_k);

  /**
   * @brief This function sets up the object vector g
   * of the QP problem
   *
   * @param grad 	Gradient vector from nlp class
   * @param rho  	Penalty Parameter
   */
  void set_gradient(std::shared_ptr<const Vector> grad, double rho);

  void set_gradient(double rho);

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

  void set_hessian(std::shared_ptr<const SpTripletMat> hessian);

  /** @brief setup the matrix A for the QP subproblems according to the
   * information from current iterate*/
  void set_jacobian(std::shared_ptr<const SpTripletMat> jacobian);

  //@}

  /** @name Update QPdata */

  //@{
  void update_delta(double delta, std::shared_ptr<const Vector> x_l,
                    std::shared_ptr<const Vector> x_u,
                    std::shared_ptr<const Vector> x_k);

  /**
   * @brief This function updates the bounds on x if there is any changes to the
   * values of trust-region or the iterate
   */
  void update_bounds(double delta, std::shared_ptr<const Vector> x_l,
                     std::shared_ptr<const Vector> x_u,
                     std::shared_ptr<const Vector> x_k,
                     std::shared_ptr<const Vector> c_l,
                     std::shared_ptr<const Vector> c_u,
                     std::shared_ptr<const Vector> c_k);

  /**
   * @brief This function updates the vector g in the QP subproblem when there
   * are any change to the values of penalty parameter
   *
   * @param rho		penalty parameter
   */
  void update_penalty(double rho);

  /**
   * @brief This function updates the vector g in the QP subproblem when there
   * are any change to the values of gradient in NLP
   *
   * @param grad		the gradient vector from NLP
   */
  void update_grad(std::shared_ptr<const Vector> grad);

  /*  @brief Update the SparseMatrix H of the QP
   *  problems when there is any change to the
   *  true function Hessian
   *
   *  */
  void update_H(std::shared_ptr<const SpTripletMat> Hessian);

  /**
   * @brief Update the Matrix H of the QP problems
   * when there is any change to the Jacobian to the constraints.
   */
  void update_A(std::shared_ptr<const SpTripletMat> Jacobian);

  //@}

  /**
   * @brief Write data for the current QP to file.
   *
   * This creates a file that can be used to run the solver stand-alone for
   * debugging.
   */
  void write_qp_data(const std::string filename);

  /**
   * @brief Test the KKT conditions for the certain qpsolver
   */
  bool test_optimality(std::shared_ptr<QPSolverInterface> qpsolverInterface,
                       QpSolver qpSolver, ActivityStatus* W_b,
                       ActivityStatus* W_c);

  const OptimalityStatus& get_QpOptimalStatus() const;

#ifdef DEBUG
#ifdef COMPARE_QP_SOLVER

  void set_bounds_debug(double delta, std::shared_ptr<const Vector> x_k,
                        std::shared_ptr<const Vector> x_l,
                        std::shared_ptr<const Vector> x_u,
                        std::shared_ptr<const Vector> c_k,
                        std::shared_ptr<const Vector> c_l,
                        std::shared_ptr<const Vector> c_u);

  bool testQPsolverDifference();

#endif
#endif

  ///////////////////////////////////////////////////////////
  //                      PRIVATE METHODS                  //
  //////////////////////////////////////////////////////////

private:
  //@{

  /** Default constructor */
  QPhandler();

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

  ///////////////////////////////////////////////////////////
  //                      PRIVATE MEMBERS                  //
  ///////////////////////////////////////////////////////////

private:
  std::shared_ptr<QPSolverInterface>
      solverInterface_; /**<an interface to the standard
                             QP solver specified by the user*/

  IdentityMatrixPositions identity_matrix_positions_;

  /** Container with the sizes of the NLP */
  std::shared_ptr<const SqpNlpSizeInfo> nlp_sizes_;

  /** Journalist for output. */
  Ipopt::SmartPtr<Ipopt::Journalist> jnlst_;

  /** Number of variables in the QP. */
  int num_qp_variables_;
  /** Number of constraints in the QP. */
  int num_qp_constraints_;

  /** Working set for the bounds. */
  ActivityStatus* bounds_working_set_; // working set for bounds;
  /** Working set for the constraints. */
  ActivityStatus* constraints_working_set_; // working set for constraints;

  /** Indicates which solver is used to solve the QPs. */
  QpSolver qp_solver_choice_;

  OptimalityStatus qpOptimalStatus_;

#ifdef DEBUG
#ifdef COMPARE_QP_SOLVER
  std::shared_ptr<qpOASESInterface> qpOASESInterface_;
  std::shared_ptr<QOREInterface> QOREInterface_;
  ActivityStatus* W_c_qpOASES_; // working set for constraints;
  ActivityStatus* W_b_qpOASES_; // working set for bounds;
  ActivityStatus* W_c_qore_;    // working set for constraints;
  ActivityStatus* W_b_qore_;    // working set for bounds;
#endif
#endif
};

} // namespace SQPhotstart
#endif // SQPHOTSTART_QPHANDLER_HPP_
