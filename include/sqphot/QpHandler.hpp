#ifndef SQPHOTSTART_QPHANDLER_HPP_
#define SQPHOTSTART_QPHANDLER_HPP_

#include "sqphot/SqpNlpBase.hpp"
#include "sqphot/Statistics.hpp"
#include "sqphot/Utils.hpp"

#include "sqphot/QpSolverInterface.hpp"

// TODO: We should get ride of this
#include "sqphot/SQPDebug.hpp"

#include "IpOptionsList.hpp"

#ifdef DEBUG
#ifdef COMPARE_QP_SOLVER
#include "sqphot/QoreInterface.hpp"
#include "sqphot/QpOasesInterface.hpp"
#endif
#endif

namespace RestartSqp {

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
  bool need_trust_region_radius_decrease() const
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
  void trigger_trust_region_radius_decrease()
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
 *  Formulation without slacks for the variables:
 *
 * 	minimize  1/2 p^T H_k p + g_k^T p+rho*e^T *(u+v)
 * 	subject to c_l<=c_k+J_k p+u-v<=c_u,
 * 		       max(-delta,x_l-x_k) <= p <= min(delta,x_u-x_k)
 *		       u,v>=0
 *
 *  Order of variables: x = (p, u, v)
 *
 *  Formulation with slacks for the variables:
 *
 * 	minimize  1/2 p^T H_k p + g_k^T p+rho*[e^T *(u+v) + e^T(w+t)]
 * 	subject to c_l<=c_k+J_k p+u-v<=c_u,
 * 		       x_l<=x_k+p+w-t<=x_u,
 *		       -delta<=p<=delta,   (these are bounds)
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

class QpHandler
{

  ///////////////////////////////////////////////////////////
  //                      PUBLIC METHODS                   //
  ///////////////////////////////////////////////////////////
public:
  QpHandler(std::shared_ptr<const SqpNlpSizeInfo>, QPType qptype,
            Ipopt::SmartPtr<Ipopt::Journalist> jnlst,
            Ipopt::SmartPtr<const Ipopt::OptionsList> options);

  /** Default destructor */
  virtual ~QpHandler();

  /**
   * @brief Solve the QP/LP subproblem according to the data set up before.
   * */
  QpSolverExitStatus solve(std::shared_ptr<Statistics> stats);

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
    return qp_solver_interface_->get_primal_solution();
  }

  /**
   * @brief Get the infeasibility measure of the quadratic model
   * @return The one norm of the last (2*nConstr) varaibles of the QP solution
   */
  double get_model_infeasibility() const;
  /**
   *@brief Get the multipliers corresponding to the NLP bound variables
   */
  std::shared_ptr<const Vector> get_bounds_multipliers() const;

  /**
   * @brief Get the multipliers corresponding to the NLP constraints
   */
  std::shared_ptr<const Vector> get_constraint_multipliers() const;

  /**
   * @brief Get the objective value of the QP
   */
  double get_qp_objective() const;

  /**
   * @brief Get the most recent KKT error
   */
  double get_qp_kkt_error() const;

#if 0
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
#endif
#if 0
  /**
   * @brief Get the return status of QPsolver
   */
  Exitflag get_status();
#endif
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
  void set_bounds(double trust_region_radius, std::shared_ptr<const Vector> x_l,
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
  void set_linear_qp_objective_coefficients(std::shared_ptr<const Vector> gradient, double penalty_parameter);
#if 0
  void set_linear_qp_objective_coefficients(double rho);
#endif

#if 0
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

  void set_hessian(std::shared_ptr<const SparseTripletMatrix> hessian);

  /** @brief setup the matrix A for the QP subproblems according to the
   * information from current iterate*/
  void set_jacobian(std::shared_ptr<const SparseTripletMatrix> jacobian);

  //@}
#endif
  /** @name Update QPdata */

  //@{
  void decrease_trust_region(double trust_region);

#if 0
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
#endif
  /**
   * @brief This function updates the vector g in the QP subproblem when there
   * are any change to the values of penalty parameter
   *
   * @param rho		penalty parameter
   */
  void update_penalty_parameter(double penalty_parameter);

#if 1
  /**
   * @brief This function updates the vector g in the QP subproblem when there
   * are any change to the values of gradient in NLP
   *
   * @param grad		the gradient vector from NLP
   */
  void set_linear_qp_objective_coefficients(std::shared_ptr<const Vector> grad);
#endif

  void set_linear_qp_objective_coefficients_to_zero()
  {
    qp_solver_interface_->get_linear_objective_coefficients_nonconst()->set_to_zero();
  }

  /*  @brief Update the SparseMatrix H of the QP
   *  problems when there is any change to the
   *  true function Hessian
   *
   *  */
  void set_hessian(std::shared_ptr<const SparseTripletMatrix> Hessian);

  /**
   * @brief Update the Matrix H of the QP problems
   * when there is any change to the Jacobian to the constraints.
   */
  void set_jacobian(std::shared_ptr<const SparseTripletMatrix> Jacobian);

  //@}

  /**
   * @brief Write data for the current QP to file.
   *
   * This creates a file that can be used to run the solver stand-alone for
   * debugging.
   */
  void write_qp_data(const std::string& filename);

#if 0
  /**
   * @brief Test the KKT conditions for the certain qpsolver
   */
  bool test_optimality(std::shared_ptr<QpSolverInterface> qpsolverInterface,
                       QpSolver qpSolver, ActivityStatus* W_b,
                       ActivityStatus* W_c);
#endif

  const KktError& get_QpOptimalStatus() const;

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
  QpHandler();

  /** Copy Constructor */
  QpHandler(const QpHandler&);

  /** Overloaded Equals Operator */
  void operator=(const QpHandler&);
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
  std::shared_ptr<QpSolverInterface>
      qp_solver_interface_; /**<an interface to the standard
                             QP solver specified by the user*/

  /** Flag indicating whether this is the formulation that permits the bounds to be violated. */
  bool slack_formulation_;

  /** Object that stores some meta-structure of the constraint Jacobian, namely
   * cthe position of identiy matrices for slack variables. */
  IdentityMatrixPositions identity_matrix_positions_;

  /** Journalist for output. */
  Ipopt::SmartPtr<Ipopt::Journalist> jnlst_;

  /** Number of variables in the NLP. */
  int num_nlp_variables_;
  /** Number of constraints in the NLP. */
  int num_nlp_constraints_;

  /** Number of variables in the QP. */
  int num_qp_variables_;
  /** Number of constraints in the QP. */
  int num_qp_constraints_;

  /** Indicates which solver is used to solve the QPs. */
  QpSolver qp_solver_choice_;

  /** Store the value for the penalty parameter used most recently */
  double last_penalty_parameter_;

  /** KKT error (worst violation) from the most recent QP/LP solve. */
  double last_kkt_error_;

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
