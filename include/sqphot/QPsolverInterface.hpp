/* Copyright (C) 2019
 * All Rights Reserved.
 *
 * Authors: Xinyi Luo      2019-07
 */
#ifndef SQPHOTSTART_QPSOLVER_INTERFACE_HPP
#define SQPHOTSTART_QPSOLVER_INTERFACE_HPP

#include "IpException.hpp"
#include "sqphot/SQPDebug.hpp"
#include "sqphot/SparseHbMatrix.hpp"
#include "sqphot/Statistics.hpp"
#include "sqphot/Vector.hpp"
#include <memory>

namespace SQPhotstart {

enum QpSolverExitStatus
{
  QPEXIT_OPTIMAL = 0,
  QPEXIT_UNBOUNDED = 1,
  QPEXIT_INFEASIBLE = 2,
  QPEXIT_UNKNOWN_STATUS = -96,
  QPEXIT_FAILED = -97,
  QPEXIT_NOT_SOLVED = -98,
  QPEXIT_INTERNAL_ERROR = -99
};

#if 0
DECLARE_STD_EXCEPTION(QP_NOT_OPTIMAL);

DECLARE_STD_EXCEPTION(LP_NOT_OPTIMAL);

DECLARE_STD_EXCEPTION(QPSOLVER_INTERNAL_ERROR);

DECLARE_STD_EXCEPTION(INVALID_WORKING_SET);
#endif
/**
 * @brief Virtual base class for all standard QP solvers that use standard
 * triplet matrix
 * form and dense vectors.
 *
 * It can optimize QP problem in the following format
 *
 *  minimize 1/2 x^T H x + g^T x
 *  subject  lb_A<=Ax<=ub_A
 *              lb<=x<=ub
 */
class QPSolverInterface
{
public:
  /**@name Getters for current QP data. */
  virtual std::shared_ptr<const Vector> get_lower_variable_bounds() const = 0;

  virtual std::shared_ptr<const Vector> get_upper_variable_bounds() const = 0;

  virtual std::shared_ptr<const Vector> get_lower_constraint_bounds() const = 0;

  virtual std::shared_ptr<const Vector> get_upper_constraint_bounds() const = 0;

  virtual std::shared_ptr<const Vector> get_linear_objective_coefficients() const = 0;

  virtual std::shared_ptr<const SparseHbMatrix> get_objective_hessian() const = 0;

  virtual std::shared_ptr<const SparseHbMatrix> get_constraint_jacobian() const = 0;
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
  virtual QpSolverExitStatus optimize_qp(std::shared_ptr<Statistics> stats) = 0;

  /**
   * @brief Solve a regular LP with given data and options
   *
   * overload this method to optimize a LP with the data specified, update the
   * stats by adding the iteration number used to solve this LP to stats.qp_iter
   */

  virtual QpSolverExitStatus optimize_lp(std::shared_ptr<Statistics> stats) = 0;

  /**-------------------------------------------------------**/
  /**                    Getters                            **/
  /**-------------------------------------------------------**/
  /**@name Getters*/
  //@{
  /**
   * @return the pointer to the optimal solution
   *
   */
  virtual std::shared_ptr<const Vector> get_primal_solution() const = 0;

  /**
   *@brief get the objective value from the QP solvers
   *
   * @return the objective function value of the QP problem
   */
  virtual double get_obj_value() = 0;

  /**
   * @brief get the pointer to the multipliers to the bounds constraints.
   */
  virtual std::shared_ptr<const Vector> get_bounds_multipliers() const = 0;

  /**
   * @brief get the pointer to the multipliers to the regular constraints.
   */
  virtual std::shared_ptr<const Vector> get_constraints_multipliers() const = 0;

  /**
   * @brief copy the working set information
   * @param W_constr a pointer to an array of length (nCon_QP_) which will store
   * the
   * working set for constraints
   * @param W_bounds a pointer to an array of length (nVar_QP_) which will store
   * the
   * working set for bounds
   *
   * overload this method by assign each entry of W_constr and W_bounds to be
   * ACTIVE_ABOVE, ACTIVE_BELOW, INACTIVE, or ACTIVE_BOTH_SIDE based on the
   * working
   * set information from the QP solver
   */
  virtual void get_working_set(ActivityStatus* W_constr,
                               ActivityStatus* W_bounds) = 0;

  virtual Exitflag get_status() = 0;

  virtual bool test_optimality(ActivityStatus* W_c = NULL,
                               ActivityStatus* W_b = NULL) = 0;

  virtual OptimalityStatus get_optimality_status() = 0;

  //@}

  /**-------------------------------------------------------**/
  /**                    Setters                            **/
  /**-------------------------------------------------------**/
  /**@name Setters, by location and value*/
  //@{
  virtual void set_lower_variable_bounds(int location, double value) = 0;

  virtual void set_upper_variable_bounds(int location, double value) = 0;

  virtual void set_lower_constraint_bounds(int location, double value) = 0;

  virtual void set_upper_constraint_bounds(int location, double value) = 0;

  virtual void set_linear_objective_coefficients(int location, double value) = 0;
  //@}

  /**@name Setters for dense vector, by vector value*/
  //@{
  virtual void set_upper_variable_bounds(std::shared_ptr<const Vector> rhs) = 0;

  virtual void set_lower_variable_bounds(std::shared_ptr<const Vector> rhs) = 0;

  virtual void set_lower_constraint_bounds(std::shared_ptr<const Vector> rhs) = 0;

  virtual void set_upper_constraint_bounds(std::shared_ptr<const Vector> rhs) = 0;

  virtual void set_linear_objective_coefficients(std::shared_ptr<const Vector> rhs) = 0;
  //@}

  /**@name Setters for matrix*/
  //@{
  virtual void set_objective_hessian(std::shared_ptr<const SpTripletMat> rhs) = 0;

  virtual void
  set_constraint_jacobian(std::shared_ptr<const SpTripletMat> rhs,
               IdentityMatrixPositions& identity_matrix_positions) = 0;
  //@}

  virtual void reset_constraints() = 0;

  /**-------------------------------------------------------**/
  /**                  Data Writer                          **/
  /**-------------------------------------------------------**/

  virtual void WriteQPDataToFile(Ipopt::EJournalLevel level,
                                 Ipopt::EJournalCategory category,
                                 const std::string filename) = 0;

protected:
private:
  /** Copy Constructor */
  QPSolverInterface(const QPSolverInterface&);

  /** Overloaded Equals Operator */
  void operator=(const QPSolverInterface&);
};

} // SQPHOTSTART
#endif
