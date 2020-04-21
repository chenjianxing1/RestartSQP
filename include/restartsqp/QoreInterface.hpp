/* Copyright (C) 2019
 * All Rights Reserved.
 *
 * Authors: Xinyi Luo
 * Date:    2019-08
 */
#ifndef __SQPHOTSTART_QOREINTERFACE_HPP__
#define __SQPHOTSTART_QOREINTERFACE_HPP__

extern "C" {
#include "qpsolver.h"
}

#include "IpOptionsList.hpp"
#include "restartsqp/MessageHandling.hpp"
#include "restartsqp/QpSolverInterface.hpp"

DECLARE_STD_EXCEPTION(INVALID_RETURN_TYPE);

namespace RestartSqp {

class QoreInterface : public QpSolverInterface
{

  ///////////////////////////////////////////////////////////
  //                      PUBLIC METHODS                   //
  ///////////////////////////////////////////////////////////
public:
  /**Constructor*/

  QoreInterface(int num_qp_variables, int num_qp_constraints, QPType qp_type,
                Ipopt::SmartPtr<const Ipopt::OptionsList> options,
                Ipopt::SmartPtr<Ipopt::Journalist> jnlst);

  /** Default destructor*/
  ~QoreInterface() override;

  /** Set the Hessian matrix from a triplet matrix. */
  void set_objective_hessian(
      std::shared_ptr<const SparseTripletMatrix> triplet_matrix) override;

  /** Set the Jacobian matrix from a triplet matrix. */
  void set_constraint_jacobian(
      std::shared_ptr<const SparseTripletMatrix> triplet_matrix,
      IdentityMatrixPositions& identity_matrix_positions) override;

  ///////////////////////////////////////////////////////////
  //                      PRIVATE METHODS                  //
  ///////////////////////////////////////////////////////////
private:
  /** @name Compiler generated methods to hide. */
  //@{
  /** Default Constructor*/
  QoreInterface();

  /** Copy Constructor */
  QoreInterface(const QoreInterface&);

  /** Overloaded Equals Operator */
  void operator=(const QoreInterface&);
  //@}

  /** Create a QORE solver object. */
  void create_qore_solver_(int num_nnz_jacobian, int num_nnz_hessian);

  /** Implementation of the optimziation method. */
  QpSolverExitStatus
  optimize_impl(std::shared_ptr<Statistics> stats = nullptr) override;

  /** Method for computing the working set from the most recent solve. */
  void retrieve_working_set_() override;

  /** Extract the options from the Ipopt options list. */
  void get_option_values_(Ipopt::SmartPtr<const Ipopt::OptionsList> options);

  /** Translate QORE exit status into our QP status. */
  QpSolverExitStatus get_qore_exit_status_();

  /**
   * @brief set options of QP solver based on the user-defined values
   */
  void set_qp_solver_options_();

  /**
   * @brief Handle errors based on current status
   */
  void handle_error(QPType qptype, std::shared_ptr<Statistics> stats = nullptr);

  ///////////////////////////////////////////////////////////
  //                      PRIVATE MEMBERS                  //
  ///////////////////////////////////////////////////////////
private:
  int status_;

  /** QORE solver object. */
  QoreProblem* qore_solver_;

  /** Flag indicating whether an instance of the QP has already been solved.  It
   *  is true if none has been solved. */
  bool first_qp_solved_;

  /** Flag indicating whether the matrices (Jacobian or Hessian) have changed
   *  since the last call to qpOASES. */
  bool qp_matrices_changed_;

  /** Most recent primal solution of QORE. */
  double* qore_primal_solution_;

  /** Most recent dual solution of QORE. */
  double* qore_dual_solution_;

  /** @name Algorithmic options */
  //@{
  // Maximum number of QP iterations per QP solve
  int qp_solver_max_num_iterations_;
  // Print level for QP solver
  int qp_solver_print_level_;
  // Maximum number of LP iterations per LP solve
  int lp_solver_max_num_iterations_;
  // Flag indicating whether the primal variables should be initialized
  // to zero after the first iteration
  bool qore_init_primal_variables_;
  // Regularization factor for the QP Hessian
  double qore_hessian_regularization_;
  //@}
};
}

#endif
