/* Copyright (C) 2019
 * All Rights Reserved.
 *
 * Authors: Xinyi Luo
 * Date:    2019-07
 */

#ifndef __QPOASES_INTERFACE_HPP__
#define __QPOASES_INTERFACE_HPP__

#include "IpOptionsList.hpp"
#include "qpOASES.hpp"
#include "sqphot/QpSolverInterface.hpp"
#include "sqphot/SqpNlpBase.hpp"

namespace RestartSqp {

/**
 * @brief This is a derived class of QPsolverInterface.
 * It uses qpOASES as the QP solver  which features the hotstart option. It is
 * used
 * as the default QP solver for SQPhostart.
 */
class QpOasesInterface : public QpSolverInterface
{

  ///////////////////////////////////////////////////////////
  //                      PUBLIC METHODS                   //
  ///////////////////////////////////////////////////////////
public:
  /**
   * @brief Constructor which also initializes the qpOASES SQProblem objects
   * @param nlp_info the struct that stores simple nlp dimension info
   * @param qptype  is the problem to be solved QP or LP?
   */
  QpOasesInterface(int num_qp_variables, int num_qp_constraints, QPType qptype,
                   Ipopt::SmartPtr<const Ipopt::OptionsList> options,
                   Ipopt::SmartPtr<Ipopt::Journalist> jnlst);

  /** This is only in unit test, TODO: remove later. */
  QpOasesInterface(std::shared_ptr<SparseHbMatrix> H,
                   std::shared_ptr<SparseHbMatrix> A, std::shared_ptr<Vector> g,
                   std::shared_ptr<Vector> lb, std::shared_ptr<Vector> ub,
                   std::shared_ptr<Vector> lbA, std::shared_ptr<Vector> ubA,
                   Ipopt::SmartPtr<const Ipopt::OptionsList> options = nullptr);

  /** Defualt Destructor */
  ~QpOasesInterface() override;

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
  QpOasesInterface();

  /** Copy Constructor */
  QpOasesInterface(const QpOasesInterface&);

  /** Overloaded Equals Operator */
  void operator=(const QpOasesInterface&);
  //@}

  /** Implementation of the optimziation method. */
  QpSolverExitStatus
  optimize_impl(std::shared_ptr<Statistics> stats = nullptr) override;

  /** Method for computing the working set from the most recent solve. */
  void retrieve_working_set_() override;

  /** Extract the options from the options list */
  void get_option_values_(Ipopt::SmartPtr<const Ipopt::OptionsList> options);

  /** Translate the qpOASES status into our exit status definition. */
  QpSolverExitStatus get_qpoases_exit_status_();

  void handle_error(QPType qptype, std::shared_ptr<Statistics> stats = nullptr);

  void set_qp_solver_options_();

  ///////////////////////////////////////////////////////////
  //                      PRIVATE MEMBERS                  //
  ///////////////////////////////////////////////////////////
private:
  /** qpOASES solver instance. */
  std::shared_ptr<qpOASES::SQProblem> qpoases_solver_;

  /** Flag indicating whether an instance of the QP has already been solved.  It
   *  is true if none has been solved. */
  bool first_qp_solved_;

  /** Flag indicating whether the matrices (Jacobian or Hessian) have changed
   *  since the last call to qpOASES. */
  bool qp_matrices_changed_;

  /** Matrix data for the QP. */
  //@{
  /** qpOASES version of the Hessian matrix.
   *
   *  IMPORTANT: This is a SHALLOW copy of our hessian_ matrix!
   */
  std::shared_ptr<qpOASES::SymSparseMat> qpoases_hessian_;

  /** qpOASES version of the Jacobian matrix.
   *
   *  IMPORTANT: This is a SHALLOW copy of our jacobian_ matrix!
   */
  std::shared_ptr<qpOASES::SparseMatrix> qpoases_jacobian_;

  /** For debugging, we store the number of nonzeros of the Jacobian matrix to
   * make sure it is not changed between solves. */
  int num_elements_jacobian_;

  /** For debugging, we store the number of nonzeros of the Hessian matrix to
   * make sure it is not changed between solves. */
  int num_elements_hessian_;
  //@}

  /** @name Algorithmic options */
  //@{
  /** Maximum number of QP iterations per QP solve. */
  int qp_solver_max_num_iterations_;
  /** Print level for QP solver */
  int qp_solver_print_level_;
  /** Maximum number of LP iterations per LP solve. */
  int lp_solver_max_num_iterations_;
  //@}
};
}
#endif
