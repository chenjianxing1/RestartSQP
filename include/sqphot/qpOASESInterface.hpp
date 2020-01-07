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
#include "sqphot/QPsolverInterface.hpp"
#include "sqphot/SqpNlpBase.hpp"

namespace SQPhotstart {

enum QPMatrixType
{
  UNDEFINED,
  FIXED,
  VARIED
};

/**
 * @brief This is a derived class of QPsolverInterface.
 * It uses qpOASES as the QP solver  which features the hotstart option. It is
 * used
 * as the default QP solver for SQPhostart.
 */
class qpOASESInterface : public QPSolverInterface
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
  qpOASESInterface(std::shared_ptr<const SqpNlpSizeInfo> nlp_sizes,
                   QPType qptype,
                   Ipopt::SmartPtr<const Ipopt::OptionsList> options,
                   Ipopt::SmartPtr<Ipopt::Journalist>
                       jnlst); // number of constraints in the QP problem

  /** This is only in unit test, TODO: remove later. */
  qpOASESInterface(std::shared_ptr<SparseHbMatrix> H,
                   std::shared_ptr<SparseHbMatrix> A, std::shared_ptr<Vector> g,
                   std::shared_ptr<Vector> lb, std::shared_ptr<Vector> ub,
                   std::shared_ptr<Vector> lbA, std::shared_ptr<Vector> ubA,
                   Ipopt::SmartPtr<const Ipopt::OptionsList> options = nullptr);

  /** Defualt Destructor */
  ~qpOASESInterface() override;

  QpSolverExitStatus
  optimize_qp(std::shared_ptr<Statistics> stats = nullptr) override;

  /**
   * @brief optimize the LP problem whose objective and constraints are defined
   * in the class members.
   */

  QpSolverExitStatus
  optimize_lp(std::shared_ptr<Statistics> stats = nullptr) override;

  /** @name Getters*/
  //@{
  /**
  * @brief copy the optimal solution of the QP to the input pointer
  */

  std::shared_ptr<const Vector> get_primal_solution() const override;

  /**
   * @brief get the pointer to the multipliers to the bounds constraints.
   */
  std::shared_ptr<const Vector> get_bounds_multipliers() const override;
  /**
   * @brief get the pointer to the multipliers to the regular constraints.
   */
  std::shared_ptr<const Vector> get_constraints_multipliers() const override;

  void get_working_set(ActivityStatus* W_constr,
                       ActivityStatus* W_bounds) override;

  /**
   *@brief get the objective value from the QP solvers
   *
   * @return the objective function value of the QP problem
   */

  double get_obj_value() override;

  /**
   * @brief get the final return status of the QP problem
   */

  Exitflag get_status() override;

  OptimalityStatus get_optimality_status() override
  {
    return qpOptimalStatus_;
  }
  //@}

  /** @name Setters */
  //@{
  void set_lb(int location, double value) override;

  void set_lb(std::shared_ptr<const Vector> rhs) override;

  void set_ub(int location, double value) override;

  void set_ub(std::shared_ptr<const Vector> rhs) override;

  void set_lbA(int location, double value) override;

  void set_lbA(std::shared_ptr<const Vector> rhs) override;

  void set_ubA(int location, double value) override;

  void set_ubA(std::shared_ptr<const Vector> rhs) override;

  void set_gradient(int location, double value) override;

  void set_gradient(std::shared_ptr<const Vector> rhs) override;

  void set_hessian(std::shared_ptr<const SpTripletMat> rhs) override;

  void
  set_jacobian(std::shared_ptr<const SpTripletMat> rhs,
               IdentityMatrixPositions& identity_matrix_positions) override;

  //@}

  void WriteQPDataToFile(Ipopt::EJournalLevel level,
                         Ipopt::EJournalCategory category,
                         const std::string filename) override;

  void reset_constraints() override;

  //@{
  const std::shared_ptr<Vector>& getLb() const override
  {
    return lb_;
  };

  const std::shared_ptr<Vector>& getUb() const override
  {
    return ub_;
  };

  const std::shared_ptr<Vector>& getLbA() const override
  {
    return lbA_;
  };

  const std::shared_ptr<Vector>& getUbA() const override
  {
    return ubA_;
  };

  const std::shared_ptr<Vector>& getG() const override
  {
    return g_;
  };

  std::shared_ptr<const SparseHbMatrix> getH() const override
  {
    return H_;
  };

  std::shared_ptr<const SparseHbMatrix> getA() const override
  {
    return A_;
  };
  //@}

  bool test_optimality(ActivityStatus* W_c = NULL,
                       ActivityStatus* W_b = NULL) override;

  ///////////////////////////////////////////////////////////
  //                      PRIVATE METHODS                  //
  ///////////////////////////////////////////////////////////
private:
  /** default constructor*/
  qpOASESInterface();

  /** Extract the options from the options list */
  void get_option_values_(Ipopt::SmartPtr<const Ipopt::OptionsList> options);

  /** Translate the qpOASES status into our exit status definition. */
  QpSolverExitStatus
  get_qpoases_exit_status_(std::shared_ptr<qpOASES::SQProblem> qpoases_solver);

  void handle_error(QPType qptype, std::shared_ptr<Statistics> stats = nullptr);

  void reset_flags();

  /**
   * @brief Allocate memory for the class members
   * @param nlp_info  the struct that stores simple nlp dimension info
   * @param qptype is the problem to be solved QP or LP?
   */
  void allocate_memory(std::shared_ptr<const SqpNlpSizeInfo>, QPType qptype);

  void set_qp_solver_options_();

  void get_Matrix_change_status();

  /** Copy Constructor */
  qpOASESInterface(const qpOASESInterface&);

  /** Overloaded Equals Operator */
  void operator=(const qpOASESInterface&);

  ///////////////////////////////////////////////////////////
  //                      PRIVATE MEMBERS                  //
  ///////////////////////////////////////////////////////////
private:
  /** Journalist for output. */
  Ipopt::SmartPtr<Ipopt::Journalist> jnlst_;

  QPMatrixType new_QP_matrix_status_ = UNDEFINED;
  QPMatrixType old_QP_matrix_status_ = UNDEFINED;

  /** Flag indicating whether the matrices (Jacobian or Hessian) have changed
   *  since the last call to qpOASES. */
  bool qp_matrices_changed_;

  /** Flag indicating whether an instance of the QP has already been solved.  It
   * is true if none has been solved. */
  bool first_qp_solved_;

  /** Number of variables in the QP. */
  int num_qp_variables_;
  /** Number of constraints in the QP. */
  int num_qp_constraints_;

  /** Number of variables in the NLP. */
  int num_nlp_variables_;
  /** Number of constraints in the NLP. */
  int num_nlp_constraints_;

  OptimalityStatus qpOptimalStatus_;

  std::shared_ptr<qpOASES::SymSparseMat>
      H_qpOASES_; /**< the Matrix object that qpOASES
                   * taken in, it only contains the
                   * pointers to array stored in
                   * the class members of H_*/

  /** qpOASES solver instance. */
  std::shared_ptr<qpOASES::SQProblem> qpoases_solver_;

  std::shared_ptr<qpOASES::SparseMatrix>
      A_qpOASES_;                /**< the Matrix object that qpOASES
                                  * taken in, it only contains the
                                  * pointers to array stored in
                                  * the class members of A_*/
  std::shared_ptr<Vector> g_;    /**< the grad used for QPsubproblem*/
  std::shared_ptr<Vector> lbA_;  /**< lower bounds of Ax */
  std::shared_ptr<Vector> lb_;   /**< lower bounds of x */
  std::shared_ptr<Vector> ubA_;  /**< upper bounds of Ax */
  std::shared_ptr<Vector> ub_;   /**< upper bounds of x */
  std::shared_ptr<Vector> x_qp_; /** the qp solution */

  /** the bounds and constraint multipliers corresponding to the
                  optimal solution.  The first num_qp_variables_ entries are for
     the
                  bound multipliers, and the remaining num_qp_constraints_
     entries
                  contain the constraint multipliers. */
  std::shared_ptr<Vector> y_qp_;
  std::shared_ptr<SparseHbMatrix>
      H_; /**< the Matrix object stores the QP data H in
             * Harwell-Boeing Sparse Matrix format*/
  std::shared_ptr<SparseHbMatrix>
      A_; /**< the Matrix object stores the QP data A in
             * Harwell-Boeing Sparse Matrix format*/

  /** Algorithmic options */
  //@{
  // Maximum number of QP iterations per QP solve
  int qp_solver_max_num_iterations_;
  // Print level for QP solver
  int qp_solver_print_level_;
  // Maximum number of LP iterations per LP solve
  int lp_solver_max_num_iterations_;
  //@}
};
}
#endif
