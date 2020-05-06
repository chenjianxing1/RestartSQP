/* Copyright (C) 2019
 * All Rights Reserved.
 *
 * Authors: Xinyi Luo      2019-07
 */
#ifndef SQPHOTSTART_QPSOLVER_INTERFACE_HPP
#define SQPHOTSTART_QPSOLVER_INTERFACE_HPP

#include "restartsqp/SparseHbMatrix.hpp"
#include "restartsqp/SqpTNlp.hpp"
#include "restartsqp/Statistics.hpp"
#include "restartsqp/Vector.hpp"
#include <memory>

namespace RestartSqp {

enum QpSolverExitStatus
{
  QPEXIT_UNINITIALIZED = -1,
  QPEXIT_OPTIMAL = 0,
  QPEXIT_UNBOUNDED = 1,
  QPEXIT_INFEASIBLE = 2,
  QPEXIT_ITERLIMIT = 3,
  QPEXIT_UNKNOWN_STATUS = -96,
  QPEXIT_FAILED = -97,
  QPEXIT_NOT_SOLVED = -98,
  QPEXIT_INTERNAL_ERROR = -99
};

// TODO: Make class and put KKT error computation as method
typedef struct
{
  double primal_infeasibility = 0.;
  double dual_infeasibility = 0.;
  double complementarity_violation = 0.;
  double working_set_error = 0.;
  double worst_violation = 0.;
  bool Second_order_opt = false;
} KktError;

/** Function for computing the different parts of the KKT error, given data and
 *  values.  This is used for the NLP and QP error.
 *
 *  For the computation fo the QP error, we need to give the Hessian matrix,
 * and for the NLP error, we need to give the constraint violation.  The method
 * figures out which to use by is_nlp flag. */
KktError
calc_kkt_error_(Ipopt::SmartPtr<Ipopt::Journalist> jnlst,
                Ipopt::EJournalLevel level, bool is_nlp,
                std::shared_ptr<const Vector> lower_variable_bounds,
                std::shared_ptr<const Vector> upper_variable_bounds,
                std::shared_ptr<const Vector> lower_constraint_bounds,
                std::shared_ptr<const Vector> upper_constraint_bounds,
                std::shared_ptr<const Vector> linear_objective_coefficients,
                std::shared_ptr<const Vector> constraint_values,
                std::shared_ptr<const Matrix> jacobian,
                std::shared_ptr<const Matrix> hessian,
                std::shared_ptr<const Vector> primal_solution,
                std::shared_ptr<const Vector> bound_multipliers,
                std::shared_ptr<const Vector> constraint_multipliers,
                const ActivityStatus* bounds_working_set,
                const ActivityStatus* constraints_working_set);

/**
 *
 * Base class for QP solvers to solving problems of the form
 *
 *  minimize 1/2 x^T H x + g^T x
 *  subject  lb_A<=Ax<=ub_A
 *              lb<=x<=ub
 *
 *
 */
class QpSolverInterface
{
public:
  /** Constructor.  Needs to be given the dimensions of the problem. */
  QpSolverInterface(QPType qp_type, int num_qp_variables,
                    int num_qp_constraints,
                    Ipopt::SmartPtr<Ipopt::Journalist> jnlst);

  /** Destructor */
  virtual ~QpSolverInterface();

  /** Return number of variables. */
  int get_num_variables() const
  {
    return num_qp_variables_;
  }
  /** Return number of constraints. */
  int get_num_constraints() const
  {
    return num_qp_constraints_;
  }

  /** Const getter functions of QP data. */
  //@{
  std::shared_ptr<const Vector> get_lower_variable_bounds() const
  {
    return lower_variable_bounds_;
  }

  std::shared_ptr<const Vector> get_upper_variable_bounds() const
  {
    return upper_variable_bounds_;
  }

  double get_lower_variable_bound(int index) const
  {
    return lower_variable_bounds_->get_value(index);
  }

  double get_upper_variable_bound(int index) const
  {
    return upper_variable_bounds_->get_value(index);
  }

  std::shared_ptr<const Vector> get_lower_constraint_bounds() const
  {
    return lower_constraint_bounds_;
  }

  std::shared_ptr<const Vector> get_upper_constraint_bounds() const
  {
    return upper_constraint_bounds_;
  }

  std::shared_ptr<const Vector> get_linear_objective_coefficients() const
  {
    return linear_objective_coefficients_;
  }
  //@}

  /** @name Non-const getter functions of QP data.  These need to be used to set
   * the data. */
  //@{
  std::shared_ptr<Vector> get_lower_variable_bounds_nonconst()
  {
    return lower_variable_bounds_;
  }

  std::shared_ptr<Vector> get_upper_variable_bounds_nonconst()
  {
    return upper_variable_bounds_;
  }

  std::shared_ptr<Vector> get_lower_constraint_bounds_nonconst()
  {
    return lower_constraint_bounds_;
  }

  std::shared_ptr<Vector> get_upper_constraint_bounds_nonconst()
  {
    return upper_constraint_bounds_;
  }

  std::shared_ptr<Vector> get_linear_objective_coefficients_nonconst()
  {
    return linear_objective_coefficients_;
  }
  //@}

  /** @name Methods for setting elements in vector data. */
  //@{
  void set_lower_variable_bound(int index, double value)
  {
    lower_variable_bounds_->set_value(index, value);
  }
  void set_upper_variable_bound(int index, double value)
  {
    upper_variable_bounds_->set_value(index, value);
  }
  void set_lower_constraint_bound(int index, double value)
  {
    lower_constraint_bounds_->set_value(index, value);
  }
  void set_upper_constraint_bound(int index, double value)
  {
    upper_constraint_bounds_->set_value(index, value);
  }
  void set_linear_objective_coefficient(int index, double value)
  {
    linear_objective_coefficients_->set_value(index, value);
  }
  //@}

  /**@name Methods for giving matrices to the QP solver. */
  //@{
  /** Set the Hessian matrix */
  virtual void set_objective_hessian(
      std::shared_ptr<const SparseTripletMatrix> triplet_matrix) = 0;

  /** Set the Jacobian matrix, including multiples of the identify. */
  virtual void set_constraint_jacobian(
      std::shared_ptr<const SparseTripletMatrix> triplet_matrix,
      IdentityMatrixPositions& identity_matrix_positions) = 0;
  //@}

  /**
   * @brief Solve the QP or LP with the current data.
   *
   *  This returns the QP solver_exit status of the solve.  If the solve was
   * successful,
   *  this also sets the primal and dual solution, optimal objective value.
   */
  QpSolverExitStatus optimize(std::shared_ptr<Statistics> stats);

  /** Exit flag from the most recent solve. */
  QpSolverExitStatus get_solver_status() const
  {
    return solver_status_;
  }

  /** Methods for accessing the optimal solution.  Those are only valid if
   * solver status is QPEXIT_OPTIMAL. */
  //@{
  /** Most recent primal solution */
  std::shared_ptr<const Vector> get_primal_solution() const
  {
    assert(solver_status_ == QPEXIT_OPTIMAL);
    return primal_solution_;
  }
  /** Most recent constraint multipliers */
  std::shared_ptr<const Vector> get_constraint_multipliers() const
  {
    assert(solver_status_ == QPEXIT_OPTIMAL);
    return constraint_multipliers_;
  }
  /** Most recent variable multipliers */
  std::shared_ptr<const Vector> get_bounds_multipliers() const
  {
    assert(solver_status_ == QPEXIT_OPTIMAL);
    return bound_multipliers_;
  }
  /** Most recent optimal objective value (as reported by the QP solver) */
  double get_optimal_objective_value() const
  {
    assert(solver_status_ == QPEXIT_OPTIMAL);
    std::shared_ptr<Vector> Hx = std::make_shared<Vector>(num_qp_variables_);

    Hx->set_to_zero();
    if (hessian_) {
      hessian_->multiply(primal_solution_, Hx);
    }

    double optimal_objective_value =
        Hx->calc_inner_product(primal_solution_) * 0.5 +
        linear_objective_coefficients_->calc_inner_product(primal_solution_);
    return optimal_objective_value;
  }
  //@}

  /** Set the initial working set for the bounds and the constraints.  This is
   * only allowed if no QP
   * has been solved yet. */
  void set_initial_working_sets(const ActivityStatus* bounds_working_set,
                                const ActivityStatus* constraints_working_set);

  /**
   *  Working set for the bounds for most recent solve.
   */
  const ActivityStatus* get_bounds_working_set() const
  {
    assert(solver_status_ == QPEXIT_OPTIMAL);
    assert(working_set_up_to_date_);
    return bounds_working_set_;
  }

  /**
   *  Working set for the constraints for most recent solve.
   */
  const ActivityStatus* get_constraints_working_set() const
  {
    assert(solver_status_ == QPEXIT_OPTIMAL);
    assert(working_set_up_to_date_);
    return constraints_working_set_;
  }

  /** Return the number of QP solver iterations from the most recent solve. */
  int get_num_qp_iterations() const
  {
    assert(solver_status_ == QPEXIT_OPTIMAL);
    return qp_solver_iterations_;
  }

  /** This method checks the KKT conditions for the most solution of the most
   *  recent solve.  It returns the maximum violation (KKT error). */
  KktError calc_kkt_error(Ipopt::EJournalLevel level = Ipopt::J_NONE) const;

  /** Write the current QP data to a file. */
  void write_qp_data_to_file(const std::string& filename);

protected:
  /** Implementation of the optimization method. */
  virtual QpSolverExitStatus
  optimize_impl(std::shared_ptr<Statistics> stats = nullptr) = 0;

  /** Method for computing the working set from the most recent solve. */
  virtual void retrieve_working_set_() = 0;

  /** QP type, indicating whether this is for solving QPs or LPs. */
  QPType qp_type_;

  /** Number of variables. */
  int num_qp_variables_;
  /** Number of constraints. */
  int num_qp_constraints_;

  /** QP data */
  //@{
  /** linear part of the objective function. */
  std::shared_ptr<Vector> linear_objective_coefficients_;
  /** lower bounds for the constraints. */
  std::shared_ptr<Vector> lower_constraint_bounds_;
  /** upper bounds for the constraints. */
  std::shared_ptr<Vector> upper_constraint_bounds_;
  /** lower bounds for the variables. */
  std::shared_ptr<Vector> lower_variable_bounds_;
  /** upper bounds for the variables. */
  std::shared_ptr<Vector> upper_variable_bounds_;

  /** Sparse objective Hessian matrix.  Depending on the QP solver, this is in
   * compressed column (qpOASES) or compressd row (QORE) format. */
  std::shared_ptr<SparseHbMatrix> hessian_;

  /** Sparse constraint Jacobian matrix.  Depending on the QP solver, this is in
   * compressed column (qpOASES) or compressd row (QORE) format. */
  std::shared_ptr<SparseHbMatrix> jacobian_;
  //@}

  /** Exit status from most recent solve. */
  QpSolverExitStatus solver_status_;

  /** QP solution. */
  //@{
  /** Most recent primal solution. */
  std::shared_ptr<Vector> primal_solution_;
  /** Most recent constraint multipliers (dual solution). */
  std::shared_ptr<Vector> constraint_multipliers_;
  /** Most recent variable multipliers (dual solution). */
  std::shared_ptr<Vector> bound_multipliers_;

  /** Flag indicating whether the working set information is up-to-date. */
  bool working_set_up_to_date_;
  /** Working set for the bound constraints.  If NULL, no QP has been solved and
   * no initial working set has been provided. */
  ActivityStatus* bounds_working_set_;
  /** Working set for the regular constraints.  If NULL, no QP has been solved
   * and no initial working set has been provided. */
  ActivityStatus* constraints_working_set_;

  /** Number or QP solver iterations. */
  int qp_solver_iterations_;
  //@}

  /** Journalist for output. */
  Ipopt::SmartPtr<Ipopt::Journalist> jnlst_;

  /** Counter for naming the QP output files. */
  int file_counter_;

private:
  /** Default constructor*/
  QpSolverInterface();

  /** Copy Constructor */
  QpSolverInterface(const QpSolverInterface&);

  /** Overloaded Equals Operator */
  void operator=(const QpSolverInterface&);
};

} // SQPHOTSTART
#endif
