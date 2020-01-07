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
#include "sqphot/MessageHandling.hpp"
#include "sqphot/QPsolverInterface.hpp"
#include "sqphot/SqpNlpBase.hpp"

DECLARE_STD_EXCEPTION(INVALID_RETURN_TYPE);

namespace SQPhotstart {

class QOREInterface : public QPSolverInterface
{

  ///////////////////////////////////////////////////////////
  //                      PUBLIC METHODS                   //
  ///////////////////////////////////////////////////////////
public:
  /**Constructor*/

  QOREInterface(std::shared_ptr<const SqpNlpSizeInfo> nlp_sizes, QPType qptype,
                Ipopt::SmartPtr<const Ipopt::OptionsList> options,
                Ipopt::SmartPtr<Ipopt::Journalist> jnlst);

  QOREInterface(std::shared_ptr<SparseHbMatrix> H,
                std::shared_ptr<SparseHbMatrix> A, std::shared_ptr<Vector> g,
                std::shared_ptr<Vector> lb, std::shared_ptr<Vector> ub,
                Ipopt::SmartPtr<const Ipopt::OptionsList> options = nullptr);

  /** Default destructor*/
  ~QOREInterface();

  /**
   * @brief Solve a regular QP with given data and options.
   */

  QpSolverExitStatus
  optimize_qp(std::shared_ptr<Statistics> stats = nullptr) override;

  /**
   * @brief Solve a regular LP with given data and options
   */

  QpSolverExitStatus
  optimize_lp(std::shared_ptr<Statistics> stats = nullptr) override;

  bool test_optimality(ActivityStatus* W_c = NULL,
                       ActivityStatus* W_b = NULL) override;

  /**@name Getters */
  //@{
  /**
   *@brief get the objective value from the QP solvers
   *
   * @return the objective function value of the QP problem
   */
  double get_obj_value() override;

  /**
   * @brief copy the optimal solution of the QP to the input pointer
   *
   * @param x_optimal a pointer to an empty array with allocated memory euqal to
   * sizeof(double)*number_variables
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

  Exitflag get_status() override;

  void get_working_set(ActivityStatus* W_constr,
                       ActivityStatus* W_bounds) override;

  OptimalityStatus get_optimality_status() override
  {
    return qpOptimalStatus_;
  }
  //@}
  //

  /**@name Getters */
  //@{
  std::shared_ptr<const Vector> get_linear_objective_coefficients() const override
  {
    return g_;
  }

  std::shared_ptr<const Vector> get_lower_variable_bounds() const override
  {
    return lb_;
  }

  std::shared_ptr<const Vector> get_upper_variable_bounds() const override
  {
    return ub_;
  }

  std::shared_ptr<const Vector> get_lower_constraint_bounds() const override
  {
    THROW_EXCEPTION(INVALID_RETURN_TYPE, INVALID_RETURN_TYPE_MSG);
  }

  std::shared_ptr<const Vector> get_upper_constraint_bounds() const override
  {
    THROW_EXCEPTION(INVALID_RETURN_TYPE, INVALID_RETURN_TYPE_MSG);
  }

  std::shared_ptr<const SparseHbMatrix> get_objective_hessian() const override
  {
    return H_;
  };

  std::shared_ptr<const SparseHbMatrix> get_constraint_jacobian() const override
  {
    return A_;
  };
  //@}
  /** @name Setters */
  //@{
  void set_linear_objective_coefficients(int location, double value) override
  {
    value = value < INF ? value : INF;
    g_->set_value(location, value);
  };

  void set_lower_variable_bounds(int location, double value) override
  {
    value = value > -INF ? value : -INF;
    lb_->set_value(location, value);
  };

  void set_upper_variable_bounds(int location, double value) override
  {
    value = value < INF ? value : INF;
    ub_->set_value(location, value);
  };

  void
  set_constraint_jacobian(std::shared_ptr<const SpTripletMat> rhs,
               IdentityMatrixPositions& identity_matrix_positions) override;

  void set_objective_hessian(std::shared_ptr<const SpTripletMat> rhs) override;

  //@}
  void WriteQPDataToFile(Ipopt::EJournalLevel level,
                         Ipopt::EJournalCategory category,
                         const std::string filename) override;

  //@{
  void set_linear_objective_coefficients(std::shared_ptr<const Vector> rhs) override
  {
    g_->copy_vector(rhs);
  };

  void set_lower_variable_bounds(std::shared_ptr<const Vector> rhs) override
  {
    lb_->copy_vector(rhs);
  };

  void set_upper_variable_bounds(std::shared_ptr<const Vector> rhs) override
  {
    ub_->copy_vector(rhs);
  };

  void set_lower_constraint_bounds(int location, double value) override{};
  void set_lower_constraint_bounds(std::shared_ptr<const Vector> rhs) override{};
  void set_upper_constraint_bounds(int location, double value) override{};
  void set_upper_constraint_bounds(std::shared_ptr<const Vector> rhs) override{};

  void reset_constraints() override
  {
    lb_->set_to_zero();
    ub_->set_to_zero();
  };
  //@}

  ///////////////////////////////////////////////////////////
  //                      PRIVATE METHODS                  //
  ///////////////////////////////////////////////////////////
private:
  /** Default Constructor*/
  QOREInterface();

  /** Copy Constructor */
  QOREInterface(const QOREInterface&);

  /** Overloaded Equals Operator */
  void operator=(const QOREInterface&);

  void get_option_values_(Ipopt::SmartPtr<const Ipopt::OptionsList> options);

  /**
   * @brief set options of QP solver based on the user-defined values
   */
  void set_qp_solver_options_();

  /**
   * @brief Handle errors based on current status
   */
  void handle_error(QPType qptype, std::shared_ptr<Statistics> stats = nullptr);
  /**
   * @brief Allocate memory for the class members
   * @param nlp_index_info  the struct that stores simple nlp dimension info
   * @param qptype is the problem to be solved QP or LP?
   */
  void allocate_memory(std::shared_ptr<const SqpNlpSizeInfo> nlp_sizes,
                       QPType qptype);

  ///////////////////////////////////////////////////////////
  //                      PRIVATE MEMBERS                  //
  ///////////////////////////////////////////////////////////
private:
  int status_;
  QoreProblem* solver_;
  bool firstQPsolved_ = false;
  int nConstr_QP_;
  int nVar_QP_;
  bool matrix_change_flag_ = false;
  OptimalityStatus qpOptimalStatus_;
  std::shared_ptr<SparseHbMatrix> A_;
  std::shared_ptr<SparseHbMatrix> H_;
  std::shared_ptr<Vector> g_;
  std::shared_ptr<Vector> lb_;
  std::shared_ptr<Vector> ub_;
  std::shared_ptr<Vector> x_qp_;
  /** the bounds and constraint multipliers corresponding to the
      optimal solution.  The first nVar_QP_ entries are for the bound
      multipliersr, and the remaining nConstr_QP_ entries contain the
      constraint multipliers. */
  std::shared_ptr<Vector> y_qp_;

  int qpiter_[1];
  Ipopt::SmartPtr<Ipopt::Journalist> jnlst_;
  int rv_; // temporarily placed here, for recording the return value from the
           // solver
  int* working_set_;

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
