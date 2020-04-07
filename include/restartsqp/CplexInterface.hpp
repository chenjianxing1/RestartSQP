/* Copyright (C) 2019
 * All Rights Reserved.
 *
 * Authors: Xinyi Luo
 * Data: 2019-09-09
 */
#ifndef SQPHOTSTART_CPLEX_INTERFACE_HPP
#define SQPHOTSTART_CPLEX_INTERFACE_HPP

#ifdef USE_CPLEX
#include <ilcplex/ilocplex.h>
#endif

#include "restartsqp/QPsolverInterface.hpp"

namespace RestartSqp {
DECLARE_STD_EXCEPTION(CPLEX_SOLVER_FAILS);

class CplexInterface : public QPSolverInterface
{
  /**-------------------------------------------------------**/
  /**                  Public Methods                       **/
  /**-------------------------------------------------------**/
public:
  /**@name Getters for private members*/
  //@{
  const std::shared_ptr<Vector>& get_lower_variable_bounds() const override;

  const std::shared_ptr<Vector>& get_upper_variable_bounds() const override;

  const std::shared_ptr<Vector>& get_lower_constraint_bounds() const override;

  const std::shared_ptr<Vector>& get_upper_constraint_bounds() const override;

  const std::shared_ptr<Vector>&
  get_linear_objective_coefficients() const override;

  std::shared_ptr<const SparseHbMatrix> get_objective_hessian() const override;

  std::shared_ptr<const SparseHbMatrix>
  get_constraint_jacobian() const override;
  //@}

  CplexInterface(NLPInfo nlp_info, QPType qptype,
                 std::shared_ptr<const Options> options,
                 Ipopt::SmartPtr<Ipopt::Journalist> jnlst);

  /** Default destructor*/
  ~CplexInterface() override;

  /**
   * @brief Solve a regular QP with given data and options.
   *
   */
  void optimize_qp(std::shared_ptr<Statistics> stats) override;

  /**
   * @brief Solve a regular LP with given data and options
   */

  void optimize_lp(std::shared_ptr<Statistics> stats) override;

  bool test_optimality(ActivityStatus* W_c = NULL,
                       ActivityStatus* W_b = NULL) override;

  /**-------------------------------------------------------**/
  /**                    Getters                            **/
  /**-------------------------------------------------------**/
  /**@name Getters*/
  //@{
  /**
   * @return the pointer to the optimal solution
   *
   */
  std::shared_ptr<const Vector> get_primal_solution() const override
  {
    return x_qp_;
  }

  /**
   *@brief get the objective value from the QP solvers
   *
   * @return the objective function value of the QP problem
   */
  double get_obj_value() override;

  /**
   * @brief get the pointer to the multipliers to the bounds constraints.
   */
  std::shared_ptr<const Vector> get_bounds_multipliers() const override
  {
    return nullptr;
  }

  /**
   * @brief get the pointer to the multipliers to the regular constraints.
   */
  std::shared_ptr<const Vector> get_constraints_multipliers() const override
  {
    return y_qp_;
  }

  /**
   * @brief copy the working set information
   * @param W_constr a pointer to an array of length (nCon_QP_) which will store
   * the
   * working set for constraints
   * @param W_bounds a pointer to an array of length (nVar_QP_) which will store
   * the
   * working set for bounds
   *
   */
  void get_working_set(ActivityStatus* W_constr,
                       ActivityStatus* W_bounds) override;

  Exitflag get_status() override;

  //@}

  /**-------------------------------------------------------**/
  /**                    Setters                            **/
  /**-------------------------------------------------------**/
  /**@name Setters, by location and value*/
  //@{
  void set_lower_variable_bounds(int location, double value) override;

  void set_upper_variable_bounds(int location, double value) override;

  void set_lower_constraint_bounds(int location, double value) override;

  void set_upper_constraint_bounds(int location, double value) override;

  void set_linear_objective_coefficients(int location, double value) override;
  //@}

  /**@name Setters for dense vector, by vector value*/
  //@{
  void set_upper_variable_bounds(std::shared_ptr<const Vector> rhs) override;

  void set_lower_variable_bounds(std::shared_ptr<const Vector> rhs) override;

  void set_lower_constraint_bounds(std::shared_ptr<const Vector> rhs) override;

  void set_upper_constraint_bounds(std::shared_ptr<const Vector> rhs) override;

  void
  set_linear_objective_coefficients(std::shared_ptr<const Vector> rhs) override;
  //@}

  /**@name Setters for matrix*/
  //@{
  void set_objective_hessian(std::shared_ptr<const SpTripletMat> rhs) override;

  void set_constraint_jacobian(std::shared_ptr<const SpTripletMat> rhs,
                               IdentityMatrixPositions I_info) override;
  //@}

  void reset_constraints() override;

  /**-------------------------------------------------------**/
  /**                  Data Writer                          **/
  /**-------------------------------------------------------**/

  void WriteQPDataToFile(Ipopt::EJournalLevel level,
                         Ipopt::EJournalCategory category,
                         const std::string filename) override
  {
  }

  /**-------------------------------------------------------**/
  /**                  Private Methods                      **/
  /**-------------------------------------------------------**/
private:
  /** Default constructor*/
  CplexInterface();

  /** Copy Constructor */
  CplexInterface(const CplexInterface&);

  /** Overloaded Equals Operator */
  void operator=(const CplexInterface&);

  void set_solver_options();
  /**-------------------------------------------------------**/
  /**                  Private Members                      **/
  /**-------------------------------------------------------**/

private:
  IdentityMatrixPositions I_info_;
  Ipopt::SmartPtr<Ipopt::Journalist> jnlst_;
  Exitflag status_;
  QPType qptype_;
  bool firstQPsolved_;
  double final_obj_;
  int nConstr_QP_;
  int nVar_QP_;
  std::shared_ptr<Vector> lb_;
  std::shared_ptr<Vector> ub_;
  std::shared_ptr<Vector> x_qp_;
  std::shared_ptr<Vector> y_qp_;
  std::shared_ptr<const Options> options_;
  std::shared_ptr<const SpTripletMat> A_;
#ifdef USE_CPLEX
  std::shared_ptr<IloEnv> cplex_env_;
  std::shared_ptr<IloModel> cplex_model_;
  IloNumExpr qobj_; /**< quadratic part of the objecitve*/
  vector<IloNumExpr> constraints_;
  vector<IloNumExpr> lterm_; /**< linear part of the objective */
  vector<IloNumVar> cplex_vars_;
#endif
};

} // SQPHOTSTART
#endif
