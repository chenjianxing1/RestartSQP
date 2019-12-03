/* Copyright (C) 2019
 * All Rights Reserved.
 *
 * Authors: Xinyi Luo
 * Date:2019-09-06
 */
#ifndef _SQPHOTSTART_GUROBI_INTERFACE_
#define _SQPHOTSTART_GUROBI_INTERFACE_

#include "sqphot/QPsolverInterface.hpp"

#ifdef USE_GUROBI
#include <gurobi_c++.h>
#endif

namespace SQPhotstart {
DECLARE_STD_EXCEPTION(GRB_SOLVER_FAILS);
class GurobiInterface : public QPSolverInterface
{

public:
  GurobiInterface(NLPInfo nlp_info, QPType qptype,
                  std::shared_ptr<const Options> options,
                  Ipopt::SmartPtr<Ipopt::Journalist> jnlst);

  /**@name Getters for private members*/
  //@{
  const std::shared_ptr<Vector>& getLb() const override;

  const std::shared_ptr<Vector>& getUb() const override;

  const std::shared_ptr<Vector>& getLbA() const override;

  const std::shared_ptr<Vector>& getUbA() const override;

  const std::shared_ptr<Vector>& getG() const override;

  std::shared_ptr<const SpHbMat> getH() const override;

  std::shared_ptr<const SpHbMat> getA() const override;
  //@}

  /** Default destructor*/
  ~GurobiInterface();

  /**
   * @brief Solve a regular QP with given data and options.
   */
  void optimizeQP(std::shared_ptr<Stats> stats) override;

  /**
   * @brief Solve a regular LP with given data and options
   *
   */

  void optimizeLP(std::shared_ptr<Stats> stats) override;

  bool test_optimality(ActiveType* W_c = NULL, ActiveType* W_b = NULL) override;

  /**-------------------------------------------------------**/
  /**                    Getters                            **/
  /**-------------------------------------------------------**/
  /**@name Getters*/
  //@{
  /**
   * @return the pointer to the optimal solution
   *
   */
  std::shared_ptr<const Vector> get_optimal_solution() const override
  {
    return x_qp;
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
  std::shared_ptr<const Vector> get_bounds_multipliers() const override;

  /**
   * @brief get the pointer to the multipliers to the regular constraints.
   */
  std::shared_ptr<const Vector> get_constraints_multipliers() const override;

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
  void get_working_set(ActiveType* W_constr, ActiveType* W_bounds) override;

  Exitflag get_status() override;

  //@}

  /**-------------------------------------------------------**/
  /**                    Setters                            **/
  /**-------------------------------------------------------**/
  /**@name Setters, by location and value*/
  //@{
  void set_lb(int location, double value) override;

  void set_ub(int location, double value) override;

  void set_lbA(int location, double value) override;

  void set_ubA(int location, double value) override;

  void set_g(int location, double value) override;
  //@}

  /**@name Setters for dense vector, by vector value*/
  //@{
  void set_ub(std::shared_ptr<const Vector> rhs) override;

  void set_lb(std::shared_ptr<const Vector> rhs) override;

  void set_lbA(std::shared_ptr<const Vector> rhs) override;

  void set_ubA(std::shared_ptr<const Vector> rhs) override;

  void set_g(std::shared_ptr<const Vector> rhs) override;
  //@}

  /**@name Setters for matrix*/
  //@{
  void set_H(std::shared_ptr<const SpTripletMat> rhs) override;

  void set_A(std::shared_ptr<const SpTripletMat> rhs,
             IdentityInfo I_info) override;
  //@}

  void reset_model();
  /**-------------------------------------------------------**/
  /**                  Data Writer                          **/
  /**-------------------------------------------------------**/

  void WriteQPDataToFile(Ipopt::EJournalLevel level,
                         Ipopt::EJournalCategory category,
                         const std::string filename) override{};

private:
  /** Default constructor*/
  GurobiInterface();

  /** Copy Constructor */
  GurobiInterface(const GurobiInterface&);

  /** Overloaded Equals Operator */
  void operator=(const GurobiInterface&);

  void set_solver_options();

  void reset_constraints() override;
  /**-------------------------------------------------------**/
  /**                  Private Members                      **/
  /**-------------------------------------------------------**/

private:
#ifdef USE_GUROBI
  GRBEnv* grb_env_;
  GRBLinExpr lterm_;
  GRBModel* grb_mod_;
  GRBQuadExpr qobj_;
  GRBVar* grb_vars_;
  vector<GRBConstr> grb_constr;
#endif
  IdentityInfo I_info_;
  Ipopt::SmartPtr<Ipopt::Journalist> jnlst_;
  Exitflag status_;
  QPType qptype_;
  int nConstr_QP_;
  int nVar_QP_;
  std::shared_ptr<Vector> x_qp;
  std::shared_ptr<Vector> y_qp;
  std::shared_ptr<const Options> options_;
  std::shared_ptr<const SpTripletMat> A_;
};
}

#endif
