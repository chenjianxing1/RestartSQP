/* Copyright (C) 2019
 * All Rights Reserved.
 *
 * Authors: Xinyi Luo
 * Date:2019-06
 */
#ifndef SQPHOTSTART_IpoptSqpNlp_HPP
#define SQPHOTSTART_IpoptSqpNlp_HPP

#include "IpTNLP.hpp"
#include "restartsqp/SqpTNlp.hpp"

namespace RestartSqp {
/**
 * This is part of SQPhotstart
 *
 * This class is used in the SqpRestart algorithm.  It is used to
 * post an SqpTNlp as an IpoptTNLP so that Ipopt can solve it.
 * At the end, the solution is not given to the SqpTNLP, though,
 * instead it is stored here.
 */
class IpoptSqpTNlp : public Ipopt::TNLP
{

public:
  /** @brief Constructor that an instance of Ipopt's TNLP */
  IpoptSqpTNlp(std::shared_ptr<SqpTNlp> sqp_tnlp);

  /** Destructor*/
  ~IpoptSqpTNlp();

  /**
   *@brief Get problem size information
   */
  bool get_nlp_info(Ipopt::Index& n, Ipopt::Index& m, Ipopt::Index& nnz_jac_g,
                    Ipopt::Index& nnz_h_lag,
                    TNLP::IndexStyleEnum& index_style) override;

  /**
   *@brief get the bounds information from the NLP object
   */
  bool get_bounds_info(Ipopt::Index n, Ipopt::Number* x_l, Ipopt::Number* x_u,
                       Ipopt::Index m, Ipopt::Number* g_l,
                       Ipopt::Number* g_u) override;

  /*
   * @brief Get the starting point from the NLP object.
   * TODO: add options_ to enable user to choose if to use default input or not
   */
  bool get_starting_point(Ipopt::Index n, bool init_x, Ipopt::Number* x,
                          bool init_z, Ipopt::Number* z_L, Ipopt::Number* z_U,
                          Ipopt::Index m, bool init_lambda,
                          Ipopt::Number* lambda) override;

  /**
   *@brief Evaluate the objective value
   */
  bool eval_f(Ipopt::Index n, const Ipopt::Number* x, bool new_x,
              Ipopt::Number& obj_value) override;

  /**
   * @brief Evaluate the constraints at point x
   *
   */
  bool eval_g(Ipopt::Index n, const Ipopt::Number* x, bool new_x,
              Ipopt::Index m, Ipopt::Number* g) override;

  /**
   *@brief Evaluate gradient at point x
   */
  bool eval_grad_f(Ipopt::Index n, const Ipopt::Number* x, bool new_x,
                   Ipopt::Number* grad_f) override;

  /**
   * @brief Get the matrix structure of the Jacobian
   * Always call this before the first time using @Eval_Jacobian
   */
  bool eval_jac_g(Ipopt::Index n, const Ipopt::Number* x, bool new_x,
                  Ipopt::Index m, Ipopt::Index nele_jac, Ipopt::Index* iRow,
                  Ipopt::Index* jCol, Ipopt::Number* values) override;

  /**
   *@brief Evaluate Hessian of Lagragian function at  (x, lambda)
   */
  bool eval_h(Ipopt::Index n, const Ipopt::Number* x, bool new_x,
              Ipopt::Number obj_factor, Ipopt::Index m,
              const Ipopt::Number* lambda, bool new_lambda,
              Ipopt::Index nele_hess, Ipopt::Index* iRow, Ipopt::Index* jCol,
              Ipopt::Number* values) override;

  /**
   * @brief Return the results of the optimization run to the user.
   */
  void finalize_solution(Ipopt::SolverReturn status, Ipopt::Index n,
                         const Ipopt::Number* x, const Ipopt::Number* z_L,
                         const Ipopt::Number* z_U, Ipopt::Index m,
                         const Ipopt::Number* g, const Ipopt::Number* lambda,
                         Ipopt::Number obj_value,
                         const Ipopt::IpoptData* ip_data,
                         Ipopt::IpoptCalculatedQuantities* ip_cq) override;

  /** @name Accessor methods for the Ipopt solution. */
  //@{
  const double* x_sol() const
  {
    return x_sol_;
  }
  const double* z_L_sol() const
  {
    return z_L_sol_;
  }
  const double* z_U_sol() const
  {
    return z_U_sol_;
  }
  const double* lambda_sol() const
  {
    return lambda_sol_;
  }
  const double* g_sol() const
  {
    return g_sol_;
  }
  double obj_value() const
  {
    return obj_value_;
  }
  //@}

private:
  /** @name Hide unused default methods. */
  //@{
  /** Default constructor*/
  IpoptSqpTNlp();

  /** Copy Constructor */
  IpoptSqpTNlp(const IpoptSqpTNlp&);

  /** Overloaded Equals Operator */
  void operator=(const IpoptSqpTNlp&);
  //@}

  /** Ipopt's TNLP object that will be called for all evaluations. */
  std::shared_ptr<SqpTNlp> sqp_tnlp_;

  /** Store the solution. */
  //@{
  /** Ipopt status after most recent solve. */
  Ipopt::SolverReturn ipopt_status_;
  /** Primal solution. */
  double* x_sol_;
  /** Multipliers for lower bounds. */
  double* z_L_sol_;
  /** Multipliers for upper bounds. */
  double* z_U_sol_;
  /** Multipliers for the constraints.  This is with reversed sign from Ipopt,
   * so consistent with SqpNLP. */
  double* lambda_sol_;
  /** Constraint values at final iterate. */
  double* g_sol_;
  /** Objective function at final iterate. */
  double obj_value_;
  //@}
};
} // namespace RestartSqp

#endif // SQPHOTSTART_IpoptSqpNlp_HPP
