/** Copyright (C) 2019
 * All Rights Reserved.
 *
 * Authors: Xinyi Luo
 * Date:2019-07
 */
#include "IpOptionsList.hpp"
#include "sqphot/QPhandler.hpp"
#include "sqphot/QOREInterface.hpp"
#include "sqphot/qpOASESInterface.hpp"

using namespace std;

namespace SQPhotstart {

QPhandler::QPhandler(std::shared_ptr<const SqpNlpSizeInfo> nlp_sizes, QPType qptype,
                     Ipopt::SmartPtr<Ipopt::Journalist> jnlst,
                     Ipopt::SmartPtr<const Ipopt::OptionsList> options)
 : nlp_sizes_(nlp_sizes)
 , jnlst_(jnlst)
{
  // Determine which solver should be used
  int enum_int;
  options->GetEnumValue("qp_solver_choice", enum_int, "");
  qp_solver_choice_ = Solver(enum_int);

  int num_vars = nlp_sizes->get_num_variables();
  int num_cons = nlp_sizes->get_num_constraints();
#ifdef NEW_FORMULATION
  nConstr_QP_ = num_cons + num_vars;
  nVar_QP_ = num_vars * 3 + 2 * num_cons;
  I_info_A_.length = 5;
  I_info_A_.irow = new int[5];
  I_info_A_.jcol = new int[5];
  I_info_A_.size = new int[5];
  I_info_A_.value = new double[5];
  I_info_A_.irow[0] = I_info_A_.irow[1] = 1;
  I_info_A_.irow[2] = I_info_A_.irow[3] = I_info_A_.irow[4] = nlp_info.nCon + 1;
  I_info_A_.jcol[0] = num_vars + 1;
  I_info_A_.jcol[1] = num_vars + num_cons + 1;
  I_info_A_.jcol[2] = 1;
  I_info_A_.jcol[3] = num_vars + num_cons * 2 + 1;
  I_info_A_.jcol[4] = num_vars * 2 + num_cons * 2 + 1;
  I_info_A_.size[0] = I_info_A_.size[1] = num_cons;
  I_info_A_.size[2] = I_info_A_.size[3] = I_info_A_.size[4] = num_vars;
  I_info_A_.value[0] = I_info_A_.value[2] = I_info_A_.value[3] = 1.0;
  I_info_A_.value[1] = I_info_A_.value[4] = -1.0;

#else
  nConstr_QP_ = num_cons;
  nVar_QP_ = num_vars + 2 * num_cons;
  I_info_A_.length = 2;
  I_info_A_.irow = new int[2];
  I_info_A_.jcol = new int[2];
  I_info_A_.size = new int[2];
  I_info_A_.value = new double[2];
  I_info_A_.irow[0] = I_info_A_.irow[1] = 1;
  I_info_A_.jcol[0] = num_vars + 1;
  I_info_A_.jcol[1] = num_vars + num_cons + 1;
  I_info_A_.size[0] = I_info_A_.size[1] = num_cons;
  I_info_A_.value[0] = 1.0;
  I_info_A_.value[1] = -1.0;
#endif

  W_b_ = new ActiveType[nVar_QP_];
  W_c_ = new ActiveType[nConstr_QP_];

  switch (qp_solver_choice_) {
    case QPOASES:
      solverInterface_ =
          make_shared<qpOASESInterface>(nlp_sizes, qptype, options, jnlst);
      break;
    case QORE:
      solverInterface_ =
          make_shared<QOREInterface>(nlp_sizes, qptype, options, jnlst);
      break;
    case GUROBI:
#ifdef USE_GUROBI
      isolverInterface_ =
          make_shared<GurobiInterface>(nlp_sizes, qptype, options, jnlst);
#endif
      break;
    case CPLEX:
#ifdef USE_CPLEX
      solverInterface_ =
          make_shared<CplexInterface>(nlp_sizes, qptype, options, jnlst);
#endif
      break;
  }

#ifdef DEBUG
#ifdef COMPARE_QP_SOLVER
  qpOASESInterface_ = make_shared<qpOASESInterface>(nlp_sizes, qptype, options);
  QOREInterface_ = make_shared<QOREInterface>(nlp_sizes, qptype, options, jnlst);
  W_b_qpOASES_ = new ActiveType[nVar_QP_];
  W_c_qpOASES_ = new ActiveType[nConstr_QP_];
  W_b_qore_ = new ActiveType[nVar_QP_];
  W_c_qore_ = new ActiveType[nConstr_QP_];
#endif
#endif
}

/**
 *Default destructor
 */
QPhandler::~QPhandler()
{
  delete[] W_b_;
  W_b_ = NULL;
  delete[] W_c_;
  W_c_ = NULL;
  delete[] I_info_A_.irow;
  I_info_A_.irow = NULL;
  delete[] I_info_A_.jcol;
  I_info_A_.jcol = NULL;
  delete[] I_info_A_.size;
  I_info_A_.size = NULL;
  delete[] I_info_A_.value;
  I_info_A_.value = NULL;
#ifdef DEBUG
#ifdef COMPARE_QP_SOLVER
  delete[] W_b_qpOASES_;
  delete[] W_c_qpOASES_;
  delete[] W_b_qore_;
  delete[] W_c_qore_;
#endif
#endif
}

/**
 * Setup the bounds for the QP subproblems according to the information from
 * current
 * iterate. We have
 * 	c_l -c_k <=J_p+ u-v<=c_u-c_k
 * The bound is formulated as
 *   max(-delta, x_l-x_k)<=p<=min(delta,x_u-x_k)
 * and  u,v>=0
 * @param delta      trust region radius
 * @param x_k        current iterate point
 * @param c_k        current constraint value evaluated at x_k
 * @param x_l        the lower bounds for variables
 * @param x_u        the upper bounds for variables
 * @param c_l        the lower bounds for constraints
 * @param c_u        the upper bounds for constraints
 */

void QPhandler::set_bounds(double delta, shared_ptr<const Vector> x_l,
                           shared_ptr<const Vector> x_u,
                           shared_ptr<const Vector> x_k,
                           shared_ptr<const Vector> c_l,
                           shared_ptr<const Vector> c_u,
                           shared_ptr<const Vector> c_k)
{

#ifdef DEBUG
#ifdef COMPARE_QP_SOLVER
  set_bounds_debug(delta, x_l, x_u, x_k, c_l, c_u, c_k);
#endif
#endif

  /*-------------------------------------------------------------*/
  /* Set lbA, ubA as well as lb and ub as qpOASES differentiates */
  /*the bound constraints from the linear constraints            */
  /*-------------------------------------------------------------*/
  int num_vars = nlp_sizes_->get_num_variables();
  int num_cons = nlp_sizes_->get_num_constraints();
  if (qp_solver_choice_ != QORE) {
#ifndef NEW_FORMULATION
    for (int i = 0; i < num_cons; i++) {
      solverInterface_->set_lbA(i,
                                c_l->get_value(i) - c_k->get_value(i)); // must
      // place before set_ubA
      solverInterface_->set_ubA(i, c_u->get_value(i) - c_k->get_value(i));
    }
    for (int i = 0; i < num_vars; i++) {
      solverInterface_->set_lb(
          i, std::max(x_l->get_value(i) - x_k->get_value(i), -delta));
      solverInterface_->set_ub(
          i, std::min(x_u->get_value(i) - x_k->get_value(i), delta));
    }
    /**
     * only set the upper bound for the last half to be infinity(those are slack
     * variables).
     * The lower bounds are initialized as 0
     */
    for (int i = 0; i < num_cons * 2; i++)
      solverInterface_->set_ub(num_vars + i, INF);
#else
    for (int i = 0; i < num_cons; i++) {
      solverInterface_->set_lbA(i, c_l->value(i) - c_k->value(i)); // must
      // place before set_ubA
      solverInterface_->set_ubA(i, c_u->value(i) - c_k->value(i));
    }
    for (int i = 0; i < num_vars; i++) {
      solverInterface_->set_lbA(num_cons + i,
                                x_l->value(i) - x_k->value(i)); // must
      // place before set_ubA
      solverInterface_->set_ubA(num_cons + i,
                                x_u->value(i) - x_k->value(i));
    }

    for (int i = 0; i < num_vars; i++) {
      solverInterface_->set_lb(i, -delta);
      solverInterface_->set_ub(i, delta);
    }
    for (int i = num_vars; i < nVar_QP_; i++)
      solverInterface_->set_ub(i, INF);
#endif
  }
  /*-------------------------------------------------------------*/
  /* Only set lb and ub, where lb = [lbx;lbA]; and ub=[ubx; ubA] */
  /*-------------------------------------------------------------*/
  else {
#ifndef NEW_FORMULATION
    for (int i = 0; i < num_vars; i++) {
      solverInterface_->set_lb(
          i, std::max(x_l->get_value(i) - x_k->get_value(i), -delta));
      solverInterface_->set_ub(
          i, std::min(x_u->get_value(i) - x_k->get_value(i), delta));
    }

    for (int i = 0; i < num_cons * 2; i++) {
      solverInterface_->set_lb(num_vars + i, 0.);
      solverInterface_->set_ub(num_vars + i, INF);
    }
#endif

    for (int i = 0; i < num_cons; i++) {
      solverInterface_->set_lb(nVar_QP_ + i,
                               c_l->get_value(i) - c_k->get_value(i));
      solverInterface_->set_ub(nVar_QP_ + i,
                               c_u->get_value(i) - c_k->get_value(i));
    }

#ifdef NEW_FORMULATION
    for (int i = 0; i < num_vars; i++) {
      solverInterface_->set_lb(i, -delta);
      solverInterface_->set_ub(i, delta);
      solverInterface_->set_lb(nVar_QP_ + num_cons + i,
                               x_l->value(i) - x_k->value(i));
      solverInterface_->set_ub(nVar_QP_ + num_cons + i,
                               x_u->value(i) - x_k->value(i));
    }
    for (int i = num_vars; i < nVar_QP_; i++)
      solverInterface_->set_ub(i, INF);
// DEBUG
//        solverInterface_->getLb()->print("lb");
//        solverInterface_->getUb()->print("ub");
#endif
  }
}

/**
 * This function sets up the object vector g of the QP problem
 * The (2*nCon+nVar) vector g_^T in QP problem will be the same as
 * [grad_f^T, rho* e^T], where the unit vector is of length (2*nCon).
 * @param grad      Gradient vector from nlp class
 * @param rho       Penalty Parameter
 */

void QPhandler::set_g(shared_ptr<const Vector> grad, double rho)
{
  int num_vars = nlp_sizes_->get_num_variables();
#ifdef DEBUG
#ifdef COMPARE_QP_SOLVER
  for (int i = 0; i < nVar_QP_; i++)
    if (i < num_vars) {
      qpOASESInterface_->set_g(i, grad->value(i));
      QOREInterface_->set_g(i, grad->value(i));
    } else {
      qpOASESInterface_->set_g(i, rho);
      QOREInterface_->set_g(i, rho);
    }
#endif
#endif
  for (int i = 0; i < nVar_QP_; i++)
    if (i < num_vars)
      solverInterface_->set_g(i, grad->get_value(i));

    else
      solverInterface_->set_g(i, rho);

  // DEBUG
  //    solverInterface_->getG()->print("G");
}

/**
 * @brief Set up the H for the first time in the QP problem.
 * It will be concatenated as [H_k 0]
 *          		         [0   0]
 * where H_k is the Lagragian hessian evaluated at x_k and lambda_k.
 *
 * This method should only be called for once.
 *
 * @param hessian the Lagragian hessian evaluated at x_k and lambda_k from nlp
 * readers.
 */
void QPhandler::set_H(shared_ptr<const SpTripletMat> hessian)
{
#ifdef DEBUG
#ifdef COMPARE_QP_SOLVER
  qpOASESInterface_->set_H_values(hessian);
  QOREInterface_->set_H_values(hessian);
#endif
#endif
  solverInterface_->set_H(hessian);
}

/**
 * @brief This function sets up the matrix A in the QP subproblem
 * The matrix A in QP problem will be concatenate as [J I -I]
 * @param jacobian  the Matrix object for Jacobian from c(x)
 */
void QPhandler::set_A(shared_ptr<const SpTripletMat> jacobian)
{
#ifdef DEBUG
#ifdef COMPARE_QP_SOLVER
  qpOASESInterface_->set_A_values(jacobian, I_info_A_);
  QOREInterface_->set_A_values(jacobian, I_info_A_);
#endif
#endif
  solverInterface_->set_A(jacobian, I_info_A_);
}

/**
 * @brief This function updates the constraint if there is any changes to
 * the iterates
 */
void QPhandler::update_bounds(double delta, shared_ptr<const Vector> x_l,
                              shared_ptr<const Vector> x_u,
                              shared_ptr<const Vector> x_k,
                              shared_ptr<const Vector> c_l,
                              shared_ptr<const Vector> c_u,
                              shared_ptr<const Vector> c_k)
{
#ifdef DEBUG
#ifdef COMPARE_QP_SOLVER
  set_bounds_debug(delta, x_l, x_u, x_k, c_l, c_u, c_k);
#endif
#endif
  int num_vars = nlp_sizes_->get_num_variables();
  int num_cons = nlp_sizes_->get_num_constraints();
#ifndef NEW_FORMULATION
  if (qp_solver_choice_ != QORE) {
    if (qp_solver_choice_ == GUROBI || qp_solver_choice_ == CPLEX)
      solverInterface_->reset_constraints();

    for (int i = 0; i < num_cons; i++) {
      solverInterface_->set_lbA(i, c_l->get_value(i) - c_k->get_value(i));
    }

    for (int i = 0; i < num_vars; i++) {
      solverInterface_->set_lb(
          i, std::max(x_l->get_value(i) - x_k->get_value(i), -delta));
      solverInterface_->set_ub(
          i, std::min(x_u->get_value(i) - x_k->get_value(i), delta));
    }
  } else {
    for (int i = 0; i < num_vars; i++) {
      solverInterface_->set_lb(
          i, std::max(x_l->get_value(i) - x_k->get_value(i), -delta));
      solverInterface_->set_ub(
          i, std::min(x_u->get_value(i) - x_k->get_value(i), delta));
    }
    for (int i = 0; i < num_cons; i++) {
      solverInterface_->set_lb(num_vars + 2 * num_cons + i,
                               c_l->get_value(i) - c_k->get_value(i));
      solverInterface_->set_ub(num_vars + 2 * num_cons + i,
                               c_u->get_value(i) - c_k->get_value(i));
    }
  }
#else
  if (QPsolverChoice_ == QORE) {
    for (int i = 0; i < num_vars; i++) {
      solverInterface_->set_lb(i, -delta);
      solverInterface_->set_ub(i, delta);
    }
    for (int i = 0; i < num_cons; i++) {
      solverInterface_->set_lb(nVar_QP_ + i, c_l->value(i) - c_k->value(i));
      solverInterface_->set_ub(nVar_QP_ + i, c_u->value(i) - c_k->value(i));
    }

    for (int i = 0; i < num_vars; i++) {
      solverInterface_->set_lb(nVar_QP_ + num_cons + i,
                               x_l->value(i) - x_k->value(i));
      solverInterface_->set_ub(nVar_QP_ + num_cons + i,
                               x_u->value(i) - x_k->value(i));
    }
  } else {
    for (int i = 0; i < num_vars; i++) {
      solverInterface_->set_lb(i, -delta);
      solverInterface_->set_ub(i, delta);
    }

    for (int i = 0; i < num_cons; i++) {
      solverInterface_->set_lbA(i, c_l->value(i) - c_k->value(i));
      solverInterface_->set_ubA(i, c_u->value(i) - c_k->value(i));
    }
    for (int i = 0; i < num_vars; i++) {
      solverInterface_->set_lbA(num_cons + i,
                                x_l->value(i) - x_k->value(i));
      solverInterface_->set_ubA(num_cons + i,
                                x_u->value(i) - x_k->value(i));
    }
  }
#endif
}

/**
 * @brief This function updates the vector g in the QP subproblem when there are
 * any
 * change to the values of penalty parameter
 *
 * @param rho               penalty parameter
 * @param nVar              number of variables in NLP
 */

void QPhandler::update_penalty(double rho)
{
  int num_vars = nlp_sizes_->get_num_variables();
#ifdef DEBUG
#ifdef COMPARE_QP_SOLVER
  int num_cons = nlp_sizes_->get_num_constraints();
  for (int i = num_vars; i < num_vars + num_cons * 2; i++) {
    qpOASESInterface_->set_g(i, rho);
    QOREInterface_->set_g(i, rho);
  }
#endif
#endif
  for (int i = num_vars; i < nVar_QP_; i++)
    solverInterface_->set_g(i, rho);
}

/**
 * @brief This function updates the vector g in the QP subproblem when there are
 * any
 * change to the values of gradient in NLP
 *
 * @param grad              the gradient vector from NLP
 */
void QPhandler::update_grad(shared_ptr<const Vector> grad)
{
  int num_vars = nlp_sizes_->get_num_variables();
#ifdef DEBUG
#ifdef COMPARE_QP_SOLVER

  for (int i = 0; i < num_vars; i++) {
    qpOASESInterface_->set_g(i, grad->value(i));
    QOREInterface_->set_g(i, grad->value(i));
  }
#endif
#endif
  for (int i = 0; i < num_vars; i++)
    solverInterface_->set_g(i, grad->get_value(i));
}

/**
 *@brief Solve the QP with objective and constraints defined by its class
 *members
 */
void QPhandler::solveQP(shared_ptr<Stats> stats, Ipopt::SmartPtr<const Ipopt::OptionsList> options)
{

//    solverInterface_->getA()->print_full("A");
//    solverInterface_->getH()->print_full("H");
//    solverInterface_->getLb()->print("Lb");
//    solverInterface_->getUb()->print("Ub");
//    solverInterface_->getG()->print("G");

#ifdef DEBUG
#ifdef COMPARE_QP_SOLVER
  QOREInterface_->optimizeQP(stats);
  qpOASESInterface_->optimizeQP(stats);
  bool qpOASES_optimal =
      OptimalityTest(qpOASESInterface_, QPOASES, W_b_qpOASES_, W_c_qpOASES_);
  bool qore_optimal =
      OptimalityTest(QOREInterface_, QORE, W_b_qore_, W_c_qore_);
  if (!qpOASES_optimal || !qore_optimal)
    testQPsolverDifference();

#endif
#endif

  solverInterface_->optimizeQP(stats);

  // manually check if the optimality condition is satisfied
  bool isOptimal =
      test_optimality(solverInterface_, qp_solver_choice_, W_b_, W_c_);
  if (!isOptimal) {
    THROW_EXCEPTION(QP_NOT_OPTIMAL, QP_NOT_OPTIMAL_MSG);
  }
}

double QPhandler::get_objective()
{

  return solverInterface_->get_obj_value();
}

void QPhandler::update_H(shared_ptr<const SpTripletMat> Hessian)
{

#ifdef DEBUG
#ifdef COMPARE_QP_SOLVER
  QOREInterface_->set_H_values(Hessian);
  qpOASESInterface_->set_H_values(Hessian);
#endif
#endif
  solverInterface_->set_H(Hessian);
}

void QPhandler::update_A(shared_ptr<const SpTripletMat> Jacobian)
{

#ifdef DEBUG
#ifdef COMPARE_QP_SOLVER
  QOREInterface_->set_A_values(Jacobian, I_info_A_);
  qpOASESInterface_->set_A_values(Jacobian, I_info_A_);
#endif
#endif
  solverInterface_->set_A(Jacobian, I_info_A_);
}

void QPhandler::update_delta(double delta, shared_ptr<const Vector> x_l,
                             shared_ptr<const Vector> x_u,
                             shared_ptr<const Vector> x_k)
{
  int num_vars = nlp_sizes_->get_num_variables();
#ifdef DEBUG
#ifdef COMPARE_QP_SOLVER
  for (int i = 0; i < num_vars; i++) {
    qpOASESInterface_->set_lb(i,
                              std::max(x_l->value(i) - x_k->value(i), -delta));
    qpOASESInterface_->set_ub(i,
                              std::min(x_u->value(i) - x_k->value(i), delta));
    QOREInterface_->set_lb(i, std::max(x_l->value(i) - x_k->value(i), -delta));
    QOREInterface_->set_ub(i, std::min(x_u->value(i) - x_k->value(i), delta));
  }
#endif
#endif

#ifdef NEW_FORMULATION
  for (int i = 0; i < num_vars; i++) {
    solverInterface_->set_lb(i, -delta);
    solverInterface_->set_ub(i, delta);
  }
#else

  for (int i = 0; i < num_vars; i++) {
    solverInterface_->set_lb(
        i, std::max(x_l->get_value(i) - x_k->get_value(i), -delta));
    solverInterface_->set_ub(
        i, std::min(x_u->get_value(i) - x_k->get_value(i), delta));
  }

#endif
}

void QPhandler::WriteQPData(const string filename)
{

  solverInterface_->WriteQPDataToFile(Ipopt::J_LAST_LEVEL, Ipopt::J_USER1,
                                      filename);
}

Exitflag QPhandler::get_status()
{
  return (solverInterface_->get_status());
}

bool QPhandler::test_optimality(shared_ptr<QPSolverInterface> qpsolverInterface,
                                Solver qpSolver, ActiveType* W_b,
                                ActiveType* W_c)
{
  qpOptimalStatus_ = qpsolverInterface->get_optimality_status();
  return (qpsolverInterface->test_optimality(W_c, W_b));
}

double QPhandler::get_infea_measure_model() const
{
  int num_vars = nlp_sizes_->get_num_variables();
  return solverInterface_->get_optimal_solution()->calc_one_norm(
      num_vars, nVar_QP_ - num_vars);
}

const OptimalityStatus& QPhandler::get_QpOptimalStatus() const
{
  return qpOptimalStatus_;
}

void QPhandler::get_active_set(ActiveType* A_c, ActiveType* A_b,
                               shared_ptr<Vector> x, shared_ptr<Vector> Ax)
{
  // use the class member to get the qp problem information
  auto lb = solverInterface_->getLb();
  auto ub = solverInterface_->getUb();
  if (x == nullptr) {
    x = make_shared<Vector>(nVar_QP_);
    x->copy_vector(get_optimal_solution());
  }
  if (Ax == nullptr) {
    Ax = make_shared<Vector>(nConstr_QP_);
    auto A = solverInterface_->getA();
    A->multiply(x, Ax);
  }

  for (int i = 0; i < nVar_QP_; i++) {
    if (abs(x->get_value(i) - lb->get_value(i)) < sqrt_m_eps) {
      if (abs(ub->get_value(i) - x->get_value(i)) < sqrt_m_eps) {
        A_b[i] = ACTIVE_BOTH_SIDE;
      } else
        A_b[i] = ACTIVE_BELOW;
    } else if (abs(ub->get_value(i) - x->get_value(i)) < sqrt_m_eps)
      A_b[i] = ACTIVE_ABOVE;
    else
      A_b[i] = INACTIVE;
  }
  if (qp_solver_choice_ == QORE) {
    // if no x and Ax are input
    for (int i = 0; i < nConstr_QP_; i++) {
      if (abs(Ax->get_value(i) - lb->get_value(i + nVar_QP_)) < sqrt_m_eps) {
        if (abs(ub->get_value(i + nVar_QP_) - Ax->get_value(i)) < sqrt_m_eps) {
          A_c[i] = ACTIVE_BOTH_SIDE;
        } else
          A_c[i] = ACTIVE_BELOW;
      } else if (abs(ub->get_value(i + nVar_QP_) - Ax->get_value(i)) <
                 sqrt_m_eps)
        A_c[i] = ACTIVE_ABOVE;
      else
        A_c[i] = INACTIVE;
    }
  } else {
    auto lbA = solverInterface_->getUbA();
    auto ubA = solverInterface_->getUbA();
    for (int i = 0; i < nConstr_QP_; i++) {
      if (abs(Ax->get_value(i) - lbA->get_value(i)) < sqrt_m_eps) {
        if (abs(ubA->get_value(i) - Ax->get_value(i)) < sqrt_m_eps) {
          A_c[i] = ACTIVE_BOTH_SIDE;
        } else
          A_c[i] = ACTIVE_BELOW;
      } else if (abs(ubA->get_value(i) - Ax->get_value(i)) < sqrt_m_eps)
        A_c[i] = ACTIVE_ABOVE;
      else
        A_c[i] = INACTIVE;
    }
  }
}

void QPhandler::set_g(double rho)
{
  int num_vars = nlp_sizes_->get_num_variables();
  for (int i = 0; i < num_vars; i++) {
    solverInterface_->set_g(i, rho);
  }
  for (int i = num_vars; i < nVar_QP_; i++) {
    solverInterface_->set_g(i, rho);
  }
}

#ifdef DEBUG
#ifdef COMPARE_QP_SOLVER
void QPhandler::set_bounds_debug(double delta, shared_ptr<const Vector> x_l,
                                 shared_ptr<const Vector> x_u,
                                 shared_ptr<const Vector> x_k,
                                 shared_ptr<const Vector> c_l,
                                 shared_ptr<const Vector> c_u,
                                 shared_ptr<const Vector> c_k)
{

  for (int i = 0; i < num_cons; i++) {
    qpOASESInterface_->set_lbA(i, c_l->value(i) - c_k->value(i));
    qpOASESInterface_->set_ubA(i, c_u->value(i) - c_k->value(i));
  }

  for (int i = 0; i < num_vars; i++) {
    qpOASESInterface_->set_lb(i,
                              std::max(x_l->value(i) - x_k->value(i), -delta));
    qpOASESInterface_->set_ub(i,
                              std::min(x_u->value(i) - x_k->value(i), delta));
  }
  /**
   * only set the upper bound for the last 2*nCon entries (those are slack
   * variables).
   * The lower bounds are initialized as 0
   */
  for (int i = 0; i < num_cons * 2; i++)
    qpOASESInterface_->set_ub(num_vars + i, INF);

  /*-------------------------------------------------------------*/
  /* Only set lb and ub, where lb = [lbx;lbA]; and ub=[ubx; ubA] */
  /*-------------------------------------------------------------*/
  for (int i = 0; i < num_vars; i++) {
    QOREInterface_->set_lb(i, std::max(x_l->value(i) - x_k->value(i), -delta));
    QOREInterface_->set_ub(i, std::min(x_u->value(i) - x_k->value(i), delta));
  }

  for (int i = 0; i < nConstr_QP_ * 2; i++)
    QOREInterface_->set_ub(num_vars + i, INF);

  for (int i = 0; i < nConstr_QP_; i++) {
    QOREInterface_->set_lb(nVar_QP_ + i, c_l->value(i) - c_k->value(i));
    QOREInterface_->set_ub(nVar_QP_ + i, c_u->value(i) - c_k->value(i));
  }
}

bool QPhandler::testQPsolverDifference()
{
  shared_ptr<Vector> qpOASESsol = make_shared<Vector>(num_vars);
  shared_ptr<Vector> QOREsol = make_shared<Vector>(num_vars);
  shared_ptr<Vector> difference = make_shared<Vector>(num_vars);
  qpOASESsol->print("qpOASESsol", jnlst_);
  QOREsol->print("QOREsol", jnlst_);
  qpOASESsol->copy_vector(qpOASESInterface_->get_optimal_solution());
  QOREsol->copy_vector(QOREInterface_->get_optimal_solution());
  difference->copy_vector(qpOASESsol);
  difference->add_vector(-1., QOREsol);
  double diff_norm = difference->one_norm();
  if (diff_norm > 1.0e-8) {
    printf("difference is %10e\n", diff_norm);
    qpOASESsol->print("qpOASESsol");
    QOREsol->print("QOREsol");
    QOREInterface_->WriteQPDataToFile(jnlst_, J_ALL, J_DBG);
  }
  assert(diff_norm < 1.0e-8);
  return true;
}

#endif
#endif
} // namespace SQPhotstart
