/* Copyright (C) 2019
* All Rights Reserved.
*
* Authors: Xinyi Luo
* Date:2019-06
*/
#ifndef SQPHOTSTART_TYPES_HPP
#define SQPHOTSTART_TYPES_HPP

#include "IpException.hpp"

/** Declaration of all non-built in types (except for classes) */

namespace RestartSqp {

enum PrintLevel
{
  NO_PRINT = 0,
  PRINT_TO_FILE = 1,
  PRINT_TO_CONSOLE = 2,
  PRINT_TO_BOTH = 3
};

// enum QPReturnType {
//    UNSOLVED = -99,
//    QP_OPTIMAL = 0,
//    QP_NOTINITIALISED = -1,
//    QP_EXCEED_MAX_ITER = 1,
//    QP_INFEASIBLE = -2,
//    QP_UNBOUNDED = -3,
//    QP_PREPARINGAUXILIARYQP = -6,
//    QP_AUXILIARYQPSOLVED = -7,
//    QP_PERFORMINGHOMOTOPY = -8,
//    QP_HOMOTOPYQPSOLVED = -9,
//    QP_UNKNOWN_ERROR = 99
//};

enum QPType
{
  QP_TYPE_LP = 1, /** solving a linear program*/
  QP_TYPE_QP = 2  /** solving a regular quadratic subproblem **/
};

enum SqpSolverExitStatus
{
  OPTIMAL = 0,
  INVALID_NLP = -1,
  CONVERGE_TO_NONOPTIMAL = 1,
  EXCEED_MAX_ITERATIONS = -2, // exceeds the maximum number of iteration
  PRED_REDUCTION_NEGATIVE = -3,
  TRUST_REGION_TOO_SMALL = -4,
  EXCEED_MAX_CPU_TIME = -6,
  EXCEED_MAX_WALLCLOCK_TIME = -7,
  PENALTY_TOO_LARGE = -8,
  QPERROR_INTERNAL_ERROR = -21, // QP solver internal error
  QPERROR_INFEASIBLE =
      -22, // QP solver error: conclude QP formulation infeasible
  QPERROR_UNBOUNDED = -23,       // QP solver error: unbounded QP
  QPERROR_EXCEED_MAX_ITER = -24, // QP solver error: Exceed maximum iteration,
  QPERROR_NOTINITIALISED = -25,
  QPERROR_PREPARINGAUXILIARYQP = -26,
  QPERROR_AUXILIARYQPSOLVED = -27,
  QPERROR_PERFORMINGHOMOTOPY = -28,
  QPERROR_HOMOTOPYQPSOLVED = -29,
  QPERROR_UNKNOWN = -30,
  UNKNOWN = -99 // unknown error
};

/** Exceptions for different fatal errors. */
//@{
DECLARE_STD_EXCEPTION(SQP_EXCEPTION_INFEASIBLE);
DECLARE_STD_EXCEPTION(SQP_EXCEPTION_UNBOUNDED);
DECLARE_STD_EXCEPTION(SQP_EXCEPTION_MAXITER);
DECLARE_STD_EXCEPTION(SQP_EXCEPTION_INTERNAL_ERROR);
DECLARE_STD_EXCEPTION(SQP_EXCEPTION_UNKNOWN);
DECLARE_STD_EXCEPTION(SQP_EXCEPTION_PENALTY_TOO_LARGE);
DECLARE_STD_EXCEPTION(SQP_EXCEPTION_TRUST_REGION_TOO_SMALL);
//@}

enum ConstraintType
{
  BOUNDED_BELOW_AND_ABOVE = 5,
  IS_EQUALITY = -5,
  BOUNDED_ABOVE = 9,
  BOUNDED_BELOW = 1,
  UNBOUNDED = 0
};

enum ActivityStatus
{
  ACTIVE_ABOVE = 1,
  ACTIVE_BELOW = -1,
  ACTIVE_EQUALITY = -99,
  INACTIVE = 0
};

enum QpSolver
{
  QPOASES,
  QORE,
  GUROBI,
  CPLEX,
  SOLVER_UNDEFINED
};

typedef struct
{
  bool primal_feasibility = false;
  double primal_violation = 0.0;
  bool dual_feasibility = false;
  double dual_violation = 0.0;
  bool complementarity = false;
  double compl_violation = 0.0;
  bool stationarity = false;
  double stationarity_violation = 0.0;
  bool first_order_opt = false;
  double KKT_error = 0.0;
  bool Second_order_opt = false;
} OptimalityStatus;

typedef struct
{
  bool Update_A = false;
  bool Update_delta = false;
  bool Update_H = false;
  bool Update_penalty = false;
  bool Update_g = false;
  bool Update_bounds = false;
} UpdateFlags;
}
#endif /* SQPHOTSTART_TYPES_HPP*/
