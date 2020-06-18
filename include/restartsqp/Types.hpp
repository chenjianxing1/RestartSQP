/* Copyright (C) 2019
 * All Rights Reserved.
 *
 * Authors: Xinyi Luo
 * Date:2019-06
 */
#ifndef SQPHOTSTART_TYPES_HPP
#define SQPHOTSTART_TYPES_HPP

// This is necessary for AW's current version of QtCreator (otherwise the error
// parser does not understand ostream correctly
#include "IpException.hpp"
#include <ostream>

/** Declaration of all non-built in types (except for classes) */

namespace RestartSqp {

#include "restartsqp/ReturnCodes.h"

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

/** Exceptions for different fatal errors. */
//@{
DECLARE_STD_EXCEPTION(SQP_EXCEPTION_INFEASIBLE);
DECLARE_STD_EXCEPTION(SQP_EXCEPTION_UNBOUNDED);
DECLARE_STD_EXCEPTION(SQP_EXCEPTION_MAXITER);
DECLARE_STD_EXCEPTION(SQP_EXCEPTION_INTERNAL_ERROR);
DECLARE_STD_EXCEPTION(SQP_EXCEPTION_UNKNOWN);
DECLARE_STD_EXCEPTION(SQP_EXCEPTION_PENALTY_TOO_LARGE);
DECLARE_STD_EXCEPTION(SQP_EXCEPTION_TRUST_REGION_TOO_SMALL);
DECLARE_STD_EXCEPTION(SQP_EXCEPTION_QP_SOLVER_FAILS);
DECLARE_STD_EXCEPTION(SQP_INVALID_WORKING_SET);
DECLARE_STD_EXCEPTION(SQP_INVALID_STARING_POINT);
DECLARE_STD_EXCEPTION(SQP_NLP);
//@}

enum ConstraintType
{
  BOUNDED_BELOW_AND_ABOVE = 5,
  IS_EQUALITY = -5,
  BOUNDED_ABOVE = 9,
  BOUNDED_BELOW = 1,
  UNBOUNDED = 0
};

enum QpSolver
{
  QPOASES,
  QORE,
  GUROBI,
  CPLEX,
  SOLVER_UNDEFINED
};

enum ActivityStatus
{
  ACTIVE_ABOVE = 1,
  ACTIVE_BELOW = -1,
  ACTIVE_EQUALITY = -2,
  INACTIVE = 0,
  ACTIVE_INVALID = -3
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
} // namespace RestartSqp
#endif /* SQPHOTSTART_TYPES_HPP*/
