/* Copyright (C) 2019
 * All Rights Reserved.
 *
 * Authors: Xinyi Luo
 * Date:2019-06
 */

#ifndef SQP_RETURNCODES_H
#define SQP_RETURNCODES_H
// Include file for the return codes (for C and C++)
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
  EXCEED_MAX_LAZY_NLP_SOLVES = -9,
  ERROR_IN_LAZY_NLP_UPDATE = -10,
  INVALID_INITIAL_WORKING_SET = -11,
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
  UNKNOWN_EXIT_STATUS = -99 // unknown error
};
#endif
