/* Copyright (C) 2019
* All Rights Reserved.
*
* Authors: Xinyi Luo
* Date:2019-06
*/
#ifndef SQPHOTSTART_TYPES_HPP
#define SQPHOTSTART_TYPES_HPP

/* Declaration of all non-built in types (except for classes) */
namespace SQPhotstart {

//
//enum TestOption {
//    NO_TEST = -99,
//    TEST_1ST_ORDER = 1,
//    TEST_2ND_ORDER = 2,
//    TEST_STATIONARITY = -1,
//    TEST_COMPLEMENTARITY = -2,
//    TEST_FEASIBILITY = -3,
//    TEST_DUAL_FEASIBILITY = -4,
//    TEST_ALL = 0
//};

enum QPReturnType {
    UNSOLVED = -99,
    QP_OPTIMAL = 0,
    QP_NOTINITIALISED = -1,
    QP_EXCEED_MAX_ITER = 1,
    QP_INFEASIBLE = -2,
    QP_UNBOUNDED = -3,
    QP_PREPARINGAUXILIARYQP = -6,
    QP_AUXILIARYQPSOLVED = -7,
    QP_PERFORMINGHOMOTOPY = -8,
    QP_HOMOTOPYQPSOLVED = -9,
    QP_UNKNOWN_ERROR = 99
};


enum QPType {
    LP = 1,/** solving a linear program*/
    QP = 2/**solving a regular qp subproblem **/
};


enum Exitflag {
    OPTIMAL,
    INVALID_NLP,
    CONVERGE_TO_NONOPTIMAL,
    EXCEED_MAX_ITER, //exceeds the maximum number of iteration
    QPERROR_INTERNAL_ERROR, //QP solver internal error
    QPERROR_INFEASIBLE,//QP solver error: conclude QP formulation infeasible
    QPERROR_UNBOUNDED,  //QP solver error: unbounded QP
    QPERROR_EXCEED_MAX_ITER, //QP solver error: Exceed maximum iteration,
    QPERROR_NOTINITIALISED,
    QPERROR_PREPARINGAUXILIARYQP,
    QPERROR_AUXILIARYQPSOLVED,
    QPERROR_PERFORMINGHOMOTOPY,
    QPERROR_HOMOTOPYQPSOLVED,
    TRUST_REGION_TOO_SMALL,
    STEP_LARGER_THAN_TRUST_REGION,
    AUXINPUT_NOT_OPTIMAL = 99,//The input point in auxInput is not optimal when hotstart is enabled.
    UNKNOWN = -99//unknown error
};


enum ConstraintType {
    BOUNDED = 5,
    EQUAL = -5,
    BOUNDED_ABOVE = 9,
    BOUNDED_BELOW = 1,
    UNBOUNDED = 0
};

enum ActiveType {
    ACTIVE_ABOVE = 1,
    ACTIVE_BELOW = -1,
    ACTIVE_BOTH_SIDE = -99,
    INACTIVE = 0
};

enum Solver {
    QPOASES,
    QORE,
    GUROBI =2
};


typedef struct {
    int  nCon;
    int  nVar;
    int  nnz_jac_g;
    int  nnz_h_lag;
} Index_info;

typedef struct {
    bool primal_feasibility = false;
    double primal_violation = 0.0;
    bool dual_feasibility = false;
    double dual_violation = 0.0;
    bool complementarity = false;
    double compl_violation = 0.0;
    bool stationarity = false;
    double stationarity_violation = 0.0;
    bool first_order_opt = false;
    double KKT_error=0.0;
    bool Second_order_opt = false;
} OptimalityStatus;

typedef struct {
    bool Update_A = false;
    bool Update_delta = false;
    bool Update_H = false;
    bool Update_penalty = false;
    bool Update_g = false;
    bool Update_bounds= false;
} UpdateFlags;

}
#endif /* SQPHOTSTART_TYPES_HPP*/
