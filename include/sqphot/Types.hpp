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


enum TestOption {
    NO_TEST = -99,
    TEST_1ST_ORDER = 1,
    TEST_2ND_ORDER = 2,
    TEST_STATIONARITY = -1,
    TEST_COMPLEMENTARITY = -2,
    TEST_FEASIBILITY = -3,
    TEST_DUAL_FEASIBILITY = -4,
    TEST_ALL = 0
};

enum QPReturnType {
    UNSOLVED = -99,
    QP_OPTIMAL = 0,
    QP_NOTINITIALISED = -1,
    QP_EXCEED_MAX_ITER = 1,
    QP_INFEASIBLE = -2,
    QP_UNBOUNDED = -3,
    QP_UNKNOWN_ERROR = 99
};


enum QPType {
    LP = 1,/** solving a linear program*/
    SOC = 2,/** solving second order correction**/
    QP = 3/**solving a regular qp subproblem **/
};


enum Exitflag {
    OPTIMAL = 0,
    INVALID_NLP = 3,
    CONVERGE_TO_NONOPTIMAL = 2,
    EXCEED_MAX_ITER = 1,//exceeds the maximum number of iteration
    QPERROR_INTERNAL_ERROR = -1,//QP solver internal error
    QPERROR_INFEASIBLE = -2,//QP solver error: conclude QP formulation infeasible
    QPERROR_UNBOUNDED = -3, //QP solver error: unbounded QP
    QPERROR_EXCEED_MAX_ITER = -4,//QP solver error: Exceed maximum iteration,
    QPERROR_NOTINITIALISED = -5,
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


/** Type of all numbers */
typedef double Number;
/** Type of all indices of vectors, matrices etc */
typedef int Index;

typedef struct {
    Index nCon;
    Index nVar;
    Index nnz_jac_g;
    Index nnz_h_lag;
} Index_info;

typedef struct {
    bool primal_feasibility = false;
    bool dual_feasibility = false;
    bool complementarity = false;
    bool stationarity = false;
    bool first_order_opt = false;
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
