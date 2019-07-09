#ifndef SQPHOTSTART_TYPES_HPP
#define SQPHOTSTART_TYPES_HPP

/* Declaration of all non-built in types (except for classes) */
namespace SQPhotstart {


    enum QPType {
        LP = 1,/** solving a linear program*/
        SOC = 2,/** solving second order correction**/
        QP = 3/**solving a regular qp subproblem **/
    };


    enum Exitflag {
        OPTIMAL = 0,
        EXCEED_MAX_ITER = 1,//exceeds the maximum number of iteration
        QPERROR_INTERNAL_ERROR = -1,//QP solver internal error
        QPERROR_INFEASIBLE = -2,//QP solver error: conclude QP formulation infeasible
        QPERROR_UNBOUNDED = -3, //QP solver error: unbounded QP
        QPERROR_EXCEED_MAX_ITER = -4,//QP solver error: Exceed maximum iteration,
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
    /** Type of default integer */
    typedef int Int;

    typedef struct {
        Index nCon;
        Index nVar;
        Index nnz_jac_g;
        Index nnz_h_lag;
    } Index_info;


    typedef struct {
        bool Update_A = true;
        bool Update_H = true;
        bool Update_penalty = true;
        bool Update_grad = true;
        bool Update_bounds = true;
    } UpdateFlags;

//    typedef struct{
//        Index* tmpRow = 0;
//        Index* tmpCol = 0;
//        Number* tmpEntry = 0;
//    }tmpMatrixData;
}
#endif /* SQPHOTSTART_TYPES_HPP*/