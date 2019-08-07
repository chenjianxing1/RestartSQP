
#include <qpOASES.hpp>

#include "hs067_unbounded.hpp"



int main( )
{
    USING_NAMESPACE_QPOASES

    long i;
    int_t nWSR;
    real_t err, tic, toc;
    real_t *x = new real_t[52];
    real_t *y = new real_t[52+21];
    /* create sparse matrices */
    SymSparseMat *H = new SymSparseMat(52, 52, H_ir, H_jc, H_val);
    SparseMatrix *A = new SparseMatrix(21, 52, A_ir, A_jc, A_val);
    H->createDiagInfo();


    /* solve with sparse matrices */
    nWSR = 1000;
    QProblem QP(52, 21);
    Options options;
    options.printLevel = PL_DEBUG_ITER;
    QP.setOptions(options);
    tic = getCPUtime();
    QP.init(H, g, A, lb, ub, lbA, ubA, nWSR, 0);
    toc = getCPUtime();
    QP.getPrimalSolution(x);
    QP.getDualSolution(y);

    fprintf(stdFile, "Solved sparse problem in %d iterations, %.3f seconds.\n", (int)nWSR, toc-tic);


    delete H;
    delete A;

    delete[] y;
    delete[] x;

    return 0;
}


/*
 *	end of file
 */
