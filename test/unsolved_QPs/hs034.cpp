
#include <qpOASES.hpp>

#include "hs034.hpp"



int main( )
{
    USING_NAMESPACE_QPOASES

    long i;
    int_t nWSR;
    real_t err, tic, toc;
    real_t *x = new real_t[13];
    real_t *y = new real_t[13+5];
    /* create sparse matrices */
    SymSparseMat *H = new SymSparseMat(13, 13, H_ir, H_jc, H_val);
    SparseMatrix *A = new SparseMatrix(5, 13, A_ir, A_jc, A_val);
    H->createDiagInfo();


    /* solve with sparse matrices */
    nWSR = 1000;
    QProblem QP(13, 5);
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
