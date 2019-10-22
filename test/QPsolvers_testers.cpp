#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

extern "C" {
#include <qpsolver.h>
}
#include <iostream>
#include <memory>
#include <qpOASES.hpp>
#include <sqphot/SpHbMat.hpp>
#include <sqphot/Vector.hpp>
#include <sqphot/MessageHandling.hpp>

using namespace SQPhotstart;
shared_ptr<SpHbMat> convert_csr_to_csc(                           shared_ptr<const SpHbMat> csr_matrix) {

    double* dense_matrix;
    dense_matrix = new double[csr_matrix->RowNum()*csr_matrix->ColNum()]();//get the dense matrix out

    csr_matrix->get_dense_matrix(dense_matrix);

    shared_ptr<SpHbMat> csc_matrix = make_shared<SpHbMat>(dense_matrix, csr_matrix->RowNum(), csr_matrix->ColNum(), true,false);
    delete[] dense_matrix;

    return csc_matrix;
}

int main(int argc, char* argv[]) {

    FILE* file = fopen(argv[1], "r");
    char buffer[100];
    qp_int tmp_int;
    double tmp_double;
    int nVar, nCon, Hnnz, Annz;// number of variables, number of constraints, number of nonzeros
    string exitflag_qore, exitflag_qpOASES;
    int i;
    // for H and A;
    if(argc < 2) {
        printf("Missing Filename\n");
        return 1;
    }
    else if(file == NULL) perror("Error opening file");
    else {
        //begin Reading Data
        if(fgets(buffer,100,file)!=NULL) {
            sscanf(buffer,"%d",&nVar);
        }
        if(fgets(buffer,100,file)!=NULL) {
            sscanf(buffer,"%d",&nCon);
        }
        if(fgets(buffer,100,file)!=NULL) {
            sscanf(buffer,"%d",&Annz);
        }
        if(fgets(buffer,100,file)!=NULL) {
            sscanf(buffer,"%d",&Hnnz);
        }
    }

    shared_ptr<Vector> x_qore = make_shared<Vector>(nCon+nVar);
    shared_ptr<Vector> y_qore = make_shared<Vector>(nCon+nVar);

    shared_ptr<Vector> x_qpOASES = make_shared<Vector>(nVar);
    shared_ptr<Vector> y_qpOASES = make_shared<Vector>(nCon+nVar);

    shared_ptr<Vector> lb = make_shared<Vector>(nCon+nVar);
    shared_ptr<Vector> ub = make_shared<Vector>(nCon+nVar);
    shared_ptr<Vector> g = make_shared<Vector>(nVar);

    shared_ptr<SpHbMat> A_qore = make_shared<SpHbMat>(Annz, nCon, nVar, true);
    shared_ptr<SpHbMat> H_qore = make_shared<SpHbMat>(Hnnz, nVar, nVar, true);


    qp_int A_ir_qore[nCon+1];//row index
    qp_int A_jc_qore[Annz];//column index
    double A_val_qore[Annz];//matrix value

    qp_int H_ir_qore[nVar+1];//row index
    qp_int H_jc_qore[Hnnz];//column index
    double H_val_qore[Hnnz];//matrix value


    for(i = 0 ; i <nCon+nVar; i++) {
        if(fgets(buffer,100,file)!=NULL) {
            sscanf(buffer,"%lf",&tmp_double);
            lb->setValueAt(i,tmp_double);
        }
    }

    for(i = 0 ; i <nCon+nVar; i++) {
        if(fgets(buffer,100,file)!=NULL) {
            sscanf(buffer,"%lf",&tmp_double);
            ub->setValueAt(i,tmp_double);
        }
    }

    for(i = 0 ; i <nVar; i++) {
        if(fgets(buffer,100,file)!=NULL) {
            sscanf(buffer,"%lf",&tmp_double);
            g->setValueAt(i,tmp_double);
        }
    }

    for(i = 0 ; i <nCon+1; i++) {
        if(fgets(buffer,100,file)!=NULL) {
            sscanf(buffer,"%d",&tmp_int);
            A_qore->setRowIndexAt(i,tmp_int);
        }

    }


    for(i = 0 ; i <Annz; i++) {
        if(fgets(buffer,100,file)!=NULL) {
            sscanf(buffer,"%d",&tmp_int);
            A_qore->setColIndexAt(i,tmp_int);
        }
    }

    for(i = 0 ; i <Annz; i++) {
        if(fgets(buffer,100,file)!=NULL) {
            sscanf(buffer,"%lf",&tmp_double);
            A_qore->setMatValAt(i,tmp_double);
        }
    }

    for(i = 0 ; i <nVar+1; i++) {
        if(fgets(buffer,100,file)!=NULL) {
            sscanf(buffer,"%d",&tmp_int);
            H_qore->setRowIndexAt(i,tmp_int);
        }
    }

    for(i = 0 ; i <Hnnz; i++) {
        if(fgets(buffer,100,file)!=NULL) {
            sscanf(buffer,"%d",&tmp_int);
            H_qore->setColIndexAt(i,tmp_int);
        }
    }

    for(i = 0 ; i <Hnnz; i++) {
        if(fgets(buffer,100,file)!=NULL) {
            sscanf(buffer,"%lf",&tmp_double);
            H_qore->setMatValAt(i,tmp_double);
        }
    }



    //read matrix data
    ///////////////////////////////////////////////////////////
    //                     QORE                              //
    ///////////////////////////////////////////////////////////

    qp_int rv;
    QoreProblem * qp_qore = 0;
    rv = QPNew( &qp_qore, nVar, nCon, Annz, Hnnz);
    assert( rv == QPSOLVER_OK );
    assert( qp_qore!= 0 );
    QPSetInt(qp_qore, "prtfreq", -1);

    rv = QPSetData( qp_qore, nVar, nCon, A_qore->RowIndex(), A_qore->ColIndex(), A_qore->MatVal(), H_qore->RowIndex(), H_qore->ColIndex(), H_qore->MatVal());
    assert( rv == QPSOLVER_OK );
    rv = QPOptimize( qp_qore, lb->values(), ub->values(), g->values(), 0, 0 );
    assert( rv == QPSOLVER_OK );
    qp_int status;
    qp_int iter_count_qore;
    QPGetInt(qp_qore, "status", &status);
    QPGetInt(qp_qore, "itercount", &iter_count_qore);

    switch (status) {
    case QPSOLVER_ITER_LIMIT:
        exitflag_qore = "EXCEED_ITER_LIMIT";
        break;
    case QPSOLVER_OPTIMAL:
        exitflag_qore = "OPTIMAL";
        break;
    case QPSOLVER_INFEASIBLE:
        exitflag_qore = "INFEASIBLE";
        break;
    case QPSOLVER_UNBOUNDED:
        exitflag_qore = "UNBOUNDED";
        break;
    default:
        exitflag_qore = "UNKNOWN";
    }


    rv = QPGetDblVector( qp_qore, "primalsol", x_qore->values());
    rv = QPGetDblVector( qp_qore, "dualsol", y_qore->values());

    //x_qore->print("x_qore");
    //y_qore->print("y_qore");

    shared_ptr<Vector> Hx = make_shared<Vector>(nVar);
    H_qore->times(x_qore,Hx);
    double obj_qore = Hx->times(x_qore)*0.5+g->times(x_qore);

    ///////////////////////////////////////////////////////////
    //                     QPOASES                           //
    ///////////////////////////////////////////////////////////
    auto A_qpOASES = convert_csr_to_csc(A_qore);
    auto H_qpOASES = convert_csr_to_csc(H_qore);

    /* create sparse matrices */
    qpOASES::SymSparseMat *H  = new qpOASES::SymSparseMat(nVar,nVar, H_qpOASES->RowIndex(),
            H_qpOASES->ColIndex(),H_qpOASES->MatVal());
    qpOASES::SparseMatrix *A = new qpOASES::SparseMatrix(nCon,nVar, A_qpOASES->RowIndex(),A_qpOASES->ColIndex(),
            A_qpOASES->MatVal());

    H->createDiagInfo();
//    A->print("A");
//    H->print("H");
    qpOASES::Options qp_options;
    qp_options.setToReliable();
    qp_options.printLevel = qpOASES::PL_NONE;



    int iter_qpOASES = 1000; //maximum number of iterations

    qpOASES::QProblem qp_qpOASES(nVar,nCon);
    qp_qpOASES.setOptions(qp_options);
    qp_qpOASES.init(H,g->values(),A,lb->values(),ub->values(),lb->values()+nVar,ub->values()+nVar,iter_qpOASES,0);


    qp_qpOASES.getPrimalSolution(x_qpOASES->values());
    qp_qpOASES.getDualSolution(y_qpOASES->values());


    if (qp_qpOASES.isInfeasible()) {
        exitflag_qpOASES = "INFEASIBLE";
    }
    else if (qp_qpOASES.isUnbounded()) {
        exitflag_qpOASES = "UNBOUNDED";
    }
    else if(qp_qpOASES.isSolved()) {
        exitflag_qpOASES = "OPTIMAL";
    }
    else
        switch (qp_qpOASES.getStatus()) {
        case qpOASES::QPS_NOTINITIALISED:
            exitflag_qpOASES = "QPERROR_NOTINITIALISED";
        case qpOASES::QPS_PREPARINGAUXILIARYQP:
            exitflag_qpOASES = "QPERROR_PREPARINGAUXILIARYQP";
        case qpOASES::QPS_AUXILIARYQPSOLVED:
            exitflag_qpOASES = "QPERROR_AUXILIARYQPSOLVED";
        case qpOASES::QPS_PERFORMINGHOMOTOPY:
            exitflag_qpOASES = "QPERROR_PERFORMINGHOMOTOPY";
        case qpOASES::QPS_HOMOTOPYQPSOLVED:
            exitflag_qpOASES = "QPERROR_HOMOTOPYQPSOLVED";
        }

    //x_qpOASES->print("x_qpOASES");
    //y_qpOASES->print("y_qpOASES");


    ///////////////////////////////////////////////////////////
    //                     PRINT OUT RESULTS                 //
    ///////////////////////////////////////////////////////////

    string* pname = new string(argv[1]);
    std::size_t found = pname->find_last_of("/\\");
    printf(DOUBLE_LONG_DIVIDER);
    printf("                         Solving QP Problem stored in %10s\n",
           pname->substr(found + 1).c_str());
    printf(DOUBLE_LONG_DIVIDER);
    printf("%20s    %23s    %23s\n","","QORE","QPOASES");
    printf("%20s    %23s    %23s\n","Exitflag",exitflag_qore.c_str(),exitflag_qpOASES.c_str());
    printf("%20s    %23d    %23d\n","Iteration",iter_count_qore,(int)iter_qpOASES);
    printf("%20s    %23.16e    %23.16e\n","Objective",obj_qore,qp_qpOASES.getObjVal());
//    printf("%20s    %23.16e    %23.16e\n","KKT Error","QORE","QPOASES");
    printf(DOUBLE_LONG_DIVIDER);

    delete pname;
    delete H;
    delete A;
    delete qp_qore;
    fclose(file);
    return 0;
}


