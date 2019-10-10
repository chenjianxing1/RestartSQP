#include <sqphot/Matrix.hpp>
#include <sqphot/Vector.hpp>
#include <stdlib.h>     /* srand, rand */
#include <time.h>
#include <random>
#include <iterator>

using namespace SQPhotstart;
template <typename T>
bool TEST_EQUAL(T a, T b) {
    return a == b;
}

bool TEST_EQUAL_INT_ARRAY(int* a, int* b, int length, const char* name = NULL) {
    for(int i = 0; i < length; i++) {
        if(a[i] != b[i]) {
            printf("TEST_EQUAL_INT_ARRAY for %s failed at index %d\n", name, i);
            return false;
        }
    }
    return true;

}

bool TEST_EQUAL_DOUBLE_ARRAY(double* a, double* b, int length, const char* name = NULL) {
    for(int i = 0; i < length; i++) {
        if(a[i] != b[i]) {
            printf("TEST_EQUAL_DOUBLE_ARRAY for %s failed at index %d\n", name, i);
            return false;
        }
    }
    return true;
}


bool TEST_DENSE_SPARSE_MATRIX_CONVERSION(int rowNumber, int colNumber,double* dense_matrix_in) {

    auto dense_matrix_out = new double[rowNumber*colNumber]();



    //    shared_ptr<SpHbMat> m_csr_row_oriented = make_shared<SpHbMat>(dense_matrix_in,rowNumber,colNumber,true,true);
    shared_ptr<SpHbMat> m_csc_row_oriented = make_shared<SpHbMat>(dense_matrix_in,
            rowNumber,colNumber,true,false);

    m_csc_row_oriented->get_dense_matrix(dense_matrix_out);

    if(TEST_EQUAL_DOUBLE_ARRAY(dense_matrix_in,dense_matrix_out, rowNumber*colNumber,
                               "csc_row_oriented_matrix")) {
        printf("---------------------------------------------------------\n");
        printf("    Conversion between row-oriented dense matrix \n    "
               "    and condensed column sparse matrix test passed!   \n");
        printf("---------------------------------------------------------\n");
    }
    else {

        printf("---------------------------------------------------------\n");
        printf("    Conversion between row-oriented dense matrix \n"
               "    and condensed column sparse matrix failed!   \n");
        printf("    Converting dense matrix to Sparse matrix FAILED!   \n");
        printf("\nMatrix in is \n");
        for(int i = 0; i <rowNumber; i++) {
            for(int j = 0; j <colNumber; j++)
                printf("%10e ",dense_matrix_in[i* colNumber+j]);
            printf("\n");
        }

        printf("\nMatrix out is \n");
        for(int i = 0; i <rowNumber; i++) {
            for(int j = 0; j <colNumber; j++)
                printf("%10e ",dense_matrix_out[i* colNumber+j]);

            printf("\n");
        }

        printf("\n---------------------------------------------------------\n");

    }

    shared_ptr<SpHbMat> m_csr_col_oriented = make_shared<SpHbMat>(dense_matrix_in,
            rowNumber,colNumber,false,false);
    shared_ptr<SpHbMat> m_csc_col_oriented = make_shared<SpHbMat>(dense_matrix_in,
            rowNumber,colNumber,false,true);

    delete [] dense_matrix_out;
    dense_matrix_out = NULL;
    return true;

}


bool TEST_TRIPLET_HB_MATIRX_CONVERSION() {
}


bool TEST_SPARSE_MATRIX_VECTOR_MULTIPLICATION(shared_ptr<SpHbMat> matrix,
        shared_ptr<Vector> vector) {

    assert(matrix->ColNum() == vector->Dim());
    shared_ptr<Vector> result_sparse = make_shared<Vector>(matrix->RowNum());
    matrix->times(vector,result_sparse);

    //compare it with matrix-vector multiplication in dense matrix

    double* dense_matrix;
    dense_matrix = new double[matrix->RowNum()*matrix->ColNum()]();

    matrix->get_dense_matrix(dense_matrix);

    shared_ptr<Vector> result_dense = make_shared<Vector>(matrix->RowNum());

    for(int i = 0; i<matrix->RowNum(); i ++) {
        for(int j = 0; j<matrix->ColNum(); j++) {
            result_dense->addNumberAt(i, dense_matrix[i*matrix->ColNum()+j]*
                                      vector->values(j));
        }
    }



    if(TEST_EQUAL_DOUBLE_ARRAY(result_dense->values(), result_sparse->values(),
                               result_sparse->Dim())) {

        printf("---------------------------------------------------------\n");
        printf("   Matrix-vector multiplication test passed!   \n");
        printf("---------------------------------------------------------\n");
    }
    else {
        result_sparse->print("result_sparse");
        result_dense->print("result_dense");

    }
    delete [] dense_matrix;


}

bool TEST_MATRIX_ALLOCATION_FROM_PERMUTATIONS() {

}




int main(int argc, char* argv[]) {

    /**
     * Generate a random matrix for testing
     */
    srand (time(NULL));
    std::random_device rd;
    std::mt19937 g(rd());

    int rowNumber = rand() %10+1;
    int colNumber = rand() %10+1;
    int EntryNum = rand() %(rowNumber*colNumber)+1;

    std::vector<int> isNonzero;
    isNonzero.reserve(rowNumber*colNumber);
for(int i = 0; i<rowNumber*colNumber; i++)
        isNonzero.push_back(i);

    std::shuffle(isNonzero.begin(), isNonzero.end(),g);

    double* dense_matrix_in;
    shared_ptr<Vector> vector_in = make_shared<Vector>(colNumber);



    dense_matrix_in = new double[rowNumber*colNumber]();
    for(int i = 0; i < rowNumber*colNumber; i++) {
        if(isNonzero.at(i)<EntryNum)
            //generate the entry between 1 to 10
            dense_matrix_in[i] = rand() %10+1;
    }


    for(int i = 0; i < colNumber; i++) {
        vector_in->setValueAt(i,rand() %10+1);
    }


    TEST_DENSE_SPARSE_MATRIX_CONVERSION(rowNumber,colNumber,dense_matrix_in);


    shared_ptr<SpHbMat> matrix_csc = make_shared<SpHbMat>(dense_matrix_in,rowNumber,
      colNumber,true,false);

    TEST_SPARSE_MATRIX_VECTOR_MULTIPLICATION(matrix_csc,vector_in);

    isNonzero.clear();

}

