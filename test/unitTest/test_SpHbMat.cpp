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


bool TEST_DENSE_SPARSE_MATRIX_CONVERSION(int rowNum, int colNum, double* dense_matrix_in) {

    shared_ptr<Vector> dense_matrix_out = make_shared<Vector>(rowNum * colNum);



    //    shared_ptr<SpHbMat> m_csr_row_oriented = make_shared<SpHbMat>(dense_matrix_in,rowNumber,colNumber,true,true);
    shared_ptr<SpHbMat> m_csc_row_oriented = make_shared<SpHbMat>(dense_matrix_in,
            rowNum, colNum, true, false);

    m_csc_row_oriented->get_dense_matrix(dense_matrix_out->values());

    if(TEST_EQUAL_DOUBLE_ARRAY(dense_matrix_in, dense_matrix_out->values(),
                               rowNum * colNum, "csc_row_oriented_matrix")) {
        printf("---------------------------------------------------------\n");
        printf("    Conversion between row-oriented dense matrix \n    "
               "and condensed column sparse matrix test passed!   \n");
        printf("---------------------------------------------------------\n");
    }
    else {

        printf("---------------------------------------------------------\n");
        printf("    Conversion between row-oriented dense matrix \n    "
               "and condensed column sparse matrix failed!   \n");
        printf("    Converting dense matrix to Sparse matrix FAILED!   \n");
        printf("\nMatrix in is \n");
        for(int i = 0; i < rowNum; i++) {
            for(int j = 0; j < colNum; j++)
                printf("%10e ",dense_matrix_in[i * colNum + j]);
            printf("\n");
        }

        printf("\nMatrix out is \n");
        for(int i = 0; i < rowNum; i++) {
            for(int j = 0; j < colNum; j++)
                printf("%10e ",dense_matrix_out->values(i * colNum + j));

            printf("\n");
        }

        printf("\n---------------------------------------------------------\n");

    }

    shared_ptr<SpHbMat> m_csr_row_oriented = make_shared<SpHbMat>(dense_matrix_in,
            rowNum, colNum, true, true);

    dense_matrix_out->set_zeros();
    m_csr_row_oriented->get_dense_matrix(dense_matrix_out->values());

    if(TEST_EQUAL_DOUBLE_ARRAY(dense_matrix_in, dense_matrix_out->values(),
                               rowNum * colNum, "csr_row_oriented_matrix")) {
        printf("---------------------------------------------------------\n");
        printf("    Conversion between row-oriented dense matrix \n    "
               "and condensed row sparse matrix test passed!   \n");
        printf("---------------------------------------------------------\n");
    }
    else {

        printf("---------------------------------------------------------\n");
        printf("    Conversion between row-oriented dense matrix \n    "
               "and condensed row sparse matrix failed!   \n");
        printf("    Converting dense matrix to Sparse matrix FAILED!   \n");
        printf("\nMatrix in is \n");
        for(int i = 0; i < rowNum; i++) {
            for(int j = 0; j < colNum; j++)
                printf("%10e ",dense_matrix_in[i * colNum + j]);
            printf("\n");
        }

        printf("\nMatrix out is \n");
        for(int i = 0; i < rowNum; i++) {
            for(int j = 0; j < colNum; j++)
                printf("%10e ",dense_matrix_out->values(i * colNum + j));

            printf("\n");
        }

        printf("\n---------------------------------------------------------\n");

    }


    shared_ptr<SpHbMat> m_csc_col_oriented = make_shared<SpHbMat>(dense_matrix_in,
            rowNum, colNum, false, true);

    dense_matrix_out->set_zeros();
    m_csc_col_oriented->get_dense_matrix(dense_matrix_out->values(),false);


    if(TEST_EQUAL_DOUBLE_ARRAY(dense_matrix_in, dense_matrix_out->values(),
                               rowNum * colNum, "csc_col_oriented_matrix")) {
        printf("---------------------------------------------------------\n");
        printf("    Conversion between col-oriented dense matrix \n    "
               "and condensed column sparse matrix test passed!   \n");
        printf("---------------------------------------------------------\n");
    }
    else {

        printf("---------------------------------------------------------\n");
        printf("    Conversion between col-oriented dense matrix \n    "
               "and condensed column sparse matrix failed!   \n");
        printf("    Converting dense matrix to Sparse matrix FAILED!   \n");
        printf("\nMatrix in is \n");
        for(int i = 0; i < rowNum; i++) {
            for(int j = 0; j < colNum; j++)
                printf("%10e ",dense_matrix_in[i * colNum + j]);
            printf("\n");
        }

        printf("\nMatrix out is \n");
        for(int i = 0; i < rowNum; i++) {
            for(int j = 0; j < colNum; j++)
                printf("%10e ",dense_matrix_out->values(i * colNum + j));

            printf("\n");
        }

        printf("\n---------------------------------------------------------\n");

    }

    shared_ptr<SpHbMat> m_csr_col_oriented = make_shared<SpHbMat>(dense_matrix_in,
            rowNum, colNum, false, true);

    dense_matrix_out->set_zeros();
    m_csr_col_oriented->get_dense_matrix(dense_matrix_out->values(),false);


    if(TEST_EQUAL_DOUBLE_ARRAY(dense_matrix_in, dense_matrix_out->values(),
                               rowNum * colNum, "csr_col_oriented_matrix")) {
        printf("---------------------------------------------------------\n");
        printf("    Conversion between col-oriented dense matrix \n    "
               "and condensed row sparse matrix test passed!   \n");
        printf("---------------------------------------------------------\n");
    }
    else {

        printf("---------------------------------------------------------\n");
        printf("    Conversion between col-oriented dense matrix \n    "
               "and condensed row sparse matrix failed!   \n");
        printf("    Converting dense matrix to Sparse matrix FAILED!   \n");
        printf("\nMatrix in is \n");
        for(int i = 0; i < rowNum; i++) {
            for(int j = 0; j < colNum; j++)
                printf("%10e ",dense_matrix_in[i * colNum + j]);
            printf("\n");
        }

        printf("\nMatrix out is \n");
        for(int i = 0; i < rowNum; i++) {
            for(int j = 0; j < colNum; j++)
                printf("%10e ",dense_matrix_out->values(i * colNum + j));

            printf("\n");
        }

        printf("\n---------------------------------------------------------\n");

    }

    return true;

}


bool TEST_TRIPLET_HB_MATIRX_CONVERSION() {
}


bool
TEST_SPARSE_MATRIX_VECTOR_MULTIPLICATION(int rowNum, int colNum, const double* dense_matrix_in,
        shared_ptr<const Vector> vector) {

    assert(colNum == vector->Dim());

    bool is_csc_passed, is_csr_passed;

    shared_ptr<Vector> result_dense = make_shared<Vector>(rowNum);
    shared_ptr<Vector> result_sparse = make_shared<Vector>(rowNum);

    //compare it with matrix-vector multiplication in dense matrix
    for(int i = 0; i<rowNum; i ++) {
        for(int j = 0; j<colNum; j++) {
            result_dense->addNumberAt(i, dense_matrix_in[i*colNum+j]*
                                      vector->values(j));
        }
    }


    shared_ptr<SpHbMat> m_csc = make_shared<SpHbMat>(dense_matrix_in,
                                rowNum, colNum, true, false);
    m_csc->times(vector,result_sparse);


    if(TEST_EQUAL_DOUBLE_ARRAY(result_sparse->values(),result_dense->values(),rowNum)) {
        printf("---------------------------------------------------------\n");
        printf("   Matrix-vector multiplication for condensed column \n"
               "   sparse matrix passed!   \n");
        printf("---------------------------------------------------------\n");
        is_csc_passed = true;
    }
    else {
        printf("---------------------------------------------------------\n");
        printf("   Matrix-vector multiplication for condensed column\n"
               "   sparse matrix failed!   \n");
        result_sparse->print("result_sparse");
        result_dense->print("result_dense");
        printf("\n---------------------------------------------------------\n");
        is_csc_passed = false;

    }

    shared_ptr<SpHbMat> m_csr = make_shared<SpHbMat>(dense_matrix_in,
                                rowNum, colNum, true, true);
    m_csr->times(vector,result_sparse);

    if(TEST_EQUAL_DOUBLE_ARRAY(result_sparse->values(),result_dense->values(),rowNum)) {
        printf("---------------------------------------------------------\n");
        printf("   Matrix-vector multiplication for condensed row \n"
               "   sparse matrix passed!   \n");
        printf("---------------------------------------------------------\n");
        is_csr_passed = true;
    }
    else {
        printf("---------------------------------------------------------\n");
        printf("   Matrix-vector multiplication for condensed row\n"
               "   sparse matrix failed!   \n");
        result_sparse->print("result_sparse");
        result_dense->print("result_dense");
        printf("\n---------------------------------------------------------\n");
        is_csr_passed = false;
    }

    if(is_csc_passed&&is_csr_passed)
        return true;
    else
        return false;

}


bool TEST_TRANSPOSED_MATRIX_VECTOR_MULTIPLICATION(int rowNum, int colNum, const double* dense_matrix_in,
        shared_ptr<const Vector> vector) {

    assert(rowNum == vector->Dim());

    bool is_csc_passed, is_csr_passed;

    shared_ptr<Vector> result_dense = make_shared<Vector>(colNum);
    shared_ptr<Vector> result_sparse = make_shared<Vector>(colNum);

    //compare it with matrix-vector multiplication in dense matrix
    for(int i = 0; i<colNum; i ++) {
        for(int j = 0; j<rowNum; j++) {
            result_dense->addNumberAt(i, dense_matrix_in[j*colNum+i]*
                                      vector->values(j));
        }
    }

    shared_ptr<SpHbMat> m_csc = make_shared<SpHbMat>(dense_matrix_in,
                                rowNum, colNum, true, false);
    m_csc->transposed_times(vector,result_sparse);


    if(TEST_EQUAL_DOUBLE_ARRAY(result_sparse->values(),result_dense->values(),colNum)) {
        printf("---------------------------------------------------------\n");
        printf("   Transposed matrix-vector multiplication for condensed \n"
               "   column sparse matrix passed!   \n");

        printf("---------------------------------------------------------\n");
        is_csc_passed = true;
    }
    else {
        printf("---------------------------------------------------------\n");
        printf("   Transposed matrix-vector multiplication for condensed \n"
               "   column sparse matrix failed!   \n");
        result_sparse->print("result_sparse");
        result_dense->print("result_dense");
        printf("\n---------------------------------------------------------\n");
        is_csc_passed = false;

    }

    shared_ptr<SpHbMat> m_csr = make_shared<SpHbMat>(dense_matrix_in,
                                rowNum, colNum, true, true);
    m_csr->transposed_times(vector,result_sparse);

    if(TEST_EQUAL_DOUBLE_ARRAY(result_sparse->values(),result_dense->values(),colNum)) {
        printf("---------------------------------------------------------\n");
        printf("   Transposed matrix-vector multiplication for condensed \n"
               "   row sparse matrix passed!   \n");

        printf("---------------------------------------------------------\n");
        is_csr_passed = true;
    }
    else {
        printf("---------------------------------------------------------\n");
        printf("   Transposed matrix-vector multiplication for condensed \n"
               "   row sparse matrix failed!   \n");
        result_sparse->print("result_sparse");
        result_dense->print("result_dense");
        printf("\n---------------------------------------------------------\n");
        is_csr_passed = false;
    }

    if(is_csc_passed&&is_csr_passed)
        return true;
    else
        return false;
}

bool TEST_MATRIX_ALLOCATION_FROM_PERMUTATIONS() {

}




int main(int argc, char* argv[]) {

    /**-------------------------------------------------------**/
    /**     Generate a random matrix for testing              **/
    /**-------------------------------------------------------**/
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

    dense_matrix_in = new double[rowNumber*colNumber]();
    for(int i = 0; i < rowNumber*colNumber; i++) {
        if(isNonzero.at(i)<EntryNum)
            //generate the entry between 1 to 10
            dense_matrix_in[i] = rand() %10+1;
    }


    /**-------------------------------------------------------**/
    /**            Dense-sparse Matrix Conversion             **/
    /**-------------------------------------------------------**/


    TEST_DENSE_SPARSE_MATRIX_CONVERSION(rowNumber,colNumber,dense_matrix_in);

    /**-------------------------------------------------------**/
    /**                 Matrix-vector Multiplication          **/
    /**-------------------------------------------------------**/

    shared_ptr<Vector> vector_mult = make_shared<Vector>(colNumber);

    for(int i = 0; i < colNumber; i++) {
        vector_mult->setValueAt(i,rand() %10+1);
    }

    TEST_SPARSE_MATRIX_VECTOR_MULTIPLICATION(rowNumber, colNumber, dense_matrix_in, vector_mult);


    /**-------------------------------------------------------**/
    /**      Transposed Matrix-vector Multiplication          **/
    /**-------------------------------------------------------**/

    shared_ptr<Vector> vector_trans_mult = make_shared<Vector>(rowNumber);

    for(int i = 0; i < rowNumber; i++) {
        vector_trans_mult->setValueAt(i,rand() %10+1);
    }

    TEST_TRANSPOSED_MATRIX_VECTOR_MULTIPLICATION(rowNumber, colNumber, dense_matrix_in, vector_trans_mult);


    delete[] dense_matrix_in;
    isNonzero.clear();

}

