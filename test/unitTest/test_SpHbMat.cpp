#include <stdlib.h>     /* srand, rand */
#include <time.h>
#include <random>
#include <iterator>
#include <algorithm>

#include "unit_test_utils.hpp"
#include "sqphot/SpHbMat.hpp"
#include "sqphot/Vector.hpp"

using namespace SQPhotstart;
using namespace std;


bool TEST_DENSE_SPARSE_MATRIX_CONVERSION(int rowNum, int colNum, const double* dense_matrix_in) {
    shared_ptr<Vector> dense_matrix_out = make_shared<Vector>(rowNum * colNum);

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
        printf("\nMatrix in is \n");
        for(int i = 0; i < rowNum; i++) {
            for(int j = 0; j < colNum; j++)
                printf("%10e ",dense_matrix_in[i * colNum + j]);
            printf("\n");
        }

        printf("\nMatrix out is \n");
        for(int i = 0; i < rowNum; i++) {
            for(int j = 0; j < colNum; j++)
                printf("%10e ",dense_matrix_out->value(i * colNum + j));

            printf("\n");
        }

        printf("\n---------------------------------------------------------\n");

    }

    shared_ptr<SpHbMat> m_csr_row_oriented = make_shared<SpHbMat>(dense_matrix_in,
            rowNum, colNum, true, true);

    dense_matrix_out->set_to_zero();
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
                printf("%10e ",dense_matrix_out->value(i * colNum + j));

            printf("\n");
        }

        printf("\n---------------------------------------------------------\n");

    }


    shared_ptr<SpHbMat> m_csc_col_oriented = make_shared<SpHbMat>(dense_matrix_in,
            rowNum, colNum, false, true);

    dense_matrix_out->set_to_zero();
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
                printf("%10e ",dense_matrix_out->value(i * colNum + j));

            printf("\n");
        }

        printf("\n---------------------------------------------------------\n");

    }

    shared_ptr<SpHbMat> m_csr_col_oriented = make_shared<SpHbMat>(dense_matrix_in,
            rowNum, colNum, false, true);

    dense_matrix_out->set_to_zero();
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
                printf("%10e ",dense_matrix_out->value(i * colNum + j));

            printf("\n");
        }

        printf("\n---------------------------------------------------------\n");

    }

    return true;

}



bool
TEST_SPARSE_MATRIX_VECTOR_MULTIPLICATION(int rowNum, int colNum, const double* dense_matrix_in,
        shared_ptr<const Vector> vector) {

    assert(colNum == vector->dim());

    bool is_csc_passed, is_csr_passed;

    shared_ptr<Vector> result_dense = make_shared<Vector>(rowNum);
    shared_ptr<Vector> result_sparse = make_shared<Vector>(rowNum);

    //compare it with matrix-vector multiplication in dense matrix
    result_dense->set_to_zero();
    for(int i = 0; i<rowNum; i ++) {
        for(int j = 0; j<colNum; j++) {
            result_dense->add_number_to_element(i, dense_matrix_in[i*colNum+j]*
                                      vector->value(j));
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

    assert(rowNum == vector->dim());

    bool is_csc_passed, is_csr_passed;

    shared_ptr<Vector> result_dense = make_shared<Vector>(colNum);
    shared_ptr<Vector> result_sparse = make_shared<Vector>(colNum);

    //compare it with matrix-vector multiplication in dense matrix
    result_dense->set_to_zero();
    for(int i = 0; i<colNum; i ++) {
        for(int j = 0; j<rowNum; j++) {
            result_dense->add_number_to_element(i, dense_matrix_in[j*colNum+i]*
                                      vector->value(j));
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


bool TEST_TRIPLET_HB_MATIRX_CONVERSION(int rowNum, int colNum, const double* dense_matrix_in) {
    shared_ptr<Vector> dense_matrix_out = make_shared<Vector>(rowNum * colNum);

    /**-------------------------------------------------------**/
    /**             Without Identity Submatrix                **/
    /**-------------------------------------------------------**/
    shared_ptr<SpTripletMat> m_triplet_in = make_shared<SpTripletMat>(dense_matrix_in,
                                            rowNum, colNum, true);//condensed column sparse matrix

    auto m_csc = make_shared<SpHbMat>(rowNum, colNum,false) ;
    m_csc->setStructure(m_triplet_in);
    m_csc->setMatVal(m_triplet_in);
    auto m_triplet_csc_out = m_csc->convert_to_triplet();
    m_triplet_csc_out->get_dense_matrix(dense_matrix_out->values(),true);


    /**-------------------------------------------------------**/
    /**             Condensed Column matrix                   **/
    /**-------------------------------------------------------**/

    if(TEST_EQUAL_DOUBLE_ARRAY(dense_matrix_in, dense_matrix_out->values(),
                               m_triplet_in->EntryNum())) {
        printf("---------------------------------------------------------\n");
        printf("   Testing Triplet and condensed column sparse matrix\n"
               "    conversion without submatrix passed!       \n");
        printf("---------------------------------------------------------\n");
    } else {
        printf("---------------------------------------------------------\n");
        printf("---------------------------------------------------------\n");
        printf("   Testing Triplet and condensed column sparse matrix\n"
               "    conversion without submatrix FAILED!       \n");
        m_triplet_in->print("m_triplet_in");
        m_csc->print("m_triplet_out");
        m_triplet_csc_out->print("m_triplet_csc_out");
        printf("---------------------------------------------------------\n");
    }

    /**-------------------------------------------------------**/
    /**             Condensed Row matrix                      **/
    /**-------------------------------------------------------**/

    auto m_csr = make_shared<SpHbMat>(rowNum, colNum,true) ;
    m_csr->setStructure(m_triplet_in);
    m_csr->setMatVal(m_triplet_in);
    auto m_triplet_csr_out = m_csc->convert_to_triplet();
    m_triplet_csr_out->get_dense_matrix(dense_matrix_out->values(),true);

    if(TEST_EQUAL_DOUBLE_ARRAY(dense_matrix_in, dense_matrix_out->values(),
                               m_triplet_in->EntryNum())) {
        printf("---------------------------------------------------------\n");
        printf("   Testing Triplet and condensed row sparse matrix\n"
               "    conversion without submatrix passed!       \n");
        printf("---------------------------------------------------------\n");
    } else {
        printf("---------------------------------------------------------\n");
        printf("---------------------------------------------------------\n");
        printf("   Testing Triplet and condensed column sparse matrix\n"
               "    conversion without submatrix FAILED!       \n");
        m_triplet_in->print("m_triplet_in");
        m_csr->print("m_triplet_out");
        m_triplet_csr_out->print("m_triplet_csr_out");
        printf("---------------------------------------------------------\n");
    }
}


bool TEST_MATRIX_ALLOCATION_FROM_PERMUTATIONS() {


}



bool TEST_SET_STRUCTURE(int rowNum, int colNum, const double* dense_matrix_in) {





}

bool TEST_SET_MATRIX_VALUE() {
}




int main(int argc, char* argv[]) {


    /**-------------------------------------------------------**/
    /**     Generate a random matrix for testing              **/
    /**-------------------------------------------------------**/
    srand (time(NULL));
    std::random_device rd;
    //std::mt19937 g(rd());

    int rowNum = rand() %10+1;
    int colNum = rand() %10+1;
    int EntryNum = rand() %(rowNum*colNum)+1;

    std::vector<int> isNonzero;
    isNonzero.reserve(rowNum*colNum);
    for(int i = 0; i<rowNum*colNum; i++)
        isNonzero.push_back(i);

    //std::shuffle(isNonzero.begin(), isNonzero.end(),g);

    auto dense_matrix_in = new double[rowNum*colNum](); // this constructor initializes all values to 0.
    for(int i = 0; i < rowNum*colNum; i++) {
        if(isNonzero.at(i)<EntryNum)
            //generate the entry between 1 to 10
            dense_matrix_in[i] = rand() %10+1;
    }

    printf("\n=========================================================\n");
    printf("    Testing Methods for Harwell-Boeing Sparse Matrix\n"
           "   on randomly generated data.");
    printf("\n=========================================================\n");
    /**-------------------------------------------------------**/
    /**            Dense-sparse Matrix Conversion             **/
    /**-------------------------------------------------------**/


    TEST_DENSE_SPARSE_MATRIX_CONVERSION(rowNum,colNum,dense_matrix_in);

    /**-------------------------------------------------------**/
    /**                 Matrix-vector Multiplication          **/
    /**-------------------------------------------------------**/

    shared_ptr<Vector> vector_mult = make_shared<Vector>(colNum);

    for(int i = 0; i < colNum; i++) {
        vector_mult->set_value(i,rand() %10+1);
    }

    TEST_SPARSE_MATRIX_VECTOR_MULTIPLICATION(rowNum, colNum, dense_matrix_in, vector_mult);


    /**-------------------------------------------------------**/
    /**      Transposed Matrix-vector Multiplication          **/
    /**-------------------------------------------------------**/

    shared_ptr<Vector> vector_trans_mult = make_shared<Vector>(rowNum);

    for(int i = 0; i < rowNum; i++) {
        vector_trans_mult->set_value(i,rand() %10+1);
    }

    TEST_TRANSPOSED_MATRIX_VECTOR_MULTIPLICATION(rowNum, colNum, dense_matrix_in, vector_trans_mult);

    /**-------------------------------------------------------**/
    /**            Triplet Hb Matrix Conversion               **/
    /**-------------------------------------------------------**/

    TEST_TRIPLET_HB_MATIRX_CONVERSION(rowNum, colNum, dense_matrix_in);


    delete[] dense_matrix_in;
    isNonzero.clear();

    return 0;
}

