#include <algorithm>
#include <iterator>
#include <random>
#include "sqphot/SparseTripletMatrix.hpp"
#include "sqphot/Vector.hpp"
#include <stdlib.h> /* srand, rand */
#include <time.h>
#include "unit_test_utils.hpp"

using namespace RestartSqp;
using namespace std;

bool TEST_DENSE_SPARSE_MATRIX_CONVERSION(int rowNum, int colNum,
                                         const double* dense_matrix_in)
{

  shared_ptr<Vector> dense_matrix_out = make_shared<Vector>(rowNum * colNum);
  /**-------------------------------------------------------**/
  /**     Testing on Row Oriented Input Data                **/
  /**-------------------------------------------------------**/
  shared_ptr<SparseTripletMatrix> m_row_oriented =
      make_shared<SparseTripletMatrix>(dense_matrix_in, rowNum, colNum, true);

  m_row_oriented->get_dense_matrix(dense_matrix_out->get_non_const_values(), true);

  if (TEST_EQUAL_DOUBLE_ARRAY(dense_matrix_in, dense_matrix_out->get_values(),
                              rowNum * colNum, "row_oriented_matrix")) {
    printf("---------------------------------------------------------\n");
    printf("    Conversion between row-oriented dense matrix \n    "
           "and triplet sparse matrix test passed!   \n");
    printf("---------------------------------------------------------\n");
  } else {
    printf("---------------------------------------------------------\n");
    printf("    Conversion between row-oriented dense matrix \n    "
           "and triplet sparse matrix failed!   \n");
    printf("    Converting dense matrix to Sparse matrix FAILED!   \n");
    printf("\nMatrix in is \n");
    for (int i = 0; i < rowNum; i++) {
      for (int j = 0; j < colNum; j++)
        printf("%10e ", dense_matrix_in[i * colNum + j]);
      printf("\n");
    }

    printf("\nMatrix out is \n");
    for (int i = 0; i < rowNum; i++) {
      for (int j = 0; j < colNum; j++)
        printf("%10e ", dense_matrix_out->get_value(i * colNum + j));

      printf("\n");
    }

    printf("\n---------------------------------------------------------\n");
  }
  /**-------------------------------------------------------**/
  /**     Testing on Column Oriented Input Data             **/
  /**-------------------------------------------------------**/

  shared_ptr<SparseTripletMatrix> m_col_oriented =
      make_shared<SparseTripletMatrix>(dense_matrix_in, rowNum, colNum, false);

  dense_matrix_out->set_to_zero();
  m_col_oriented->get_dense_matrix(dense_matrix_out->get_non_const_values(), false);

  if (TEST_EQUAL_DOUBLE_ARRAY(dense_matrix_in, dense_matrix_out->get_values(),
                              rowNum * colNum, "col_oriented_matrix")) {
    printf("---------------------------------------------------------\n");
    printf("    Conversion between col-oriented dense matrix \n    "
           "and triplet sparse matrix test passed!   \n");
    printf("---------------------------------------------------------\n");
  } else {
    printf("---------------------------------------------------------\n");
    printf("    Conversion between col-oriented dense matrix \n    "
           "and triplet sparse matrix failed!   \n");
    printf("\nMatrix in is \n");
    for (int i = 0; i < rowNum; i++) {
      for (int j = 0; j < colNum; j++)
        printf("%10e ", dense_matrix_in[i * colNum + j]);
      printf("\n");
    }

    printf("\nMatrix out is \n");
    for (int i = 0; i < rowNum; i++) {
      for (int j = 0; j < colNum; j++)
        printf("%10e ", dense_matrix_out->get_value(i * colNum + j));

      printf("\n");
    }
    printf("\n---------------------------------------------------------\n");
  }
}

bool TEST_SPARSE_MATRIX_VECTOR_MULTIPLICATION(int rowNum, int colNum,
                                              const double* dense_matrix_in,
                                              shared_ptr<const Vector> vector)
{

  assert(colNum == vector->get_dim());

  shared_ptr<Vector> result_dense = make_shared<Vector>(rowNum);
  shared_ptr<Vector> result_sparse = make_shared<Vector>(rowNum);

  // compare it with matrix-vector multiplication in dense matrix
  result_dense->set_to_zero();
  for (int i = 0; i < rowNum; i++) {
    for (int j = 0; j < colNum; j++) {
      result_dense->add_number_to_element(i, dense_matrix_in[i * colNum + j] *
                                                 vector->get_value(j));
    }
  }

  shared_ptr<SparseTripletMatrix> sparse_matrix =
      make_shared<SparseTripletMatrix>(dense_matrix_in, rowNum, colNum, true);

  result_sparse->set_to_zero();
  sparse_matrix->multiply(vector, result_sparse);

  if (TEST_EQUAL_DOUBLE_ARRAY(result_sparse->get_values(),
                              result_dense->get_values(), rowNum)) {
    printf("---------------------------------------------------------\n");
    printf("   Matrix-vector multiplication for triplet \n"
           "   sparse matrix passed!   \n");
    printf("---------------------------------------------------------\n");
    return true;
  } else {
    printf("---------------------------------------------------------\n");
    printf("   Matrix-vector multiplication for triplet\n"
           "   sparse matrix failed!   \n");
    result_sparse->print("result_sparse");
    result_dense->print("result_dense");
    printf("\n---------------------------------------------------------\n");
    return false;
  }
}

bool TEST_TRANSPOSED_MATRIX_VECTOR_MULTIPLICATION(
    int rowNum, int colNum, const double* dense_matrix_in,
    shared_ptr<const Vector> vector)
{

  assert(rowNum == vector->get_dim());

  shared_ptr<Vector> result_dense = make_shared<Vector>(colNum);
  shared_ptr<Vector> result_sparse = make_shared<Vector>(colNum);

  // compare it with matrix-vector multiplication in dense matrix
  result_dense->set_to_zero();
  for (int i = 0; i < colNum; i++) {
    for (int j = 0; j < rowNum; j++) {
      result_dense->add_number_to_element(i, dense_matrix_in[j * colNum + i] *
                                                 vector->get_value(j));
    }
  }

  shared_ptr<SparseTripletMatrix> sparse_matrix =
      make_shared<SparseTripletMatrix>(dense_matrix_in, rowNum, colNum, true);
  sparse_matrix->multiply_transpose(vector, result_sparse);

  if (TEST_EQUAL_DOUBLE_ARRAY(result_sparse->get_values(),
                              result_dense->get_values(), colNum)) {
    printf("---------------------------------------------------------\n");
    printf("   Transposed matrix-vector multiplication for triplet\n"
           "   sparse matrix passed!   \n");

    printf("---------------------------------------------------------\n");
    return true;
  } else {
    printf("---------------------------------------------------------\n");
    printf("   Transposed matrix-vector multiplication for triplet\n"
           "   sparse matrix failed!   \n");
    result_sparse->print("result_sparse");
    result_dense->print("result_dense");
    printf("\n---------------------------------------------------------\n");
    return false;
  }
}

bool TEST_TRIPLET_HB_MATIRX_CONVERSION()
{
}

int main(int argc, char* argv[])
{

  /**-------------------------------------------------------**/
  /**     Generate a random matrix for testing              **/
  /**-------------------------------------------------------**/
  srand(time(NULL));
  std::random_device rd;
  //    std::mt19937 g(rd());

  int rowNum = rand() % 10 + 1;
  int colNum = rand() % 10 + 1;
  int EntryNum = rand() % (rowNum * colNum) + 1;

  std::vector<int> isNonzero;
  isNonzero.reserve(rowNum * colNum);
  for (int i = 0; i < rowNum * colNum; i++)
    isNonzero.push_back(i);

  //    std::shuffle(isNonzero.begin(), isNonzero.end(),g);

  auto dense_matrix_in = new double[rowNum * colNum]();
  for (int i = 0; i < rowNum * colNum; i++) {
    if (isNonzero.at(i) < EntryNum)
      // generate the entry between 1 to 10
      dense_matrix_in[i] = rand() % 10 + 1;
  }

  printf("\n=========================================================\n");
  printf("    Testing Methods for Triplet Sparse Matrix\n"
         "   on randomly generated data.");
  printf("\n=========================================================\n");

  /**-------------------------------------------------------**/
  /**            Dense-sparse Matrix Conversion             **/
  /**-------------------------------------------------------**/

  TEST_DENSE_SPARSE_MATRIX_CONVERSION(rowNum, colNum, dense_matrix_in);

  /**-------------------------------------------------------**/
  /**                 Matrix-vector Multiplication          **/
  /**-------------------------------------------------------**/

  shared_ptr<Vector> vector_mult = make_shared<Vector>(colNum);

  for (int i = 0; i < colNum; i++) {
    vector_mult->set_value(i, rand() % 10 + 1);
  }

  TEST_SPARSE_MATRIX_VECTOR_MULTIPLICATION(rowNum, colNum, dense_matrix_in,
                                           vector_mult);

  /**-------------------------------------------------------**/
  /**      Transposed Matrix-vector Multiplication          **/
  /**-------------------------------------------------------**/

  shared_ptr<Vector> vector_trans_mult = make_shared<Vector>(rowNum);

  for (int i = 0; i < rowNum; i++) {
    vector_trans_mult->set_value(i, rand() % 10 + 1);
  }

  TEST_TRANSPOSED_MATRIX_VECTOR_MULTIPLICATION(rowNum, colNum, dense_matrix_in,
                                               vector_trans_mult);

  isNonzero.clear();
  /**-------------------------------------------------------**/
  /**           Testing on Symmetric Matrix                 **/
  /**-------------------------------------------------------**/

  // generate a random symmetric matrix
  int dim = rand() % 10 + 2;
  isNonzero.reserve(dim * (dim - 1) / 2);
  for (int i = 0; i < dim * (dim - 1) / 2; i++)
    isNonzero.push_back(i);
  // std::shuffle(isNonzero.begin(), isNonzero.end(),g);

  int EntryNum_sym = rand() % (dim * (dim - 1) / 2) + 1; // with nonzero
                                                         // diagonal

  delete[] dense_matrix_in;
  dense_matrix_in = new double[dim * dim]();

  int iterator = 0;
  for (int i = 0; i < dim; i++) {
    for (int j = i; j < dim; j++) {
      if (i == j)
        dense_matrix_in[i * dim + j] = rand() % 10 + 1;
      else if (isNonzero.at(iterator) < EntryNum_sym) {
        dense_matrix_in[i * dim + j] = rand() % 10 + 1;
        dense_matrix_in[j * dim + i] = dense_matrix_in[i * dim + j];
        iterator++;
      } else
        iterator++;
    }
  }

  printf("\n=========================================================\n");
  printf("    Testing Methods for Symmetric Triplet Sparse Matrix\n"
         "   on randomly generated data.");
  printf("\n=========================================================\n");

  /**-------------------------------------------------------**/
  /**            Dense-sparse Matrix Conversion             **/
  /**-------------------------------------------------------**/

  TEST_DENSE_SPARSE_MATRIX_CONVERSION(dim, dim, dense_matrix_in);

  /**-------------------------------------------------------**/
  /**                 Matrix-vector Multiplication          **/
  /**-------------------------------------------------------**/

  shared_ptr<Vector> vector_sym_mult = make_shared<Vector>(dim);

  for (int i = 0; i < dim; i++) {
    vector_sym_mult->set_value(i, rand() % 10 + 1);
  }

  TEST_SPARSE_MATRIX_VECTOR_MULTIPLICATION(dim, dim, dense_matrix_in,
                                           vector_sym_mult);

  /**-------------------------------------------------------**/
  /**      Transposed Matrix-vector Multiplication          **/
  /**-------------------------------------------------------**/

  TEST_TRANSPOSED_MATRIX_VECTOR_MULTIPLICATION(dim, dim, dense_matrix_in,
                                               vector_sym_mult);

  delete[] dense_matrix_in;
  isNonzero.clear();

  return 0;
}
