#include <algorithm>
#include <iterator>
#include <random>
#include <stdlib.h> /* srand, rand */
#include <time.h>

#include "restartsqp/SparseHbMatrix.hpp"
#include "restartsqp/Vector.hpp"
#include "unit_test_utils.hpp"

using namespace RestartSqp;
using namespace std;

bool TEST_DENSE_SPARSE_MATRIX_CONVERSION(int num_rows, int num_columns,
                                         const double* dense_matrix_in)
{
  shared_ptr<Vector> dense_matrix_out =
      make_shared<Vector>(num_rows * num_columns);

  bool is_compressed_row = false;
  bool row_oriented = true;
  shared_ptr<SparseHbMatrix> m_csc_row_oriented =
      make_shared<SparseHbMatrix>(num_rows, num_columns, is_compressed_row);
  m_csc_row_oriented->copy_from_dense_matrix(
      dense_matrix_in, num_rows, num_columns, row_oriented, is_compressed_row);

  m_csc_row_oriented->get_dense_matrix(
      dense_matrix_out->get_non_const_values());

  if (TEST_EQUAL_DOUBLE_ARRAY(dense_matrix_in, dense_matrix_out->get_values(),
                              num_rows * num_columns,
                              "csc_row_oriented_matrix")) {
    printf("---------------------------------------------------------\n");
    printf("    Conversion between row-oriented dense matrix \n    "
           "and condensed column sparse matrix test passed!   \n");
    printf("---------------------------------------------------------\n");
  } else {
    printf("---------------------------------------------------------\n");
    printf("    Conversion between row-oriented dense matrix \n    "
           "and condensed column sparse matrix failed!   \n");
    printf("\nMatrix in is \n");
    for (int i = 0; i < num_rows; i++) {
      for (int j = 0; j < num_columns; j++)
        printf("%10e ", dense_matrix_in[i * num_columns + j]);
      printf("\n");
    }

    printf("\nMatrix out is \n");
    for (int i = 0; i < num_rows; i++) {
      for (int j = 0; j < num_columns; j++)
        printf("%10e ", dense_matrix_out->get_value(i * num_columns + j));

      printf("\n");
    }

    printf("\n---------------------------------------------------------\n");
  }

  is_compressed_row = true;
  row_oriented = true;
  shared_ptr<SparseHbMatrix> m_csr_row_oriented =
      make_shared<SparseHbMatrix>(num_rows, num_columns, is_compressed_row);
  m_csr_row_oriented->copy_from_dense_matrix(
      dense_matrix_in, num_rows, num_columns, row_oriented, is_compressed_row);

  m_csr_row_oriented->get_dense_matrix(
      dense_matrix_out->get_non_const_values());

  if (TEST_EQUAL_DOUBLE_ARRAY(dense_matrix_in, dense_matrix_out->get_values(),
                              num_rows * num_columns,
                              "csr_row_oriented_matrix")) {
    printf("---------------------------------------------------------\n");
    printf("    Conversion between row-oriented dense matrix \n    "
           "and condensed row sparse matrix test passed!   \n");
    printf("---------------------------------------------------------\n");
  } else {

    printf("---------------------------------------------------------\n");
    printf("    Conversion between row-oriented dense matrix \n    "
           "and condensed row sparse matrix failed!   \n");
    printf("    Converting dense matrix to Sparse matrix FAILED!   \n");
    printf("\nMatrix in is \n");
    for (int i = 0; i < num_rows; i++) {
      for (int j = 0; j < num_columns; j++)
        printf("%10e ", dense_matrix_in[i * num_columns + j]);
      printf("\n");
    }

    printf("\nMatrix out is \n");
    for (int i = 0; i < num_rows; i++) {
      for (int j = 0; j < num_columns; j++)
        printf("%10e ", dense_matrix_out->get_value(i * num_columns + j));

      printf("\n");
    }

    printf("\n---------------------------------------------------------\n");
  }

  is_compressed_row = false;
  row_oriented = false;
  shared_ptr<SparseHbMatrix> m_csc_col_oriented =
      make_shared<SparseHbMatrix>(num_rows, num_columns, is_compressed_row);
  m_csc_col_oriented->copy_from_dense_matrix(
      dense_matrix_in, num_rows, num_columns, row_oriented, is_compressed_row);

  m_csc_col_oriented->get_dense_matrix(dense_matrix_out->get_non_const_values(),
                                       false);

  if (TEST_EQUAL_DOUBLE_ARRAY(dense_matrix_in, dense_matrix_out->get_values(),
                              num_rows * num_columns,
                              "csc_col_oriented_matrix")) {
    printf("---------------------------------------------------------\n");
    printf("    Conversion between col-oriented dense matrix \n    "
           "and condensed column sparse matrix test passed!   \n");
    printf("---------------------------------------------------------\n");
  } else {

    printf("---------------------------------------------------------\n");
    printf("    Conversion between col-oriented dense matrix \n    "
           "and condensed column sparse matrix failed!   \n");
    printf("    Converting dense matrix to Sparse matrix FAILED!   \n");
    printf("\nMatrix in is \n");
    for (int i = 0; i < num_rows; i++) {
      for (int j = 0; j < num_columns; j++)
        printf("%10e ", dense_matrix_in[i * num_columns + j]);
      printf("\n");
    }

    printf("\nMatrix out is \n");
    for (int i = 0; i < num_rows; i++) {
      for (int j = 0; j < num_columns; j++)
        printf("%10e ", dense_matrix_out->get_value(i * num_columns + j));

      printf("\n");
    }

    printf("\n---------------------------------------------------------\n");
  }

  is_compressed_row = true;
  row_oriented = false;
  shared_ptr<SparseHbMatrix> m_csr_col_oriented =
      make_shared<SparseHbMatrix>(num_rows, num_columns, is_compressed_row);
  m_csr_col_oriented->copy_from_dense_matrix(
      dense_matrix_in, num_rows, num_columns, row_oriented, is_compressed_row);

  m_csr_col_oriented->get_dense_matrix(dense_matrix_out->get_non_const_values(),
                                       false);

  if (TEST_EQUAL_DOUBLE_ARRAY(dense_matrix_in, dense_matrix_out->get_values(),
                              num_rows * num_columns,
                              "csr_col_oriented_matrix")) {
    printf("---------------------------------------------------------\n");
    printf("    Conversion between col-oriented dense matrix \n    "
           "and condensed row sparse matrix test passed!   \n");
    printf("---------------------------------------------------------\n");
  } else {

    printf("---------------------------------------------------------\n");
    printf("    Conversion between col-oriented dense matrix \n    "
           "and condensed row sparse matrix failed!   \n");
    printf("    Converting dense matrix to Sparse matrix FAILED!   \n");
    printf("\nMatrix in is \n");
    for (int i = 0; i < num_rows; i++) {
      for (int j = 0; j < num_columns; j++)
        printf("%10e ", dense_matrix_in[i * num_columns + j]);
      printf("\n");
    }

    printf("\nMatrix out is \n");
    for (int i = 0; i < num_rows; i++) {
      for (int j = 0; j < num_columns; j++)
        printf("%10e ", dense_matrix_out->get_value(i * num_columns + j));

      printf("\n");
    }

    printf("\n---------------------------------------------------------\n");
  }

  return true;
}

bool TEST_SPARSE_MATRIX_VECTOR_MULTIPLICATION(int num_rows, int num_columns,
                                              const double* dense_matrix_in,
                                              shared_ptr<const Vector> vector)
{

  assert(num_columns == vector->get_dim());

  bool is_csc_passed, is_csr_passed;

  shared_ptr<Vector> result_dense = make_shared<Vector>(num_rows);
  shared_ptr<Vector> result_sparse = make_shared<Vector>(num_rows);

  // compare it with matrix-vector multiplication in dense matrix
  result_dense->set_to_zero();
  for (int i = 0; i < num_rows; i++) {
    for (int j = 0; j < num_columns; j++) {
      result_dense->add_number_to_element(
          i, dense_matrix_in[i * num_columns + j] * vector->get_value(j));
    }
  }

  bool is_compressed_row = false;
  bool row_oriented = true;
  shared_ptr<SparseHbMatrix> m_csc =
      make_shared<SparseHbMatrix>(num_rows, num_columns, is_compressed_row);
  m_csc->copy_from_dense_matrix(dense_matrix_in, num_rows, num_columns,
                                row_oriented, is_compressed_row);

  result_sparse->set_to_zero();
  m_csc->multiply(vector, result_sparse);

  if (TEST_EQUAL_DOUBLE_ARRAY(result_sparse->get_values(),
                              result_dense->get_values(), num_rows)) {
    printf("---------------------------------------------------------\n");
    printf("   Matrix-vector multiplication for condensed column \n"
           "   sparse matrix passed!   \n");
    printf("---------------------------------------------------------\n");
    is_csc_passed = true;
  } else {
    printf("---------------------------------------------------------\n");
    printf("   Matrix-vector multiplication for condensed column\n"
           "   sparse matrix failed!   \n");
    result_sparse->print("result_sparse");
    result_dense->print("result_dense");
    printf("\n---------------------------------------------------------\n");
    is_csc_passed = false;
  }

  is_compressed_row = true;
  row_oriented = true;
  shared_ptr<SparseHbMatrix> m_csr =
      make_shared<SparseHbMatrix>(num_rows, num_columns, is_compressed_row);
  m_csr->copy_from_dense_matrix(dense_matrix_in, num_rows, num_columns,
                                row_oriented, is_compressed_row);

  result_sparse->set_to_zero();
  m_csr->multiply(vector, result_sparse);

  if (TEST_EQUAL_DOUBLE_ARRAY(result_sparse->get_values(),
                              result_dense->get_values(), num_rows)) {
    printf("---------------------------------------------------------\n");
    printf("   Matrix-vector multiplication for condensed row \n"
           "   sparse matrix passed!   \n");
    printf("---------------------------------------------------------\n");
    is_csr_passed = true;
  } else {
    printf("---------------------------------------------------------\n");
    printf("   Matrix-vector multiplication for condensed row\n"
           "   sparse matrix failed!   \n");
    result_sparse->print("result_sparse");
    result_dense->print("result_dense");
    printf("\n---------------------------------------------------------\n");
    is_csr_passed = false;
  }

  if (is_csc_passed && is_csr_passed)
    return true;
  else
    return false;
}

bool TEST_TRANSPOSED_MATRIX_VECTOR_MULTIPLICATION(
    int num_rows, int num_columns, const double* dense_matrix_in,
    shared_ptr<const Vector> vector)
{

  assert(num_rows == vector->get_dim());

  bool is_csc_passed, is_csr_passed;

  shared_ptr<Vector> result_dense = make_shared<Vector>(num_columns);
  shared_ptr<Vector> result_sparse = make_shared<Vector>(num_columns);

  // compare it with matrix-vector multiplication in dense matrix
  result_dense->set_to_zero();
  for (int i = 0; i < num_columns; i++) {
    for (int j = 0; j < num_rows; j++) {
      result_dense->add_number_to_element(
          i, dense_matrix_in[j * num_columns + i] * vector->get_value(j));
    }
  }

  bool is_compressed_row = false;
  bool row_oriented = true;
  shared_ptr<SparseHbMatrix> m_csc =
      make_shared<SparseHbMatrix>(num_rows, num_columns, is_compressed_row);
  m_csc->copy_from_dense_matrix(dense_matrix_in, num_rows, num_columns,
                                row_oriented, is_compressed_row);

  result_sparse->set_to_zero();
  m_csc->multiply_transpose(vector, result_sparse);

  if (TEST_EQUAL_DOUBLE_ARRAY(result_sparse->get_values(),
                              result_dense->get_values(), num_columns)) {
    printf("---------------------------------------------------------\n");
    printf("   Transposed matrix-vector multiplication for condensed \n"
           "   column sparse matrix passed!   \n");

    printf("---------------------------------------------------------\n");
    is_csc_passed = true;
  } else {
    printf("---------------------------------------------------------\n");
    printf("   Transposed matrix-vector multiplication for condensed \n"
           "   column sparse matrix failed!   \n");
    result_sparse->print("result_sparse");
    result_dense->print("result_dense");
    printf("\n---------------------------------------------------------\n");
    is_csc_passed = false;
  }

  is_compressed_row = true;
  row_oriented = true;
  shared_ptr<SparseHbMatrix> m_csr =
      make_shared<SparseHbMatrix>(num_rows, num_columns, is_compressed_row);
  m_csr->copy_from_dense_matrix(dense_matrix_in, num_rows, num_columns,
                                row_oriented, is_compressed_row);

  result_sparse->set_to_zero();
  m_csr->multiply_transpose(vector, result_sparse);

  if (TEST_EQUAL_DOUBLE_ARRAY(result_sparse->get_values(),
                              result_dense->get_values(), num_columns)) {
    printf("---------------------------------------------------------\n");
    printf("   Transposed matrix-vector multiplication for condensed \n"
           "   row sparse matrix passed!   \n");

    printf("---------------------------------------------------------\n");
    is_csr_passed = true;
  } else {
    printf("---------------------------------------------------------\n");
    printf("   Transposed matrix-vector multiplication for condensed \n"
           "   row sparse matrix failed!   \n");
    result_sparse->print("result_sparse");
    result_dense->print("result_dense");
    printf("\n---------------------------------------------------------\n");
    is_csr_passed = false;
  }

  if (is_csc_passed && is_csr_passed)
    return true;
  else
    return false;
}

bool TEST_TRIPLET_HB_MATIRX_CONVERSION(int num_rows, int num_columns,
                                       const double* dense_matrix_in)
{
  shared_ptr<Vector> dense_matrix_out =
      make_shared<Vector>(num_rows * num_columns);

  /**-------------------------------------------------------**/
  /**             Without Identity Submatrix                **/
  /**-------------------------------------------------------**/
  shared_ptr<SparseTripletMatrix> m_triplet_in =
      make_shared<SparseTripletMatrix>(dense_matrix_in, num_rows, num_columns,
                                       true); // condensed column sparse matrix
  bool is_compressed_row_format = false;
  auto m_csc = make_shared<SparseHbMatrix>(num_rows, num_columns,
                                           is_compressed_row_format);

  m_csc->set_structure(m_triplet_in);
  m_csc->set_values(m_triplet_in);

  auto m_triplet_csc_out = m_csc->convert_to_triplet();
  m_triplet_csc_out->get_dense_matrix(dense_matrix_out->get_non_const_values(),
                                      true);

  /**-------------------------------------------------------**/
  /**             Condensed Column matrix                   **/
  /**-------------------------------------------------------**/

  if (TEST_EQUAL_DOUBLE_ARRAY(dense_matrix_in, dense_matrix_out->get_values(),
                              m_triplet_in->get_num_entries())) {
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

  auto m_csr = make_shared<SparseHbMatrix>(num_rows, num_columns, true);
  m_csr->set_structure(m_triplet_in);
  m_csr->set_values(m_triplet_in);
  auto m_triplet_csr_out = m_csc->convert_to_triplet();
  m_triplet_csr_out->get_dense_matrix(dense_matrix_out->get_non_const_values(),
                                      true);

  if (TEST_EQUAL_DOUBLE_ARRAY(dense_matrix_in, dense_matrix_out->get_values(),
                              m_triplet_in->get_num_entries())) {
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

bool TEST_MATRIX_ALLOCATION_FROM_PERMUTATIONS()
{
}

bool TEST_SET_STRUCTURE(int num_rows, int num_columns,
                        const double* dense_matrix_in)
{
}

bool TEST_SET_MATRIX_VALUE()
{
}

int main(int argc, char* argv[])
{

  /**-------------------------------------------------------**/
  /**     Generate a random matrix for testing              **/
  /**-------------------------------------------------------**/
  srand(time(NULL));
  std::random_device rd;
  // std::mt19937 g(rd());

  int num_rows = rand() % 10 + 1;
  int num_columns = rand() % 10 + 1;
  int EntryNum = rand() % (num_rows * num_columns) + 1;

  std::vector<int> isNonzero;
  isNonzero.reserve(num_rows * num_columns);
  for (int i = 0; i < num_rows * num_columns; i++)
    isNonzero.push_back(i);

  // std::shuffle(isNonzero.begin(), isNonzero.end(),g);

  auto dense_matrix_in =
      new double[num_rows * num_columns](); // this constructor initializes all
                                            // values to 0.
  for (int i = 0; i < num_rows * num_columns; i++) {
    if (isNonzero.at(i) < EntryNum)
      // generate the entry between 1 to 10
      dense_matrix_in[i] = rand() % 10 + 1;
  }

  printf("\n=========================================================\n");
  printf("    Testing Methods for Harwell-Boeing Sparse Matrix\n"
         "   on randomly generated data.");
  printf("\n=========================================================\n");
  /**-------------------------------------------------------**/
  /**            Dense-sparse Matrix Conversion             **/
  /**-------------------------------------------------------**/

  TEST_DENSE_SPARSE_MATRIX_CONVERSION(num_rows, num_columns, dense_matrix_in);

  /**-------------------------------------------------------**/
  /**                 Matrix-vector Multiplication          **/
  /**-------------------------------------------------------**/

  shared_ptr<Vector> vector_mult = make_shared<Vector>(num_columns);

  for (int i = 0; i < num_columns; i++) {
    vector_mult->set_value(i, rand() % 10 + 1);
  }

  TEST_SPARSE_MATRIX_VECTOR_MULTIPLICATION(num_rows, num_columns,
                                           dense_matrix_in, vector_mult);

  /**-------------------------------------------------------**/
  /**      Transposed Matrix-vector Multiplication          **/
  /**-------------------------------------------------------**/

  shared_ptr<Vector> vector_trans_mult = make_shared<Vector>(num_rows);

  for (int i = 0; i < num_rows; i++) {
    vector_trans_mult->set_value(i, rand() % 10 + 1);
  }

  TEST_TRANSPOSED_MATRIX_VECTOR_MULTIPLICATION(
      num_rows, num_columns, dense_matrix_in, vector_trans_mult);

  /**-------------------------------------------------------**/
  /**            Triplet Hb Matrix Conversion               **/
  /**-------------------------------------------------------**/

  TEST_TRIPLET_HB_MATIRX_CONVERSION(num_rows, num_columns, dense_matrix_in);

  delete[] dense_matrix_in;
  isNonzero.clear();

  return 0;
}
