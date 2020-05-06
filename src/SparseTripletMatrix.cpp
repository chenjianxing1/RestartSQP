#include "restartsqp/SparseTripletMatrix.hpp"
#include "restartsqp/SparseHbMatrix.hpp"
#include <iostream>

using namespace std;
using namespace Ipopt;

namespace RestartSqp {
/** Constructor for an empty matrix with N non-zero
 * entries*/

SparseTripletMatrix::SparseTripletMatrix(int nnz, int num_rows, int num_columns,
                                         bool isSymmetric, bool allocate)
 : row_indices_(nullptr)
 , column_indices_(nullptr)
 , values_(nullptr)
 , is_allocated_(allocate)
 , is_symmetric_(isSymmetric)
{
  num_entries_ = nnz;
  num_rows_ = num_rows;
  num_columns_ = num_columns;
  // do nothing unless any data is to be assigned
  if (allocate) {
    row_indices_ = new int[nnz];
    column_indices_ = new int[nnz];
    values_ = new double[nnz];
  }
}

SparseTripletMatrix::SparseTripletMatrix(const SparseTripletMatrix& rhs)
  : is_allocated_(rhs.is_allocated_)
  , is_symmetric_(rhs.is_symmetric_)
  , num_columns_(rhs.num_columns_)
  , num_entries_(rhs.num_entries_)
  , num_rows_(rhs.num_rows_)
{
  assert(is_allocated_);

  // Reserve space for the arrays (TODO: To save memory we could store the structure in a separate object)
  values_ = new double[num_entries_];
  column_indices_ = new int[num_entries_];
  row_indices_ = new int[num_entries_];

  std::copy(rhs.values_, rhs.values_+num_entries_, values_);
  std::copy(rhs.column_indices_, rhs.column_indices_+num_entries_, column_indices_);
  std::copy(rhs.row_indices_, rhs.row_indices_+num_entries_, row_indices_);
}


SparseTripletMatrix::SparseTripletMatrix(const double* data, int num_rows,
                                         int num_columns, bool row_oriented)
 : row_indices_(nullptr)
 , column_indices_(nullptr)
 , values_(nullptr)
 , is_allocated_(true)
 , num_entries_(0)
 , num_rows_(num_rows)
 , num_columns_(num_columns)
{
  // tuple with nonzero entries specified by <row, col, entry, order>
  vector<tuple<int, int, double>> nonzero_entries;

  if (row_oriented) {
    // check if the input matrix is symmetric
    if (num_rows == num_columns) {
      is_symmetric_ = true;
      for (int i = 0; i < num_rows; i++) {
        for (int j = i + 1; j < num_columns; j++) {
          if (data[i * num_columns + j] != data[j * num_columns + i]) {
            is_symmetric_ = false;
            break;
          }
        }
      }
    } else
      is_symmetric_ = false;

    if (is_symmetric_) {
      for (int i = 0; i < num_rows; i++) {
        for (int j = i; j < num_columns; j++) {
          if (data[i * num_columns + j] != 0.) {
            nonzero_entries.push_back(make_tuple(
                i, j, data[i * num_columns + j]));
            num_entries_++;
          }
        }
      }
    } else {
      for (int i = 0; i < num_rows; i++)
        for (int j = 0; j < num_columns; j++) {
          if (data[i * num_columns + j] != 0.) {
            nonzero_entries.push_back(make_tuple(
                i, j, data[i * num_columns + j]));
            num_entries_++;
          }
        }
    }
  } else { // if the input dense matrix is column oriented
    // check if the input matrix is symmetric
    if (num_rows == num_columns) {
      is_symmetric_ = true;
      for (int i = 0; i < num_rows; i++) {
        for (int j = i + 1; j < num_columns; j++) {
          if (data[j * num_rows + i] != data[i * num_rows + j]) {
            is_symmetric_ = false;
            break;
          }
        }
      }
    } else
      is_symmetric_ = false;

    if (is_symmetric_) {
      for (int i = 0; i < num_rows; i++) {
        for (int j = i; j < num_columns; j++) {
          if (data[j * num_rows + i] != 0.) {
            nonzero_entries.push_back(
                make_tuple(i, j, data[j * num_rows + i]));
            num_entries_++;
          }
        }
      }
    } else {
      for (int i = 0; i < num_rows; i++)
        for (int j = 0; j < num_columns; j++) {
          if (data[j * num_rows + i] != 0.) {
            nonzero_entries.push_back(
                make_tuple(i, j, data[j * num_rows + i]));
            num_entries_++;
          }
        }
    }
  }
  // allocate memory for RowIndex_, ColIndex, MatVal, and order
  row_indices_ = new int[num_entries_];
  column_indices_ = new int[num_entries_];
  values_ = new double[num_entries_];

  for (int i = 0; i < num_entries_; i++) {
    row_indices_[i] = get<0>(nonzero_entries[i]);
    column_indices_[i] = get<1>(nonzero_entries[i]);
    values_[i] = get<2>(nonzero_entries[i]);
  }
  nonzero_entries.clear();
}

SparseTripletMatrix::SparseTripletMatrix(
    shared_ptr<const SparseHbMatrix> sparse_hb_matrix)
 : is_allocated_(true)
{
  num_entries_ = sparse_hb_matrix->get_num_entries();
  is_symmetric_ = sparse_hb_matrix->is_symmetric();
  num_rows_ = sparse_hb_matrix->get_num_rows();
  num_columns_ = sparse_hb_matrix->get_num_columns();

  row_indices_ = new int[num_entries_];
  column_indices_ = new int[num_entries_];
  values_ = new double[num_entries_];

  const int* hb_row_indices = sparse_hb_matrix->get_row_indices();
  const int* hb_column_indices = sparse_hb_matrix->get_column_indices();
  const double* hb_values = sparse_hb_matrix->get_values();

  if (sparse_hb_matrix->is_compressed_row_format()) {
    for (int row = 0; row < num_rows_; ++row) {
      int start = hb_row_indices[row];
      int end = hb_row_indices[row + 1];
      for (int i = start; i < end; ++i) {
        int col = hb_column_indices[i];
        row_indices_[i] = row;
        column_indices_[i] = col;
        values_[i] = hb_values[i];
      }
    }
  } else {
    for (int col = 0; col < num_columns_; ++col) {
      int start = hb_column_indices[col];
      int end = hb_column_indices[col + 1];
      for (int i = start; i < end; ++i) {
        int row = hb_row_indices[i];
        row_indices_[i] = row;
        column_indices_[i] = col;
        values_[i] = hb_values[i];
      }
    }
  }
}

SparseTripletMatrix::~SparseTripletMatrix()
{
  if (is_allocated_)
    free_memory_();
}

/**
 *@name print the sparse matrix in triplet form
 */
//@{
void SparseTripletMatrix::print(const string& matrix_name,
                                SmartPtr<Journalist> jnlst, EJournalLevel level,
                                EJournalCategory category) const
{
  if (IsValid(jnlst)) {
    jnlst->Printf(level, category, "Matrix %s with %d nonzero elements:\n",
                  matrix_name.c_str(), num_entries_);
    for (int i = 0; i < num_entries_; ++i) {
      jnlst->Printf(level, category, "%d %d %23.16e\n", row_indices_[i],
                    column_indices_[i]+1, values_[i]+1);
    }
  } else {
    printf("Matrix %s with %d nonzero elements:\n", matrix_name.c_str(),
           num_entries_);
    for (int i = 0; i < num_entries_; ++i) {
      printf("%d %d %23.16e\n", row_indices_[i]+1, column_indices_[i]+1,
             values_[i]);
    }
  }
}

/**
 *@name print the sparse matrix in dense form
 */
//@{
void SparseTripletMatrix::print_dense(const char* name,
                                      SmartPtr<Journalist> jnlst,
                                      EJournalLevel level,
                                      EJournalCategory category) const
{

  auto dense_matrix = new double[num_rows_ * num_columns_];

  for (int i = 0; i < num_entries_; i++) {
    dense_matrix[num_columns_ * row_indices_[i] + column_indices_[i]] = values_[i];
    if (is_symmetric_ && row_indices_[i] != column_indices_[i])
      dense_matrix[num_columns_ * column_indices_[i] + row_indices_[i]] = values_[i];
  }
  if (!IsNull(jnlst)) {
    //    if (name != nullptr) {
    //            jnlst->Printf(J_DBG,J_MATRIX,name);
    //            jnlst->Printf(J_DBG,J_MATRIX," =: {\n");
    //    }
    //    for (int i = 0; i < RowNum_; i++) {
    //        for (int j = 0; j < ColNum_; j++) {
    //            sprintf(mat_val, "%f  ", dense_matrix[i * ColNum() + j]);
    //               jnlst->Print(J_DBG,J_MATRIX,mat_val);
    //        }
    //           jnlst->Printf(J_DBG,J_MATRIX,"\n");
    //    }
    //       jnlst->Printf(J_DBG,J_MATRIX,"}\n\n");
  } else {
    if (name != nullptr)
      printf("%s =:{\n", name);

    for (int i = 0; i < num_rows_; i++) {
      for (int j = 0; j < num_columns_; j++) {
        printf("%10e  ", dense_matrix[i * get_num_columns() + j]);
      }
      printf("\n");
    }
    printf("}\n\n");
  }
  delete[] dense_matrix;
}
//@}

/** free all memory*/
void SparseTripletMatrix::free_memory_()
{
  if (is_allocated_) {
    delete[] row_indices_;
    row_indices_ = nullptr;
    delete[] column_indices_;
    column_indices_ = nullptr;
    delete[] values_;
    values_ = nullptr;
  }
}

/**
 * @brief Times a matrix with a vector p, the pointer to the matrix-vector
 * product  will be stored in the class member of another Vector class object
 * called "result"
 */

void SparseTripletMatrix::multiply(shared_ptr<const Vector> p,
                                   shared_ptr<Vector> result,
                                   double factor) const
{
  assert(p->get_dim() == num_columns_);
  assert(result->get_dim() == num_rows_);
  // For now we just implement this for factor = 1 and -1
  assert(factor == 1. || factor == -1.);

  if (is_symmetric_) {
    if (factor == 1.) {
      for (int i = 0; i < num_entries_; i++) {
        result->add_number_to_element(
            row_indices_[i],
            values_[i] * p->get_values()[column_indices_[i]]);
        if (row_indices_[i] != column_indices_[i]) {
          result->add_number_to_element(
              column_indices_[i],
              values_[i] * p->get_values()[row_indices_[i]]);
        }
      }
    } else {
      for (int i = 0; i < num_entries_; i++) {
        result->add_number_to_element(
            row_indices_[i],
            -values_[i] * p->get_values()[column_indices_[i]]);
        if (row_indices_[i] != column_indices_[i]) {
          result->add_number_to_element(
              column_indices_[i],
              -values_[i] * p->get_values()[row_indices_[i]]);
        }
      }
    }
  } else {
    if (factor == 1.) {
      for (int i = 0; i < num_entries_; i++) {
        result->add_number_to_element(
            row_indices_[i],
            values_[i] * p->get_values()[column_indices_[i]]);
      }
    } else {
      for (int i = 0; i < num_entries_; i++) {
        result->add_number_to_element(
            row_indices_[i],
            -values_[i] * p->get_values()[column_indices_[i]]);
      }
    }
  }
}

void SparseTripletMatrix::copy(shared_ptr<const SparseTripletMatrix> rhs,
                               bool deep_copy)
{

  if (!deep_copy) {
    row_indices_ = rhs->row_indices_;
    column_indices_ = rhs->column_indices_;
    values_ = rhs->values_;
  }
  for (int i = 0; i < num_entries_; i++) {
    row_indices_[i] = rhs->row_indices_[i];
    column_indices_[i] = rhs->column_indices_[i];
    values_[i] = rhs->values_[i];
  }
}

void SparseTripletMatrix::multiply_transpose(shared_ptr<const Vector> p,
                                             shared_ptr<Vector> result,
                                             double factor) const
{
  assert(p->get_dim() == num_rows_);
  assert(result->get_dim() == num_columns_);
  // For now we just implement this for factor = 1 and -1
  assert(factor == 1. || factor == -1.);

  if (is_symmetric_) {
    multiply(p, result);
  } else {
    if (factor == 1.) {
      for (int i = 0; i < num_entries_; i++) {
        result->add_number_to_element(column_indices_[i],
                                      values_[i] *
                                          p->get_values()[row_indices_[i]]);
      }
    } else {
      for (int i = 0; i < num_entries_; i++) {
        result->add_number_to_element(column_indices_[i],
                                      -values_[i] *
                                          p->get_values()[row_indices_[i]]);
      }
    }
  }
}

void SparseTripletMatrix::get_dense_matrix(double* dense_matrix,
                                           bool row_oriented) const
{
  for (int i = 0; i < num_rows_ * num_columns_; i++) {
    dense_matrix[i] = 0.;
  }
  if (row_oriented) {
    for (int i = 0; i < num_entries_; i++) {
      dense_matrix[num_columns_ * row_indices_[i] +
                   column_indices_[i]] = values_[i];
      if (is_symmetric_ && row_indices_[i] != column_indices_[i])
        dense_matrix[num_columns_ * column_indices_[i] + row_indices_[i]] = values_[i];
    }
  } else {
    for (int i = 0; i < num_entries_; i++) {
      dense_matrix[num_rows_ * column_indices_[i] + row_indices_[i]] =
          values_[i];
      if (is_symmetric_ && row_indices_[i] != column_indices_[i])
        dense_matrix[num_rows_ * row_indices_[i] + column_indices_[i]] = values_[i];
    }
  }
}

/**
 * qpOASESSparseMatrix
 */
}
