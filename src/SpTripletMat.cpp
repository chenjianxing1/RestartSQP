
#include "sqphot/SpTripletMat.hpp"
#include <iostream>

namespace SQPhotstart {
/** Constructor for an empty matrix with N non-zero
 * entries*/

SpTripletMat::SpTripletMat(int nnz, int num_rows, int num_columns,
                           bool isSymmetric, bool allocate)
 : row_indices_(nullptr)
 , column_indices_(nullptr)
 , values_(nullptr)
 , order_(nullptr)
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
    order_ = new int[nnz];
    // initialize the order to 0:N-1
    for (int i = 0; i < nnz; i++) {
      order_[i] = i;
    }
  }
}

SpTripletMat::SpTripletMat(const double* data, int num_rows, int num_columns,
                           bool row_oriented)
 : row_indices_(nullptr)
 , column_indices_(nullptr)
 , values_(nullptr)
 , order_(nullptr)
 , is_allocated_(true)
 , num_entries_(0)
 , num_rows_(num_rows)
 , num_columns_(num_columns)
{
  // tuple with nonzero entries specified by <row, col, entry, order>
  std::vector<std::tuple<int, int, double, int>> nonzero_entries;

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
          if (abs(data[i * num_columns + j]) > m_eps) {
            nonzero_entries.push_back(std::make_tuple(
                i + 1, j + 1, data[i * num_columns + j], num_entries_));
            num_entries_++;
          }
        }
      }
    } else {
      for (int i = 0; i < num_rows; i++)
        for (int j = 0; j < num_columns; j++) {
          if (abs(data[i * num_columns + j]) > m_eps) {
            nonzero_entries.push_back(std::make_tuple(
                i + 1, j + 1, data[i * num_columns + j], num_entries_));
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
          if (abs(data[j * num_rows + i]) > m_eps) {
            nonzero_entries.push_back(std::make_tuple(
                i + 1, j + 1, data[j * num_rows + i], num_entries_));
            num_entries_++;
          }
        }
      }
    } else {
      for (int i = 0; i < num_rows; i++)
        for (int j = 0; j < num_columns; j++) {
          if (abs(data[j * num_rows + i]) > m_eps) {
            nonzero_entries.push_back(std::make_tuple(
                i + 1, j + 1, data[j * num_rows + i], num_entries_));
            num_entries_++;
          }
        }
    }
  }
  // allocate memory for RowIndex_, ColIndex, MatVal, and order
  row_indices_ = new int[num_entries_];
  column_indices_ = new int[num_entries_];
  values_ = new double[num_entries_];
  order_ = new int[num_entries_];

  for (int i = 0; i < num_entries_; i++) {
    row_indices_[i] = std::get<0>(nonzero_entries[i]);
    column_indices_[i] = std::get<1>(nonzero_entries[i]);
    values_[i] = std::get<2>(nonzero_entries[i]);
    order_[i] = std::get<3>(nonzero_entries[i]);
  }
  nonzero_entries.clear();
}

/** Default destructor */
SpTripletMat::~SpTripletMat()
{
  if (is_allocated_)
    free_memory_();
}

/**
 *@name print the sparse matrix in triplet form
 */
//@{
void SpTripletMat::print(const char* name,
                         Ipopt::SmartPtr<Ipopt::Journalist> jnlst,
                         Ipopt::EJournalLevel level,
                         Ipopt::EJournalCategory category) const
{

  std::cout << "Row    Column    Entry    Order" << std::endl;
  for (int i = 0; i < num_entries_; i++) {
    printf("%d       ", row_indices_[i]);
    printf("%d       ", column_indices_[i]);
    printf("%10e     ", values_[i]);
    printf("%d     \n", order_[i]);
  }
}

/**
 *@name print the sparse matrix in dense form
 */
//@{
void SpTripletMat::print_full(const char* name,
                              Ipopt::SmartPtr<Ipopt::Journalist> jnlst,
                              Ipopt::EJournalLevel level,
                              Ipopt::EJournalCategory category) const
{

  auto dense_matrix = new double[num_rows_ * num_columns_];

  for (int i = 0; i < num_entries_; i++) {
    dense_matrix[num_columns_ * (row_indices_[i] - 1) + column_indices_[i] -
                 1] = values_[i];
    if (is_symmetric_ && row_indices_[i] != column_indices_[i])
      dense_matrix[num_columns_ * (column_indices_[i] - 1) + row_indices_[i] -
                   1] = values_[i];
  }
  if (!IsNull(jnlst)) {
    //    if (name != nullptr) {
    //            jnlst->Printf(Ipopt::J_DBG,Ipopt::J_MATRIX,name);
    //            jnlst->Printf(Ipopt::J_DBG,Ipopt::J_MATRIX," =: {\n");
    //    }
    //    for (int i = 0; i < RowNum_; i++) {
    //        for (int j = 0; j < ColNum_; j++) {
    //            sprintf(mat_val, "%f  ", dense_matrix[i * ColNum() + j]);
    //               jnlst->Print(Ipopt::J_DBG,Ipopt::J_MATRIX,mat_val);
    //        }
    //           jnlst->Printf(Ipopt::J_DBG,Ipopt::J_MATRIX,"\n");
    //    }
    //       jnlst->Printf(Ipopt::J_DBG,Ipopt::J_MATRIX,"}\n\n");
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
void SpTripletMat::free_memory_()
{
  if (is_allocated_) {
    delete[] row_indices_;
    row_indices_ = nullptr;
    delete[] column_indices_;
    column_indices_ = nullptr;
    delete[] values_;
    values_ = nullptr;
    delete[] order_;
    order_ = nullptr;
  }
}

/**
 * @brief Times a matrix with a vector p, the pointer to the matrix-vector
 * product  will be stored in the class member of another Vector class object
 * called "result"
 */

void SpTripletMat::multiply(std::shared_ptr<const Vector> p,
                            std::shared_ptr<Vector> result) const
{

  result->set_to_zero();
  assert(num_columns_ == p->get_dim());
  if (is_symmetric_) {
    result->set_to_zero();
    for (int i = 0; i < num_entries_; i++) {
      result->add_number_to_element(
          row_indices_[i] - 1,
          values_[i] * p->get_values()[column_indices_[i] - 1]);
      if (row_indices_[i] != column_indices_[i]) {
        result->add_number_to_element(column_indices_[i] - 1,
                                      values_[i] *
                                          p->get_values()[row_indices_[i] - 1]);
      }
    }
  } else {
    result->set_to_zero(); // set all entries to be 0
    for (int i = 0; i < num_entries_; i++) {
      result->add_number_to_element(
          row_indices_[i] - 1,
          values_[i] * p->get_values()[column_indices_[i] - 1]);
    }
  }
}
#if 0
double SpTripletMat::calc_inf_norm()
{
  // TODO: test it!

  std::shared_ptr<Vector> rowSums = std::make_shared<Vector>(num_rows_);
  for (int i = 0; i < num_entries_; i++) {
    rowSums->add_number_to_element(get_row_indices()[i] - 1, abs(values_[i]));
    if (is_symmetric_)
      rowSums->add_number_to_element(get_column_indices()[i] - 1, abs(values_[i]));
  }

  double InfNorm =
      rowSums->calc_inf_norm(); // same as calculating the MAX of an array

  return InfNorm;
}

double SpTripletMat::calc_one_norm()
{

  std::shared_ptr<Vector> colSums = std::make_shared<Vector>(num_columns_);
  for (int i = 0; i < num_entries_; i++) {
    colSums->add_number_to_element(get_column_indices()[i] - 1, abs(values_[i]));
    if (is_symmetric_)
      colSums->add_number_to_element(get_row_indices()[i] - 1, abs(values_[i]));
  }

  double OneNorm =
      colSums->calc_inf_norm(); // same as calculating the MAX of an array

  return OneNorm;
}
#endif
void SpTripletMat::copy(std::shared_ptr<const SpTripletMat> rhs, bool deep_copy)
{

  if (!deep_copy) {
    row_indices_ = rhs->row_indices_;
    column_indices_ = rhs->column_indices_;
    values_ = rhs->values_;
    order_ = rhs->order_;
  }
  for (int i = 0; i < num_entries_; i++) {
    row_indices_[i] = rhs->row_indices_[i];
    column_indices_[i] = rhs->column_indices_[i];
    values_[i] = rhs->values_[i];
    order_[i] = rhs->order_[i];
  }
}

void SpTripletMat::multiply_transpose(std::shared_ptr<const Vector> p,
                                      std::shared_ptr<Vector> result) const
{

  if (is_symmetric_) {
    multiply(p, result);
  } else {
    result->set_to_zero(); // set all entries to be 0
    for (int i = 0; i < num_entries_; i++) {
      result->add_number_to_element(column_indices_[i] - 1,
                                    values_[i] *
                                        p->get_values()[row_indices_[i] - 1]);
    }
  }
}

void SpTripletMat::get_dense_matrix(double* dense_matrix,
                                    bool row_oriented) const
{
  for (int i = 0; i < num_rows_ * num_columns_; i++) {
    dense_matrix[i] = 0.;
  }
  if (row_oriented) {
    for (int i = 0; i < num_entries_; i++) {
      dense_matrix[num_columns_ * (row_indices_[i] - 1) +
                   (column_indices_[i] - 1)] = values_[i];
      if (is_symmetric_ && row_indices_[i] != column_indices_[i])
        dense_matrix[num_columns_ * (column_indices_[i] - 1) + row_indices_[i] -
                     1] = values_[i];
    }
  } else {
    for (int i = 0; i < num_entries_; i++) {
      dense_matrix[num_rows_ * (column_indices_[i] - 1) + row_indices_[i] - 1] =
          values_[i];
      if (is_symmetric_ && row_indices_[i] != column_indices_[i])
        dense_matrix[num_rows_ * (row_indices_[i] - 1) + column_indices_[i] -
                     1] = values_[i];
    }
  }
}

/**
 * qpOASESSparseMatrix
 */
}
