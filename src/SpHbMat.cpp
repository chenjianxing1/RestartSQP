
#include <algorithm>
#include <iostream>

#include "sqphot/SpHbMat.hpp"

using namespace std;

namespace SQPhotstart {

/** @name constructor/destructor */
//@{
/** Default constructor*/
SpHbMat::SpHbMat(int RowNum, int ColNum, bool isCompressedRow)
 : row_indices_(NULL)
 , column_indices_(NULL)
 , values_(NULL)
 , order_(NULL)
 , num_entries_(-1)
 , is_initialized_(false)
 , num_rows_(RowNum)
 , num_columns_(ColNum)
 , is_compressed_row_format_(isCompressedRow)
{
  if (is_compressed_row_format_) {
    row_indices_ = new int[RowNum + 1]();
  } else
    column_indices_ = new int[ColNum + 1]();
}

/**
 *@brief A constructor with the number of non-zero entries, row number and
 *column
 * number specified.
 *
 * @param: nnz number of nonzero entry
 * @param RowNum: number of rows of a matrix
 * @param ColNum: number of columns of a matrix
 */
SpHbMat::SpHbMat(int nnz, int RowNum, int ColNum, bool isCompressedRow)
 : row_indices_(NULL)
 , column_indices_(NULL)
 , values_(NULL)
 , order_(NULL)
 , num_entries_(nnz)
 , num_rows_(RowNum)
 , num_columns_(ColNum)
 , is_initialized_(false)
 , is_compressed_row_format_(isCompressedRow)
{

  if (isCompressedRow) {
    column_indices_ = new int[nnz]();
    row_indices_ = new int[RowNum + 1]();
  } else {
    column_indices_ = new int[ColNum + 1]();
    row_indices_ = new int[nnz]();
  }
  values_ = new double[nnz]();
  order_ = new int[nnz]();
  for (int i = 0; i < nnz; i++)
    order_[i] = i;
}

SpHbMat::SpHbMat(const double* data, int RowNum, int ColNum, bool row_oriented,
                 bool isCompressedRow)
 : row_indices_(NULL)
 , column_indices_(NULL)
 , values_(NULL)
 , order_(NULL)
 , num_entries_(0)
 , num_rows_(RowNum)
 , num_columns_(ColNum)
 , is_compressed_row_format_(isCompressedRow)
{

  int* RowIndex_tmp = NULL;
  int* ColIndex_tmp = NULL;
  double* MatVal_tmp = NULL;
  // allocate memory
  if (isCompressedRow) {
    row_indices_ = new int[RowNum + 1]();
    ColIndex_tmp = new int[RowNum * ColNum]();
    MatVal_tmp = new double[RowNum * ColNum]();
  } else {
    column_indices_ = new int[ColNum + 1]();
    RowIndex_tmp = new int[RowNum * ColNum]();
    MatVal_tmp = new double[RowNum * ColNum]();
  }

  if (row_oriented) {
    if (isCompressedRow) {
      for (int i = 0; i < RowNum; i++) {
        for (int j = 0; j < ColNum; j++) {
          // identify nonzero entry
          if (abs(data[i * ColNum + j]) > m_eps) {
            MatVal_tmp[num_entries_] = data[i * ColNum + j];
            ColIndex_tmp[num_entries_] = j;
            num_entries_++;
          }
        }
        row_indices_[i + 1] = num_entries_;
      }
    } else { // if it is condensed column
      for (int j = 0; j < ColNum; j++) {
        for (int i = 0; i < RowNum; i++) {
          if (abs(data[j + i * ColNum]) > m_eps) {
            MatVal_tmp[num_entries_] = data[j + i * ColNum];
            RowIndex_tmp[num_entries_] = i;
            num_entries_++;
          }
        }
        column_indices_[j + 1] = num_entries_;
      }
    }
  } else {
    if (isCompressedRow) {
      for (int j = 0; j < RowNum; j++) {
        for (int i = 0; i < ColNum; i++) {
          if (abs(data[j + i * RowNum]) > m_eps) {
            MatVal_tmp[num_entries_] = data[j + i * RowNum];
            ColIndex_tmp[num_entries_] = i;
            num_entries_++;
          }
        }
        row_indices_[j + 1] = num_entries_;
      }
    } else {
      for (int i = 0; i < ColNum; i++) {
        for (int j = 0; j < RowNum; j++) {
          // identify nonzero entry
          if (abs(data[i * RowNum + j]) > m_eps) {
            MatVal_tmp[num_entries_] = data[i * RowNum + j];
            RowIndex_tmp[num_entries_] = j;
            num_entries_++;
          }
        }
        column_indices_[i + 1] = num_entries_;
      }
    }
  }
  // allocate memory for class member
  values_ = new double[num_entries_];
  if (isCompressedRow) {
    column_indices_ = new int[num_entries_];
    for (int i = 0; i < num_entries_; i++) {
      column_indices_[i] = ColIndex_tmp[i];
      values_[i] = MatVal_tmp[i];
    }
    delete[] MatVal_tmp;
    MatVal_tmp = NULL;
    delete[] ColIndex_tmp;
    ColIndex_tmp = NULL;

  } else {
    row_indices_ = new int[num_entries_];
    for (int i = 0; i < num_entries_; i++) {
      row_indices_[i] = RowIndex_tmp[i];
      values_[i] = MatVal_tmp[i];
    }
    delete[] MatVal_tmp;
    MatVal_tmp = NULL;
    delete[] RowIndex_tmp;
    RowIndex_tmp = NULL;
  }
  order_ = new int[num_entries_];
}

/**
 *Default destructor
 */
SpHbMat::~SpHbMat()
{

  freeMemory();
}

//@}

/** @name setStructure */

//@{

/**
 * @brief setup the structure of the sparse matrix for QPsolvrs
 * This method should be only called for once
 *
 * This method will convert the strucutre information from the triplet form from
 * a
 * SpMatrix object to the format required by the corresponding QPsolvers
 *
 * @param rhs a SpMatrix object whose content will be copied to the class
 * members
 * (in a different sparse matrix representations)
 * @param I_info the information of 2 identity sub matrices.
 *
 */
void SpHbMat::setStructure(std::shared_ptr<const SpTripletMat> rhs,
                           IdentityInfo I_info)
{
  set_zero();
  is_symmetric_ = false;
  assert(is_initialized_ == false);

  int counter = 0; // the counter for recording the index location
  std::vector<std::tuple<int, int, double, int>> sorted_index_info;
  for (int i = 0; i < rhs->get_num_entries(); i++) {
    sorted_index_info.push_back(std::make_tuple(
        rhs->get_row_index_at_entry(i), rhs->get_column_index_at_entry(i),
        rhs->get_value_at_entry(i), counter));
    counter++;
  }

  // adding identity matrices info to the tuple array.
  for (int i = 0; i < I_info.length; i++) {
    for (int j = 0; j < I_info.size[i]; j++) {
      sorted_index_info.push_back(std::make_tuple(
          I_info.irow[i] + j, I_info.jcol[i] + j, I_info.value[i], counter));
      counter++;
    }
  }

  assert(counter == num_entries_);

  //  for(int i = 0; i <sorted_index_info.size(); i++) {
  //      printf("%i ", std::get<0>(sorted_index_info[i]));
  //      printf("%i ", std::get<1>(sorted_index_info[i]));
  //      printf("%5.2e ", std::get<2>(sorted_index_info[i]));
  //      printf("%i \n", std::get<3>(sorted_index_info[i]));

  //  }

  if (is_compressed_row_format_) {
    std::sort(sorted_index_info.begin(), sorted_index_info.end(),
              tuple_sort_rule_compressed_row);
    for (int i = 0; i < num_entries_; i++) {
      column_indices_[i] = std::get<1>(sorted_index_info[i]) - 1;
      values_[i] = std::get<2>(sorted_index_info[i]);
      order_[std::get<3>(sorted_index_info[i])] = i;

      for (int j = std::get<0>(sorted_index_info[i]); j < num_rows_; j++) {
        row_indices_[j]++;
      }
    }
    row_indices_[num_rows_] = num_entries_;

  } else {
    std::sort(sorted_index_info.begin(), sorted_index_info.end(),
              tuple_sort_rule_compressed_column);
    for (int i = 0; i < num_entries_; i++) {
      values_[i] = std::get<2>(sorted_index_info[i]);
      row_indices_[i] = std::get<0>(sorted_index_info[i]) - 1;
      order_[std::get<3>(sorted_index_info[i])] = i;

      for (int j = std::get<1>(sorted_index_info[i]); j < num_columns_; j++) {
        column_indices_[j]++;
      }
    }
    column_indices_[num_columns_] = num_entries_;
  }
  is_initialized_ = true;
  sorted_index_info.clear();
}

/**
 * @brief setup the structure of the sparse matrix for solver qpOASES(should
 * be called only for once).
 *
 * This method will convert the strucutre information from the triplet form from
 * a
 * SpMatrix object to the format required by the QPsolver qpOASES.
 *
 * @param rhs a SpMatrix object whose content will be copied to the class
 * members
 * (in a different sparse matrix representations)
 *
 */
void SpHbMat::setStructure(std::shared_ptr<const SpTripletMat> rhs)
{

  set_zero();
  is_symmetric_ = rhs->is_symmetric();

  assert(is_initialized_ == false);
  std::vector<std::tuple<int, int, double, int>> sorted_index_info;

  // if it is symmetric, it will calculate the
  // number of entry and allocate_memory the memory
  // of RowIndex and MatVal by going through
  // all of its entries one by one
  for (int i = 0; i < rhs->get_num_entries(); i++) {
    sorted_index_info.emplace_back(
        rhs->get_row_index_at_entry(i), rhs->get_column_index_at_entry(i),
        rhs->get_value_at_entry(i), sorted_index_info.size());

    if (is_symmetric_ &&
        rhs->get_row_index_at_entry(i) != rhs->get_column_index_at_entry(i)) {
      sorted_index_info.emplace_back(
          rhs->get_column_index_at_entry(i), rhs->get_row_index_at_entry(i),
          rhs->get_value_at_entry(i), sorted_index_info.size());
    }
  }

  if (num_entries_ < 0) { // no memory allocated so far
    num_entries_ = sorted_index_info.size();
    values_ = new double[num_entries_]();
    order_ = new int[num_entries_]();

    if (is_compressed_row_format_) {
      column_indices_ = new int[num_entries_]();
    } else {
      row_indices_ = new int[num_entries_]();
    }
  }

  if (is_compressed_row_format_) {
    std::sort(sorted_index_info.begin(), sorted_index_info.end(),
              tuple_sort_rule_compressed_row);
    // copy the order information back
    for (int i = 0; i < num_entries_; i++) {
      column_indices_[i] = std::get<1>(sorted_index_info[i]) - 1;
      values_[i] = std::get<2>(sorted_index_info[i]);
      order_[std::get<3>(sorted_index_info[i])] = i;
      for (int j = std::get<0>(sorted_index_info[i]); j < num_rows_; j++) {
        row_indices_[j]++;
      }
    }
    row_indices_[num_rows_] = num_entries_;
  } else {
    std::sort(sorted_index_info.begin(), sorted_index_info.end(),
              tuple_sort_rule_compressed_column);
    // copy the order information back

    for (int i = 0; i < num_entries_; i++) {
      row_indices_[i] = std::get<0>(sorted_index_info[i]) - 1;
      values_[i] = std::get<2>(sorted_index_info[i]);
      order_[std::get<3>(sorted_index_info[i])] = i;

      for (int j = std::get<1>(sorted_index_info[i]); j < num_columns_; j++) {
        column_indices_[j]++;
      }
    }
    column_indices_[num_columns_] = num_entries_;
  }
  is_initialized_ = true;
  sorted_index_info.clear();
}

//@}

/** @name setMatVal */
//@{
/**
 * @brief set the Matrix values to the matrix, convert from triplet format to
 * Harwell-Boeing Matrix format.
 * @param rhs entry values(orders are not yet under permutation)
 * @param I_info struct which stores identity matrices information
 */
void SpHbMat::setMatVal(std::shared_ptr<const SpTripletMat> rhs,
                        IdentityInfo I_info)
{
  // adding identity submatrices  to the matrix
  int total_I_entries = 0;
  for (int i = 0; i < I_info.length; i++) {
    total_I_entries += I_info.size[i];
  }

  // assign each matrix entry to the corresponding position after permutation
  for (int i = 0; i < num_entries_ - total_I_entries; i++) {
    values_[order_[i]] = rhs->get_value_at_entry(i);
  }
}

void SpHbMat::setMatVal(std::shared_ptr<const SpTripletMat> rhs)
{
  int j = 0;
  for (int i = 0; i < rhs->get_num_entries(); i++) {
    values_[order_[j]] = rhs->get_value_at_entry(i);
    j++;
    if (is_symmetric_ &&
        (rhs->get_column_index_at_entry(i) != rhs->get_row_index_at_entry(i))) {
      values_[order_[j]] = rhs->get_value_at_entry(i);
      j++;
    }
  }
}

//@}

/**
 * Free all memory allocated
 */
void SpHbMat::freeMemory()
{
  delete[] column_indices_;
  column_indices_ = NULL;
  delete[] row_indices_;
  row_indices_ = NULL;
  delete[] values_;
  values_ = NULL;
  delete[] order_;
  order_ = NULL;
}

void SpHbMat::copy(std::shared_ptr<const SpHbMat> rhs)
{

  assert(num_entries_ == rhs->get_num_entries());
  assert(num_rows_ == rhs->get_num_rows());
  assert(num_columns_ == rhs->get_num_columns());
  for (int i = 0; i < num_entries_; i++) {
    row_indices_[i] = rhs->get_row_index_at_entry(i);
    values_[i] = rhs->get_value_at_entry(i);
    order_[i] = rhs->get_order_at_entry(i);
  }

  for (int i = 0; i < num_columns_ + 1; i++) {
    column_indices_[i] = rhs->get_column_index_at_entry(i);
  }
}

void SpHbMat::set_zero()
{
  if (is_compressed_row_format_) {
    for (int i = 0; i < num_entries_; i++) {
      values_[i] = 0;
      column_indices_[i] = 0;
      order_[i] = i;
    }
    for (int i = 0; i < num_rows_ + 1; i++) {
      row_indices_[i] = 0;
    }
  } else {
    for (int i = 0; i < num_entries_; i++) {
      values_[i] = 0;
      row_indices_[i] = 0;
      order_[i] = i;
    }
    for (int i = 0; i < num_columns_ + 1; i++) {
      column_indices_[i] = 0;
    }
  }
}

shared_ptr<SpTripletMat> SpHbMat::convert_to_triplet() const
{

  shared_ptr<SpTripletMat> result;
  int nnz; // number of non-zero entries
  int j = 1;
  if (is_symmetric_) {
    std::vector<std::tuple<int, int, double, int>> triplet_vector;
    for (int i = 0; i < num_entries_; i++) {
      if (is_compressed_row_format_) {
        while (row_indices_[j] == i)
          j++;
        triplet_vector.emplace_back(j, column_indices_[i] + 1, values_[i],
                                    order_[i]);
      }
    }
    for (int i = 0; i < num_entries_; i++) // delete repetitive entry
      if (get<0>(triplet_vector[i]) != get<1>(triplet_vector[i]))
        for (int k = i; k < triplet_vector.size(); k++) {
          if (get<0>(triplet_vector[i]) == get<1>(triplet_vector[k]) &&
              get<1>(triplet_vector[i]) == get<0>(triplet_vector[k]))
            triplet_vector.erase(triplet_vector.begin() + k);
        }

    if (is_compressed_row_format_)
      std::sort(triplet_vector.begin(), triplet_vector.end(),
                tuple_sort_rule_compressed_row);
    else
      std::sort(triplet_vector.begin(), triplet_vector.end(),
                tuple_sort_rule_compressed_column);

    result = make_shared<SpTripletMat>(triplet_vector.size(), num_rows_,
                                       num_columns_, true, true);

    for (int i = 0; i < result->get_num_entries(); i++) {
      result->set_row_index_at_entry(i, get<0>(triplet_vector[i]));
      result->set_column_index_at_entry(i, get<1>(triplet_vector[i]));
      result->set_value_at_entry(i, get<2>(triplet_vector[i]));
      //            result->setOrderAt(i,get<3>(triplet_vector[i]));
    }
  } else {
    result = make_shared<SpTripletMat>(num_entries_, num_rows_, num_columns_,
                                       false, true);
    for (int i = 0; i < num_entries_; i++) {
      if (is_compressed_row_format_) {
        while (row_indices_[j] == i)
          j++;
        result->set_row_index_at_entry(i, j);
        result->set_column_index_at_entry(i, column_indices_[i] + 1);
      } else {
        while (column_indices_[j] == i)
          j++;
        result->set_column_index_at_entry(i, j);
        result->set_row_index_at_entry(i, row_indices_[i] + 1);
      }

      result->set_value_at_entry(i, values_[i]);
    }
  }
  return result;
}

void SpHbMat::write_to_file(const char* name,
                            Ipopt::SmartPtr<Ipopt::Journalist> jnlst,
                            Ipopt::EJournalLevel level,
                            Ipopt::EJournalCategory category, Solver solver)
{
#ifdef DEBUG
#ifdef PRINT_QP_IN_CPP
  const char* var_type_int;
  const char* var_type_double;
  var_type_int = (solver == QPOASES) ? "sparse_int_t" : "qp_int";
  var_type_double = (solver == QPOASES) ? "real_t" : "double";
  jnlst->Printf(level, category, "%s %s_jc[] = \n{", var_type_int, name);
  int i;
  for (i = 0; i < ColNum_ + 1; i++) {
    if (i % 10 == 0 && i > 1)
      jnlst->Printf(level, category, "\n");
    if (i == ColNum_)
      jnlst->Printf(level, category, "%d};\n\n", ColIndex_[i]);
    else
      jnlst->Printf(level, category, "%d, ", ColIndex_[i]);
  }
  jnlst->Printf(level, category, "%s %s_ir[] = \n{", var_type_int, name);
  for (i = 0; i < EntryNum_; i++) {
    if (i % 10 == 0 && i > 1)
      jnlst->Printf(level, category, "\n");
    if (i == EntryNum_ - 1)
      jnlst->Printf(level, category, "%d};\n\n", RowIndex_[i]);
    else
      jnlst->Printf(level, category, "%d, ", RowIndex_[i]);
  }
  jnlst->Printf(level, category, "%s %s_val[] = \n{", var_type_double, name);
  for (i = 0; i < EntryNum_; i++) {
    if (i % 10 == 0 && i > 1)
      jnlst->Printf(level, category, "\n");
    if (i == EntryNum_ - 1)
      jnlst->Printf(level, category, "%23.16e};\n\n", MatVal_[i]);
    else
      jnlst->Printf(level, category, "%23.16e, ", MatVal_[i]);
  }
#else
  int i;
  if (solver == QORE) {
    for (i = 0; i < RowNum_ + 1; i++) {
      jnlst->Printf(level, category, "%d\n", RowIndex_[i]);
    }
    for (i = 0; i < EntryNum_; i++) {
      jnlst->Printf(level, category, "%d\n", ColIndex_[i]);
    }
    for (i = 0; i < EntryNum_; i++) {
      jnlst->Printf(level, category, "%23.16e\n", MatVal_[i]);
    }
  } else if (solver == QPOASES) {
    for (i = 0; i < EntryNum_; i++) {
      jnlst->Printf(level, category, "%d\n", RowIndex_[i]);
    }
    for (i = 0; i < ColNum_ + 1; i++) {
      jnlst->Printf(level, category, "%d\n", ColIndex_[i]);
    }
    for (i = 0; i < EntryNum_; i++) {
      jnlst->Printf(level, category, "%23.16e\n", MatVal_[i]);
    }
  }
#endif
#endif
}

void SpHbMat::get_dense_matrix(double* dense_matrix, bool row_oriented) const
{

  // Initialize dense matrix values to zero
  for (int i = 0; i < num_rows_ * num_columns_; i++) {
    dense_matrix[i] = 0.;
  }
  int row;
  if (row_oriented) {
    if (is_compressed_row_format_) {
      for (int i = 1; i < num_rows_ + 1; i++) {
        if (row_indices_[i] > 0) {
          row = i - 1;
          break;
        }
      }
      for (int i = 0; i < num_entries_; i++) {
        while (i == row_indices_[row + 1]) {
          row++;
        }
        dense_matrix[num_columns_ * row + column_indices_[i]] = values_[i];
      }
    } else {
      int col;
      for (int i = 1; i < num_columns_ + 1; i++) {
        if (column_indices_[i] > 0) {
          col = i - 1;
          break;
        }
      }
      for (int i = 0; i < num_entries_; i++) {
        while (i == column_indices_[col + 1]) {
          col++;
        }

        dense_matrix[num_columns_ * row_indices_[i] + col] = values_[i];
      }
    }
  } else {
    if (is_compressed_row_format_) {
      for (int i = 1; i < num_rows_ + 1; i++) {
        if (row_indices_[i] > 0) {
          row = i - 1;
          break;
        }
      }
      for (int i = 0; i < num_entries_; i++) {
        while (i == row_indices_[row + 1]) {
          row++;
        }
        dense_matrix[num_rows_ * column_indices_[i] + row] = values_[i];
      }
    } else {
      int col;
      for (int i = 1; i < num_columns_ + 1; i++) {
        if (column_indices_[i] > 0) {
          col = i - 1;
          break;
        }
      }
      for (int i = 0; i < num_entries_; i++) {
        while (i == column_indices_[col + 1]) {
          col++;
        }

        dense_matrix[num_rows_ * col + row_indices_[i]] = values_[i];
      }
    }
  }
}

void SpHbMat::multiply_transpose(shared_ptr<const Vector> p,
                                 shared_ptr<Vector> result) const
{

  result->set_to_zero();

  if (is_compressed_row_format_) {
    int row;
    for (int i = 1; i < num_rows_ + 1; i++) {
      if (row_indices_[i] > 0) {
        row = i - 1;
        break;
      }
    }
    for (int i = 0; i < num_entries_; i++) {
      while (i == row_indices_[row + 1]) {
        row++;
      }
      result->add_number_to_element(column_indices_[i],
                                    values_[i] * p->get_value(row));
    }
  } else {
    int col;
    // find the col corresponding to the first nonzero entry
    for (int i = 1; i < num_columns_ + 1; i++) {
      if (column_indices_[i] > 0) {
        col = i - 1;
        break;
      }
    }
    for (int i = 0; i < num_entries_; i++) {
      // go to the next col
      while (i == column_indices_[col + 1]) {
        col++;
      }
      result->add_number_to_element(col,
                                    values_[i] * p->get_value(row_indices_[i]));
    }
  }
}

void SpHbMat::multiply(std::shared_ptr<const Vector> p,
                       std::shared_ptr<Vector> result) const
{

  result->set_to_zero();

  if (is_compressed_row_format_) {
    int row;
    // find the row corresponding to the first nonzero entry
    for (int i = 1; i < num_rows_ + 1; i++) {
      if (row_indices_[i] > 0) {
        row = i - 1;
        break;
      }
    }
    for (int i = 0; i < num_entries_; i++) {
      // go to the next row
      while (i == row_indices_[row + 1]) {
        row++;
      }
      result->add_number_to_element(row, values_[i] *
                                             p->get_value(column_indices_[i]));
    }
  } else {
    int col;
    // find the col corresponding to the first nonzero entry
    for (int i = 1; i < num_columns_ + 1; i++) {
      if (column_indices_[i] > 0) {
        col = i - 1;
        break;
      }
    }
    for (int i = 0; i < num_entries_; i++) {
      // go to the next col
      while (i == column_indices_[col + 1]) {
        col++;
      }
      result->add_number_to_element(row_indices_[i],
                                    values_[i] * p->get_value(col));
    }
  }
}

//@{
void SpHbMat::print_full(const char* name,
                         Ipopt::SmartPtr<Ipopt::Journalist> jnlst,
                         Ipopt::EJournalLevel level,
                         Ipopt::EJournalCategory category) const
{
  char mat_val[99];
  auto dense_matrix = new double[num_rows_ * num_columns_]();
  int row;
  if (is_compressed_row_format_) {
    for (int i = 1; i < num_rows_ + 1; i++) {
      if (row_indices_[i] > 0) {
        row = i - 1;
        break;
      }
    }
    for (int i = 0; i < num_entries_; i++) {
      while (i == row_indices_[row + 1]) {
        row++;
      }
      dense_matrix[num_columns_ * row + column_indices_[i]] = values_[i];
    }
  } else {
    int col;
    for (int i = 1; i < num_columns_ + 1; i++) {
      if (column_indices_[i] > 0) {
        col = i - 1;
        break;
      }
    }
    for (int i = 0; i < num_entries_; i++) {
      while (i == column_indices_[col + 1]) {
        col++;
      }
      dense_matrix[num_columns_ * row_indices_[i] + col] = values_[i];
    }
  }

  if (!IsNull(jnlst)) {
    if (name != nullptr) {
      jnlst->Printf(level, category, name);
      jnlst->Printf(level, category, " =: {\n");
    }
    for (int i = 0; i < num_rows_; i++) {
      for (int j = 0; j < num_columns_; j++) {
        sprintf(mat_val, "%f  ", dense_matrix[i * get_num_columns() + j]);
        jnlst->Printf(level, category, mat_val);
      }
      jnlst->Printf(level, category, "\n");
    }
    jnlst->Printf(level, category, "}\n\n");
  } else {
    if (name != nullptr)
      printf("%s =:{\n", name);

    for (int i = 0; i < num_rows_; i++) {
      for (int j = 0; j < num_columns_; j++) {
        printf("%23.16e  ", dense_matrix[i * get_num_columns() + j]);
      }
      printf("\n");
    }
    printf("}\n\n");
  }
  delete[] dense_matrix;
}

void SpHbMat::print(const char* name, Ipopt::SmartPtr<Ipopt::Journalist> jnlst,
                    Ipopt::EJournalLevel level,
                    Ipopt::EJournalCategory category) const
{

  if (is_compressed_row_format_) {
    std::cout << name << "= " << std::endl;
    std::cout << "ColIndex: ";
    for (int i = 0; i < num_entries_; i++)
      std::cout << get_column_index_at_entry(i) << " ";
    std::cout << " " << std::endl;

    std::cout << "RowIndex: ";
    for (int i = 0; i < num_rows_ + 1; i++)
      std::cout << get_row_index_at_entry(i) << " ";
    std::cout << " " << std::endl;
  } else {
    // for compressed column format
    std::cout << name << "= " << std::endl;
    std::cout << "ColIndex: ";
    for (int i = 0; i < num_columns_ + 1; i++)
      std::cout << get_column_index_at_entry(i) << " ";

    std::cout << " " << std::endl;
    std::cout << "RowIndex: ";

    for (int i = 0; i < num_entries_; i++)
      std::cout << get_row_index_at_entry(i) << " ";
    std::cout << " " << std::endl;
  }
  std::cout << "MatVal:   ";

  for (int i = 0; i < num_entries_; i++)
    std::cout << get_value_at_entry(i) << " ";
  std::cout << " " << std::endl;

  std::cout << "order:    ";
  for (int i = 0; i < num_entries_; i++)
    std::cout << get_order_at_entry(i) << " ";
  std::cout << " " << std::endl;
}
#if 0
/** @name norms */
//@{
const double SpHbMat::calc_one_norm() const
{
  // TODO: test on it!
  std::shared_ptr<Vector> colSums = std::make_shared<Vector>(num_columns_);
  if (is_compressed_row_format_) {
    for (int i = 0; i < num_entries_; i++) {
      colSums->add_number_to_element(get_column_indices(i), abs(values_[i]));
    }
  } else {
    int col = 0;
    for (int i = 1; i < num_columns_ + 1; i++) {
      if (column_indices_[i] > 0) {
        col = i - 1;
        break;
      }
    }
    for (int i = 0; i < num_entries_; i++) {
      while (column_indices_[col + 1])
        col++;
      colSums->add_number_to_element(col, abs(values_[i]));
    }
  }

  double oneNorm =
      colSums->calc_inf_norm(); // same as calculating the MAX of an array
  return oneNorm;
}

const double SpHbMat::calc_inf_norm() const
{
  // TODO: test on it!
  std::shared_ptr<Vector> rowSums = std::make_shared<Vector>(num_rows_);
  if (is_compressed_row_format_) {
    for (int i = 0; i < num_entries_; i++) {
      rowSums->add_number_to_element(get_row_indices(i), abs(values_[i]));
    }
  } else {
    int row = 0;
    for (int i = 1; i < num_rows_ + 1; i++) {
      if (row_indices_[i] > 0) {
        row = i - 1;
        break;
      }
    }
    for (int i = 0; i < num_entries_; i++) {
      while (row_indices_[row + 1])
        row++;
      rowSums->add_number_to_element(row, abs(values_[i]));
    }
  }

  double InfNorm =
      rowSums->calc_inf_norm(); // same as calculating the MAX of an array
  return InfNorm;
}
#endif
//@}
} // END_OF_NAMESPACE
