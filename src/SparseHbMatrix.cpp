#include <algorithm>
#include <iostream>

#include "sqphot/SparseHbMatrix.hpp"

using namespace std;

namespace SQPhotstart {

/** @name constructor/destructor */
//@{
/** Constructor that sets dimensions and type  No memory is allocated. */
SparseHbMatrix::SparseHbMatrix(int num_rows, int num_columns,
                               bool is_compressed_row_format)
 : row_indices_(NULL)
 , column_indices_(NULL)
 , values_(NULL)
 , triplet_order_(NULL)
 , num_entries_(-1)
 , is_initialized_(false)
 , num_rows_(num_rows)
 , num_columns_(num_columns)
 , num_triplet_entries_(-1)
 , is_compressed_row_format_(is_compressed_row_format)
{
}

/**
 *@brief A constructor with the number of non-zero entries, row number and
 *column
 * number specified.
 *
 * @param: num_entries number of nonzero entry
 * @param num_rows: number of rows of a matrix
 * @param num_columns: number of columns of a matrix
 */
SparseHbMatrix::SparseHbMatrix(int num_entries, int num_rows, int num_columns,
                               bool is_compressed_row_format)
 : row_indices_(NULL)
 , column_indices_(NULL)
 , values_(NULL)
 , triplet_order_(NULL)
 , num_entries_(num_entries)
 , num_rows_(num_rows)
 , num_columns_(num_columns)
 , is_initialized_(false)
 , is_compressed_row_format_(is_compressed_row_format)
{
}

void SparseHbMatrix::allocate_memory_()
{
  if (is_compressed_row_format_) {
    column_indices_ = new int[num_entries_];
    row_indices_ = new int[num_rows_ + 1];
  } else {
    column_indices_ = new int[num_columns_ + 1];
    row_indices_ = new int[num_entries_];
  }
  values_ = new double[num_entries_];
  triplet_order_ = new int[num_entries_];
  for (int i = 0; i < num_entries_; i++)
    triplet_order_[i] = i;
}

void SparseHbMatrix::copy_from_dense_matrix(const double* data, int num_rows,
                                            int num_columns, bool row_oriented,
                                            bool is_compressed_row_format)
{
  assert(!is_initialized_);
  assert(row_indices_ == NULL);
  assert(num_rows_ == num_rows);
  assert(num_columns_ == num_columns);
  assert(is_compressed_row_format_ == is_compressed_row_format);

  // First count the number of nonzeros in the dense matrix.
  num_entries_ = 0;
  for (int i = 0; i < num_rows * num_columns; ++i) {
    if (data[i] != 0.) {
      num_entries_++;
    }
  }

  // Allocate memory
  allocate_memory_();

  // Transfer the dense matrix elements into the sparse format, both structure
  // and values
  int i_nnz = 0;
  if (row_oriented) {
    if (is_compressed_row_format) {
      for (int i = 0; i < num_rows; i++) {
        for (int j = 0; j < num_columns; j++) {
          // identify nonzero entry
          if (data[i * num_columns + j] != 0.) {
            values_[i_nnz] = data[i * num_columns + j];
            column_indices_[i_nnz] = j;
            i_nnz++;
          }
        }
        row_indices_[i + 1] = i_nnz;
      }
    } else { // if it is condensed column
      for (int j = 0; j < num_columns; j++) {
        for (int i = 0; i < num_rows; i++) {
          if (data[j + i * num_columns] != 0.) {
            values_[i_nnz] = data[j + i * num_columns];
            row_indices_[i_nnz] = i;
            i_nnz++;
          }
        }
        column_indices_[j + 1] = i_nnz;
      }
    }
  } else {
    if (is_compressed_row_format) {
      for (int j = 0; j < num_rows; j++) {
        for (int i = 0; i < num_columns; i++) {
          if (data[j + i * num_rows] != 0.) {
            values_[i_nnz] = data[j + i * num_rows];
            column_indices_[i_nnz] = i;
            i_nnz++;
          }
        }
        row_indices_[j + 1] = i_nnz;
      }
    } else {
      for (int i = 0; i < num_columns; i++) {
        for (int j = 0; j < num_rows; j++) {
          // identify nonzero entry
          if (data[i * num_rows + j] != 0.) {
            values_[i_nnz] = data[i * num_rows + j];
            row_indices_[i_nnz] = j;
            i_nnz++;
          }
        }
        column_indices_[i + 1] = i_nnz;
      }
    }
  }
  assert(i_nnz == num_entries_);

  is_initialized_ = true;
  triplet_order_ = new int[num_entries_];
}

/**
 *Default destructor
 */
SparseHbMatrix::~SparseHbMatrix()
{
  delete[] column_indices_;
  column_indices_ = NULL;
  delete[] row_indices_;
  row_indices_ = NULL;
  delete[] values_;
  values_ = NULL;
  delete[] triplet_order_;
  triplet_order_ = NULL;
}

//@}

/**
 * @brief This is the sorted rule that used to sort data, first based on
 * column
 * index then based on row index
 */
static bool
tuple_sort_rule_compressed_column(const tuple<int, int, int, double>& left,
                                  const tuple<int, int, int, double>& right)
{

  if (get<1>(left) < get<1>(right))
    return true;
  else if (get<1>(left) > get<1>(right))
    return false;
  else {
    return get<0>(left) < get<0>(right);
  }
}

/**
 * @brief This is the sorted rule that used to sort data, first based on row
 * index then based on column index
 */
static bool
tuple_sort_rule_compressed_row(const tuple<int, int, int, double>& left,
                               const tuple<int, int, int, double>& right)
{

  if (get<0>(left) < get<0>(right))
    return true;
  else if (get<0>(left) > get<0>(right))
    return false;
  else {
    return get<1>(left) < get<1>(right);
  }
}

/** @name setStructure */

//@{

void add_triplet_to_element_list_(
    shared_ptr<const SpTripletMat> triplet_matrix,
    vector<tuple<int, int, int, double>>& ele_list)
{
  bool is_symmetric = triplet_matrix->is_symmetric();

  // Populate the touples with the matrix structure of the triplet matrix
  int num_entries_triplet = triplet_matrix->get_num_entries();
  const int* trip_row_indices = triplet_matrix->get_row_indices();
  const int* trip_column_indices = triplet_matrix->get_column_indices();

  int counter = ele_list.size();
  for (int i = 0; i < num_entries_triplet; i++) {
    ele_list.emplace_back(trip_row_indices[i] - 1, trip_column_indices[i] - 1,
                          counter, 0.);
    counter++;
    if (is_symmetric && (trip_row_indices[i] != trip_column_indices[i])) {
      ele_list.emplace_back(trip_column_indices[i] - 1, trip_row_indices[i] - 1,
                            counter, 0.);
      counter++;
    }
  }
}

void SparseHbMatrix::set_structure_from_list_(
    vector<tuple<int, int, int, double>> elements_list)
{
  // Here we also need to allocated memory for the permutation vector
  triplet_order_ = new int[num_triplet_entries_];

  if (is_compressed_row_format_) {
    // Sort the triplet elements according to the compressed row format
    sort(elements_list.begin(), elements_list.end(),
         tuple_sort_rule_compressed_row);

    int current_row = 0;
    row_indices_[0] = 0;
    for (int i = 0; i < num_entries_; i++) {
      int row_idx = get<0>(elements_list[i]);
      int col_idx = get<1>(elements_list[i]);
      int pos_triplet = get<2>(elements_list[i]);
      column_indices_[i] = col_idx;
      if (pos_triplet >= 0) {
        // This entry comes from a triplet matrix
        triplet_order_[pos_triplet] = i;
      } else {
        // This entry comes from an identity matrix
        double identity_factor = pos_triplet = get<3>(elements_list[i]);
        values_[i] = identity_factor;
      }
      while (current_row != row_idx) {
        row_indices_[++current_row] = i;
      }
    }
    for (int i = current_row + 1; i <= num_rows_; i++) {
      row_indices_[i] = num_entries_;
    }
  } else {
    // Sort the triplet elements according to the compressed column format
    sort(elements_list.begin(), elements_list.end(),
         tuple_sort_rule_compressed_column);

    int current_column = 0;
    column_indices_[0] = 0;
    for (int i = 0; i < num_entries_; i++) {
      int row_idx = get<0>(elements_list[i]);
      int col_idx = get<1>(elements_list[i]);
      int pos_triplet = get<2>(elements_list[i]);
      row_indices_[i] = row_idx;
      assert(row_idx < num_rows_);
      if (pos_triplet >= 0) {
        // This entry comes from a triplet matrix
        triplet_order_[pos_triplet] = i;
      } else {
        // This entry comes from an identity matrix
        double identity_factor = pos_triplet = get<3>(elements_list[i]);
        values_[i] = identity_factor;
      }
      while (current_column != col_idx) {
        current_column++;
        column_indices_[current_column] = i;
      }
    }
    for (int i = current_column + 1; i <= num_columns_; i++) {
      column_indices_[i] = num_entries_;
    }
  }
}

/**
 * @brief setup the structure of the sparse matrix for solver qpOASES(should
 * be called only for once).
 *
 * This method will convert the strucutre information from the triplet form from
 * a
 * SpMatrix object to the format required by the QPsolver qpOASES.
 *
 * @param triplet_matrix a SpMatrix object whose content will be copied to the
 * class
 * members
 * (in a different sparse matrix representations)
 *
 */
void SparseHbMatrix::set_structure(
    shared_ptr<const SpTripletMat> triplet_matrix)
{
  assert(is_initialized_ == false);

  is_symmetric_ = triplet_matrix->is_symmetric();

  // We will need to keep track of the order of the elements in the original
  // matrix
  //
  // In the touples we store (row, column, counter, constant value)
  vector<tuple<int, int, int, double>> elements_list;

  // Put the elements of the triplet matrix into the list
  add_triplet_to_element_list_(triplet_matrix, elements_list);
  num_triplet_entries_ = elements_list.size();

  // Now we have all the information to allocate the memory
  num_entries_ = elements_list.size();
  allocate_memory_();

  set_structure_from_list_(elements_list);

  is_initialized_ = true;
}

/**
 * @brief setup the structure of the sparse matrix for QPsolvrs
 * This method should be only called for once
 *
 * This method will convert the strucutre information from the triplet form from
 *
 * SpMatrix object to the format required by the corresponding QPsolvers
 *
 * @param triplet_matrix a SpMatrix object whose content will be copied to the
 * class
 * members
 * (in a different sparse matrix representations)
 * @param I_info the information of 2 identity sub matrices.
 *
 */
void SparseHbMatrix::set_structure(
    shared_ptr<const SpTripletMat> triplet_matrix,
    IdentityMatrixPositions& identity_matrix_positions)
{
  assert(is_initialized_ == false);
  assert(!triplet_matrix->is_symmetric());
  is_symmetric_ = false;

  // We will need to keep track of the order of the elements in the original
  // matrix
  //
  // In the touples we store (row, column, counter, constant value)
  vector<tuple<int, int, int, double>> elements_list;

  // First put the elements of the triplet matrix into the list
  add_triplet_to_element_list_(triplet_matrix, elements_list);
  num_triplet_entries_ = elements_list.size();

  // Now add the elements from the identity matrices.
  int num_identiy_matrices = identity_matrix_positions.get_num_matrices();
  for (int i = 0; i < num_identiy_matrices; i++) {
    int dim = identity_matrix_positions.get_dimension(i);
    int row_offset = identity_matrix_positions.get_row_offset(i) - 1;
    int col_offset = identity_matrix_positions.get_column_offset(i) - 1;
    double factor = identity_matrix_positions.get_multiplicator(i);
    for (int j = 0; j < dim; j++) {
      elements_list.push_back(
          make_tuple(row_offset + j, col_offset + j, -1, factor));
    }
  }
  num_entries_ = elements_list.size();

  // Now we have all the information to allocate the memory
  allocate_memory_();

  // Now set the structure from the list
  set_structure_from_list_(elements_list);

  is_initialized_ = true;
}

//@}

/** @name setMatVal */
//@{
/**
 * @brief set the Matrix values to the matrix, convert from triplet format to
 * Harwell-Boeing Matrix format.
 * @param triplet_matrix entry values(orders are not yet under permutation)
 * @param I_info struct which stores identity matrices information
 */
#if 0
void SparseHbMatrix::set_values(
    shared_ptr<const SpTripletMat> triplet_matrix,
    IdentityMatrixPositions& identity_matrix_positions)
{
  // adding identity submatrices  to the matrix
  int total_I_entries = 0;
  int num_matrices = identity_matrix_positions.get_num_matrices();
  for (int i = 0; i < num_matrices; i++) {
    total_I_entries += identity_matrix_positions.get_dimension(i);
  }

  // assign each matrix entry to the corresponding position after permutation
  for (int i = 0; i < num_entries_ - total_I_entries; i++) {
    values_[triplet_order_[i]] = triplet_matrix->get_value_at_entry(i);
  }
}
#endif

void SparseHbMatrix::set_values(shared_ptr<const SpTripletMat> triplet_matrix)
{
  int num_trip_entries = triplet_matrix->get_num_entries();
  const int* trip_row_indices = triplet_matrix->get_row_indices();
  const int* trip_col_indices = triplet_matrix->get_column_indices();
  const double* trip_values = triplet_matrix->get_values();

  if (is_symmetric_) {
    int j = 0;
    for (int i = 0; i < num_trip_entries; i++) {
      values_[triplet_order_[j]] = trip_values[i];
      j++;
      if (is_symmetric_ && (trip_col_indices[i] != trip_row_indices[i])) {
        values_[triplet_order_[j]] = trip_values[i];
        j++;
      }
    }
    assert(j == num_triplet_entries_);
  } else {
    assert(num_trip_entries == num_triplet_entries_);

    for (int i = 0; i < num_trip_entries; i++) {
      values_[triplet_order_[i]] = trip_values[i];
    }
  }
}

//@}

void SparseHbMatrix::copy(shared_ptr<const SparseHbMatrix> rhs)
{

  assert(num_entries_ == rhs->get_num_entries());
  assert(num_rows_ == rhs->get_num_rows());
  assert(num_columns_ == rhs->get_num_columns());
  for (int i = 0; i < num_entries_; i++) {
    row_indices_[i] = rhs->get_row_index_at_entry(i);
    values_[i] = rhs->get_value_at_entry(i);
    triplet_order_[i] = rhs->get_order_at_entry(i);
  }

  for (int i = 0; i < num_columns_ + 1; i++) {
    column_indices_[i] = rhs->get_column_index_at_entry(i);
  }
}

void SparseHbMatrix::set_zero()
{
  if (is_compressed_row_format_) {
    for (int i = 0; i < num_entries_; i++) {
      values_[i] = 0;
      column_indices_[i] = 0;
      triplet_order_[i] = i;
    }
    for (int i = 0; i < num_rows_ + 1; i++) {
      row_indices_[i] = 0;
    }
  } else {
    for (int i = 0; i < num_entries_; i++) {
      values_[i] = 0;
      row_indices_[i] = 0;
      triplet_order_[i] = i;
    }
    for (int i = 0; i < num_columns_ + 1; i++) {
      column_indices_[i] = 0;
    }
  }
}

shared_ptr<SpTripletMat> SparseHbMatrix::convert_to_triplet() const
{
  shared_ptr<SpTripletMat> result;

  if (is_symmetric_) {
    // First count the number of nonzeros in the lower triangular part
    int num_trip_entries = 0;
    if (is_compressed_row_format_) {
      for (int row_idx = 0; row_idx < num_rows_; ++row_idx) {
        int j_start = row_indices_[row_idx];
        int j_end = row_indices_[row_idx + 1];
        for (int j = j_start; j < j_end; ++j) {
          int col_idx = column_indices_[j];
          if (row_idx <= col_idx) {
            // this is an element in the lower triangular part
            num_trip_entries++;
          }
        }
      }
    } else {
      // This is now compressed column format
      for (int col_idx = 0; col_idx < num_columns_; ++col_idx) {
        int j_start = column_indices_[col_idx];
        int j_end = column_indices_[col_idx + 1];
        for (int j = j_start; j < j_end; ++j) {
          int row_idx = row_indices_[j];
          if (row_idx <= col_idx) {
            // this is an element in the lower triangular part
            num_trip_entries++;
          }
        }
      }
    }

    // Create the Triplet matrix with allocated memory
    bool allocate = true;
    result = make_shared<SpTripletMat>(num_trip_entries, num_rows_,
                                       num_columns_, is_symmetric_, allocate);

    // Now loop over all elements again and add the lower triangular elements
    int i_ele = 0;
    int* trip_row_indices = result->get_nonconst_row_indices();
    int* trip_col_indices = result->get_nonconst_column_indices();
    double* trip_values = result->get_nonconst_values();
    if (is_compressed_row_format_) {
      for (int row_idx = 0; row_idx < num_rows_; ++row_idx) {
        int j_start = row_indices_[row_idx];
        int j_end = row_indices_[row_idx + 1];
        for (int j = j_start; j < j_end; ++j) {
          int col_idx = column_indices_[j];
          if (row_idx <= col_idx) {
            // this is an element in the lower triangular part
            trip_row_indices[i_ele] = row_idx + 1;
            trip_col_indices[i_ele] = col_idx + 1;
            trip_values[i_ele] = values_[j];
            i_ele++;
          }
        }
      }
    } else {
      // This is now compressed column format
      for (int col_idx = 0; col_idx < num_columns_; ++col_idx) {
        int j_start = column_indices_[col_idx];
        int j_end = column_indices_[col_idx + 1];
        for (int j = j_start; j < j_end; ++j) {
          int row_idx = row_indices_[j];
          if (row_idx <= col_idx) {
            // this is an element in the lower triangular part
            trip_row_indices[i_ele] = row_idx + 1;
            trip_col_indices[i_ele] = col_idx + 1;
            trip_values[i_ele] = values_[j];
            i_ele++;
          }
        }
      }
    }
    assert(i_ele == num_trip_entries);
  } else {
    // This is a non-symmetric matrix
    bool allocate = true;
    result = make_shared<SpTripletMat>(num_entries_, num_rows_, num_columns_,
                                       is_symmetric_, allocate);

    // Now loop over all elements again and add all elements
    int i_ele = 0;
    int* trip_row_indices = result->get_nonconst_row_indices();
    int* trip_col_indices = result->get_nonconst_column_indices();
    double* trip_values = result->get_nonconst_values();
    if (is_compressed_row_format_) {
      for (int row_idx = 0; row_idx < num_rows_; ++row_idx) {
        int j_start = row_indices_[row_idx];
        int j_end = row_indices_[row_idx + 1];
        for (int j = j_start; j < j_end; ++j) {
          int col_idx = column_indices_[j];
          trip_row_indices[i_ele] = row_idx + 1;
          trip_col_indices[i_ele] = col_idx + 1;
          trip_values[i_ele] = values_[j];
          i_ele++;
        }
      }
    } else {
      // This is now compressed column format
      for (int col_idx = 0; col_idx < num_columns_; ++col_idx) {
        int j_start = column_indices_[col_idx];
        int j_end = column_indices_[col_idx + 1];
        for (int j = j_start; j < j_end; ++j) {
          int row_idx = row_indices_[j];
          trip_row_indices[i_ele] = row_idx + 1;
          trip_col_indices[i_ele] = col_idx + 1;
          trip_values[i_ele] = values_[j];
          i_ele++;
        }
      }
    }
    assert(i_ele == num_entries_);
  }

  return result;
}

void SparseHbMatrix::write_to_file(const char* name,
                                   Ipopt::SmartPtr<Ipopt::Journalist> jnlst,
                                   Ipopt::EJournalLevel level,
                                   Ipopt::EJournalCategory category,
                                   QpSolver solver)
{
#ifdef DEBUG
#ifdef PRINT_QP_IN_CPP
  const char* var_type_int;
  const char* var_type_double;
  var_type_int = (solver == QPOASES) ? "sparse_int_t" : "qp_int";
  var_type_double = (solver == QPOASES) ? "real_t" : "double";
  jnlst->Printf(level, category, "%s %s_jc[] = \n{", var_type_int, name);
  int i;
  for (i = 0; i < num_columns_ + 1; i++) {
    if (i % 10 == 0 && i > 1)
      jnlst->Printf(level, category, "\n");
    if (i == num_columns_)
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
    for (i = 0; i < num_rows_ + 1; i++) {
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
    for (i = 0; i < num_columns_ + 1; i++) {
      jnlst->Printf(level, category, "%d\n", ColIndex_[i]);
    }
    for (i = 0; i < EntryNum_; i++) {
      jnlst->Printf(level, category, "%23.16e\n", MatVal_[i]);
    }
  }
#endif
#endif
}

void SparseHbMatrix::get_dense_matrix(double* dense_matrix,
                                      bool row_oriented) const
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

void SparseHbMatrix::multiply_transpose(shared_ptr<const Vector> p,
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

void SparseHbMatrix::multiply(shared_ptr<const Vector> p,
                              shared_ptr<Vector> result) const
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
void SparseHbMatrix::print_full(const char* name,
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

void SparseHbMatrix::print(const char* name,
                           Ipopt::SmartPtr<Ipopt::Journalist> jnlst,
                           Ipopt::EJournalLevel level,
                           Ipopt::EJournalCategory category) const
{

  if (is_compressed_row_format_) {
    cout << name << "= " << endl;
    cout << "ColIndex: ";
    for (int i = 0; i < num_entries_; i++)
      cout << get_column_index_at_entry(i) << " ";
    cout << " " << endl;

    cout << "RowIndex: ";
    for (int i = 0; i < num_rows_ + 1; i++)
      cout << get_row_index_at_entry(i) << " ";
    cout << " " << endl;
  } else {
    // for compressed column format
    cout << name << "= " << endl;
    cout << "ColIndex: ";
    for (int i = 0; i < num_columns_ + 1; i++)
      cout << get_column_index_at_entry(i) << " ";

    cout << " " << endl;
    cout << "RowIndex: ";

    for (int i = 0; i < num_entries_; i++)
      cout << get_row_index_at_entry(i) << " ";
    cout << " " << endl;
  }
  cout << "MatVal:   ";

  for (int i = 0; i < num_entries_; i++)
    cout << get_value_at_entry(i) << " ";
  cout << " " << endl;

  cout << "order:    ";
  for (int i = 0; i < num_entries_; i++)
    cout << get_order_at_entry(i) << " ";
  cout << " " << endl;
}
#if 0
/** @name norms */
//@{
const double SpHbMat::calc_one_norm() const
{
  // TODO: test on it!
  shared_ptr<Vector> colSums = make_shared<Vector>(num_columns_);
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
  shared_ptr<Vector> rowSums = make_shared<Vector>(num_rows_);
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
