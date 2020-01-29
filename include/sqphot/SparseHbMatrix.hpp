
/* Copyright (C) 2019
 * All Rights Reserved.
 *
 * Authors: Xinyi Luo
 * Date:2019-10
 */

#ifndef SQPHOTSTART_SPHBMAT_HPP_
#define SQPHOTSTART_SPHBMAT_HPP_

#include "sqphot/SparseTripletMatrix.hpp"

namespace RestartSqp {

/** Class that stores the positions of identity matrices in the QP Jacobian. */
class IdentityMatrixPositions
{
public:
  /** Default constructor. */
  IdentityMatrixPositions()
  {
  }

  /** Destructor. */
  ~IdentityMatrixPositions()
  {
  }

  /** Add an identity matrix. */
  void add_matrix(int row_offset, int column_offset, int dimension,
                  double multiplicator)
  {
    row_offsets_.push_back(row_offset);
    column_offsets_.push_back(column_offset);
    dimensions_.push_back(dimension);
    multiplicators_.push_back(multiplicator);
  }

  /** Getter methods. */
  //@{
  /** Return number of matrices. */
  int get_num_matrices() const
  {
    return (int)row_offsets_.size();
  }
  /** Return row offset of matrix i. Counting starts at 1. */
  int get_row_offset(int i) const
  {
    assert(i < row_offsets_.size());
    return row_offsets_[i];
  }
  /** Return column offset of matrix i. Counting starts at 1. */
  int get_column_offset(int i) const
  {
    assert(i < column_offsets_.size());
    return column_offsets_[i];
  }
  /** Return size of matrix i. */
  int get_dimension(int i) const
  {
    assert(i < dimensions_.size());
    return dimensions_[i];
  }
  /** Return multiplicator for matrix i */
  double get_multiplicator(int i) const
  {
    return multiplicators_[i];
  }
  //@}

private:
  /** Copy Constructor */
  IdentityMatrixPositions(const IdentityMatrixPositions&);

  /** Overloaded Equals Operator */
  void operator=(const IdentityMatrixPositions&);

  /** Row offsets. */
  std::vector<int> row_offsets_;

  /** Column offsets. */
  std::vector<int> column_offsets_;

  /** Dimesions of the identity matrices. */
  std::vector<int> dimensions_;

  /** Factors by which the identity matrices are multiplied. */
  std::vector<double> multiplicators_;
};

/**
 *@brief This is a derived class of Matrix.
 * It strored matrix in Harwell-Boeing format which is required by qpOASES and
 *QORE.
 * It contains method to transform matrix format from Triplet form to
 * Harwell-Boeing format and then stored to its class members
 */

class SparseHbMatrix : public Matrix
{

public:
  /** constructor/destructor */
  //@{

  /** Constructor.
   *
   *  This constructor initializes the sizes and format type, but does not
   * allocate memory.
   */
  SparseHbMatrix(int num_rows, int num_columns, bool is_compressed_row,
                 bool is_symmetric = false);

  /**
   *  This constructor initializes sizes and format type, and allocates memory.
   */
  SparseHbMatrix(int num_entries, int num_rows, int num_columns,
                 bool is_compressed_row, bool is_symmetric);

  /**
   * @brief Copy structure and data from a dense matrix.  This can only be done
   * if this matrix has not been set to any structure or values before.
   *
   */
  void copy_from_dense_matrix(const double* data, int get_num_rows,
                              int get_num_columns, bool row_oriented = true,
                              bool is_compressed_row = false);

  /**
   * @brief Default destructor
   */
  ~SparseHbMatrix();
  //@}

  //@{
  void set_structure(std::shared_ptr<const SparseTripletMatrix> rhs,
                     IdentityMatrixPositions& identity_matrix_positions);

  void set_structure(std::shared_ptr<const SparseTripletMatrix> rhs);
  //@}

  //@{
  inline void set_value_at_entry(int i, double value)
  {
    assert(i < num_entries_);
    values_[i] = value;
  }

  //
  inline void set_row_index_at_entry(int i, int value)
  {
    if (is_compressed_row_format_) {
      assert(i < num_rows_ + 1);
    } else {
      assert(i < num_entries_);
    }

    row_indices_[i] = value;
  }

  inline void set_column_index_at_entry(int i, int value)
  {
    if (is_compressed_row_format_) {
      assert(i < num_entries_);
    } else {
      assert(i < num_columns_ + 1);
    }

    column_indices_[i] = value;
  }

  //@}
  /**
   * @brief set the Matrix values to the matrix, convert from triplet format to
   * Harwell-Boeing Matrix format.
   * @param rhs entry values(orders are not yet under permutation)
   * @param I_info the 2 identity matrices information
   */
  void set_values(std::shared_ptr<const SparseTripletMatrix> triplet_matrix);

  void get_dense_matrix(double* dense_matrix, bool row_oriented = true) const;

  /**
   * @brief print the matrix information
   */
  void print(const char* name = nullptr,
             Ipopt::SmartPtr<Ipopt::Journalist> jnlst = nullptr,
             Ipopt::EJournalLevel level = Ipopt::J_ALL,
             Ipopt::EJournalCategory category = Ipopt::J_DBG) const;

  void print_dense(const char* name = nullptr,
                   Ipopt::SmartPtr<Ipopt::Journalist> jnlst = nullptr,
                   Ipopt::EJournalLevel level = Ipopt::J_ALL,
                   Ipopt::EJournalCategory category = Ipopt::J_DBG) const;

  /** This methods add the multiple of the matrix with the vector p to the
   * result vector.  The factor says which multiple of the product is added. */
  void multiply(std::shared_ptr<const Vector> p, std::shared_ptr<Vector> result,
                double factor = 1.) const override;

  /** This methods add the multiple of the transpose of this matrix with the
   * vector p to the result vector.   The factor says which multiple of the
   * product is added. */
  void multiply_transpose(std::shared_ptr<const Vector> p,
                          std::shared_ptr<Vector> result,
                          double factor = 1.) const override;

  /**
   * @brief convert the matrix data stored in the class members to a triplet
   * matrix
   * speci fied by rhs */
  std::shared_ptr<SparseTripletMatrix> convert_to_triplet() const;

  /** Extract class member information*/
  //@{

  inline int get_num_entries() const
  {
    return num_entries_;
  }

  inline int get_num_columns() const
  {
    return num_columns_;
  }

  inline int get_num_rows() const
  {
    return num_rows_;
  }

  inline int get_row_index_at_entry(int i) const
  {
    return row_indices_[i];
  }

  inline int get_column_index_at_entry(int i) const
  {
    return column_indices_[i];
  }

  inline double get_value_at_entry(int i) const
  {
    return values_[i];
  }

  inline int get_order_at_entry(int i) const
  {
    return triplet_order_[i];
  }

  inline const int* get_row_indices() const
  {
    return row_indices_;
  }

  inline const int* get_column_indices() const
  {
    return column_indices_;
  }

  inline const double* get_values() const
  {
    return values_;
  }

  inline const int* get_order() const
  {
    return triplet_order_;
  }

  inline int* get_nonconst_row_indices()
  {
    return row_indices_;
  }

  inline int* get_nonconst_column_indices()
  {
    return column_indices_;
  }

  inline double* get_nonconst_values()
  {
    return values_;
  }

  inline int* get_nonconst_order()
  {
    return triplet_order_;
  }

  inline bool is_symmetric() const
  {
    return is_symmetric_;
  }

  inline bool is_initialized() const
  {
    return is_initialized_;
  }

  inline bool is_compressed_row_format() const
  {
    return is_compressed_row_format_;
  }

  /** Write the matrix to a file, using the provided file pointer.  The output
   * is in triplet format. */
  void write_to_file(FILE* file, const std::string& matrix_name) const;

  ///////////////////////////////////////////////////////////
  //                     PRIVATE  METHODS                  //
  ///////////////////////////////////////////////////////////

private:
  /** @name Compiler generated methods to hide. */
  //@{
  /** Default constructor*/
  SparseHbMatrix();

  /** Copy Constructor */
  SparseHbMatrix(const SparseHbMatrix&);

  /** Overloaded Equals Operator */
  void operator=(const SparseHbMatrix&);
  //@}

  void allocate_memory_();

  /** Set the structure from an (unsorted) list of triplet elements. */
  void set_structure_from_list_(
      std::vector<std::tuple<int, int, int, double>>& elements_list);

  /** Add entries of the triplet matrix to the list of elements. */
  void add_triplet_to_element_list_(
      std::shared_ptr<const SparseTripletMatrix> triplet_matrix,
      std::vector<std::tuple<int, int, int, double>>& ele_list);

  ///////////////////////////////////////////////////////////
  //                     PRIVATE  MEMBERS                  //
  ///////////////////////////////////////////////////////////

  /** Flag indicating whether the nonzero patterns has been set already. */
  bool is_initialized_;

  /** Flag indicating whether this matrix is in compressed row (true) or
   * compressed colums (false) formet. */
  bool is_compressed_row_format_;

  /** Flag indicating whether this matrix is symmetric.  That that case, only
   * the upper triangular part is stored. */
  bool is_symmetric_;

  /** Nuber of rows. */
  int num_rows_;

  /** Number of columns. */
  int num_columns_;

  /** Number of nonzero entries. */
  int num_entries_;

  /** Column indices. */
  int* column_indices_;

  /** Row indices. */
  int* row_indices_;

  /** Values of the nonzero entries.
   *
   *  If this matrix is constructed with and IdentityMatrixPositions object, the
   * values corresponding to these identity matrices are set together with the
   * structure and not changed afterwards. */
  double* values_;

  /** Number of elements that come from the triplet matrix that was used to set
   * the structure (if such a matrix was used.) */
  int num_triplet_entries_;

  /** Permutation that tells u that the i-th element in the triple matrix is the
   * order[i]-th element in this matrix. */
  int* triplet_order_;

  /**
   * @brief setup the structure of the sparse matrix for QPsolvrs
   * This method should be only called for once
   *
   * This method will convert the strucutre information from the triplet form
   * from a
   * SpMatrix object to the format required by the corresponding QPsolvers
   *
   * @param rhs a SpMatrix object whose content will be copied to the class
   * members
   * (in a different sparse matrix representations)
   * @param I_info the information of 2 identity sub matrices.
   *
   */
};
}

#endif
