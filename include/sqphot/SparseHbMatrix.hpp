
/* Copyright (C) 2019
 * All Rights Reserved.
 *
 * Authors: Xinyi Luo
 * Date:2019-10
 */

#ifndef SQPHOTSTART_SPHBMAT_HPP_
#define SQPHOTSTART_SPHBMAT_HPP_

#include "sqphot/SpTripletMat.hpp"

namespace SQPhotstart {

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
   *  This constructor initializes the sizes and format type, but does not allocate memory.
   */
  SparseHbMatrix(int num_rows, int num_columns, bool is_compressed_row);

  /**
   *  This constructor initializes sizes and format type, and allocates memory.
   */
  SparseHbMatrix(int num_entries_, int num_rows, int num_columns,
                 bool is_compressed_row);

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
  ~SparseHbMatrix() override;
  //@}

  //@{
  void set_structure(std::shared_ptr<const SpTripletMat> rhs,
                     IdentityMatrixPositions& identity_matrix_positions);

  void set_structure(std::shared_ptr<const SpTripletMat> rhs);
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
#if 0
  void set_values(std::shared_ptr<const SpTripletMat> triplet_matrix,
                 IdentityMatrixPositions& identity_matrix_positions);
#endif

  void set_values(std::shared_ptr<const SpTripletMat> triplet_matrix);

  void get_dense_matrix(double* dense_matrix, bool row_oriented = true) const;

  /**
   * @brief print the matrix information
   */
  void print(const char* name = nullptr,
             Ipopt::SmartPtr<Ipopt::Journalist> jnlst = nullptr,
             Ipopt::EJournalLevel level = Ipopt::J_ALL,
             Ipopt::EJournalCategory category = Ipopt::J_DBG) const override;

  void
  print_full(const char* name = nullptr,
             Ipopt::SmartPtr<Ipopt::Journalist> jnlst = nullptr,
             Ipopt::EJournalLevel level = Ipopt::J_ALL,
             Ipopt::EJournalCategory category = Ipopt::J_DBG) const override;

  void multiply(std::shared_ptr<const Vector> p,
                std::shared_ptr<Vector> result) const override;

  void multiply_transpose(std::shared_ptr<const Vector> p,
                          std::shared_ptr<Vector> result) const override;

  // TODO: Maybe replace by something not taking arrays directly?
  void multiply_transpose(const double* p, double* result) const
  {
    for (int i = 0; i < num_columns_; i++)
      result[i] = 0.0;

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
        result[column_indices_[i]] += values_[i] * p[row];
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
        result[col] += values_[i] * p[row_indices_[i]];
      }
    }
  }

  /**
   * @brief make a deep copy of a matrix information
   */

  void copy(std::shared_ptr<const SparseHbMatrix> rhs);
#if 0
  const double calc_one_norm() const;

  const double calc_inf_norm() const;
#endif
  /**
   * @brief convert the matrix data stored in the class members to a triplet
   * matrix
   * speci fied by rhs */
  std::shared_ptr<SpTripletMat> convert_to_triplet() const;

  /** Extract class member information*/
  //@{

  inline int get_num_entries() const override
  {
    return num_entries_;
  }

  inline int get_num_columns() const override
  {
    return num_columns_;
  }

  inline int get_num_rows() const override
  {
    return num_rows_;
  }

  inline int get_row_index_at_entry(int i) const override
  {
    return row_indices_[i];
  }

  inline int get_column_index_at_entry(int i) const override
  {
    return column_indices_[i];
  }

  inline double get_value_at_entry(int i) const override
  {
    return values_[i];
  }

  inline int get_order_at_entry(int i) const override
  {
    return triplet_order_[i];
  }
#if 0
  inline int get_row_indices(int i) const
  {
    return RowIndex_[i];
  }

  inline int get_column_indices(int i) override
  {

    return ColIndex_[i];
  }

  inline double get_value_at_entry(int i) override
  {

    return MatVal_[i];
  }

  inline int get_order(int i) override
  {

    return order_[i];
  }
#endif
  inline const int* get_row_indices() const override
  {
    return row_indices_;
  }

  inline const int* get_column_indices() const override
  {
    return column_indices_;
  }

  inline const double* get_values() const override
  {
    return values_;
  }

  inline const int* get_order() const override
  {
    return triplet_order_;
  }

  inline int* get_nonconst_row_indices() override
  {
    return row_indices_;
  }

  inline int* get_nonconst_column_indices() override
  {
    return column_indices_;
  }

  inline double* get_nonconst_values() override
  {
    return values_;
  }

  inline int* get_nonconst_order() override
  {
    return triplet_order_;
  }

  inline bool is_symmetric() const override
  {
    return is_symmetric_;
  }

  inline bool is_initialized() const override
  {
    return is_initialized_;
  }

  inline bool is_compressed_row_format() const override
  {
    return is_compressed_row_format_;
  }

  /**
   *@brief write data to a file
   * Only works when DEBUG is enabled.
   */
  void write_to_file(const char* name, Ipopt::SmartPtr<Ipopt::Journalist> jnlst,
                     Ipopt::EJournalLevel level,
                     Ipopt::EJournalCategory category, QpSolver solver);

  ///////////////////////////////////////////////////////////
  //                     PRIVATE  METHODS                  //
  ///////////////////////////////////////////////////////////

private:
  /** Compiler generated methods to hide. */
  //@{
  /** Default constructor*/
  SparseHbMatrix();

  /** Copy Constructor */
  SparseHbMatrix(const SparseHbMatrix&);

  /** Overloaded Equals Operator */
  void operator=(const SparseHbMatrix&);
  //@}

  void set_zero();

  void allocate_memory_();

  /** Set the structure from an (unsorted) list of triplet elements. */
  void set_structure_from_list_(
      std::vector<std::tuple<int, int, int, double>> elements_list);

#if 0
    template <typename T>
    static void print_tuple(std::vector<std::tuple<int,int,T>> tuple) {
        for(int i = 0; i<tuple.size(); i++) {
            printf("%d %d ", std::get<0>(tuple[i]),std::get<1>(tuple[i]));
            std::cout<<std::get<2>(tuple[i])<<std::endl;
        }
    }
#endif

  ///////////////////////////////////////////////////////////
  //                     PRIVATE  MEMBERS                  //
  ///////////////////////////////////////////////////////////

  /** Flag indicating whether the nonzero patterns has been set already. */
  bool is_initialized_;

  /** Flag indicating whether this matrix is symmetric.  That that case, only
   * the upper triangular part is stored. */
  bool is_symmetric_;

  /** Flag indicating whether this matrix is in compressed row (true) or
   * compressed colums (false) formet. */
  bool is_compressed_row_format_;

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

  /** Number of elements that come from the triplet matrix that was used to set
   * the structure (if such a matrix was used.) */
  int num_triplet_entries_;

  /** Values of the nonzero entries.
   *
   *  If this matrix is constructed with and IdentityMatrixPositions object, the
   * values corresponding to these identity matrices are set together with the
   * structure and not changed afterwards. */
  double* values_;

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
