/* Copyright (C) 2019
 * All Rights Reserved.
 *
 * Authors: Xinyi Luo
 * Date:2019-10
 */

#ifndef SQPHOTSTART_SPTRIPLETMAT_HPP_
#define SQPHOTSTART_SPTRIPLETMAT_HPP_

#include "restartsqp/Utils.hpp"
#include "restartsqp/Vector.hpp"

namespace RestartSqp {

// Matrix base class.  It is very simple since we onlt need to do matrix-vector
// products with general matrices. */
class Matrix
{
public:
  /** Default constructor. */
  Matrix()
  {
  }

  /** Default destructor. */
  virtual ~Matrix()
  {
  }

  /** This methods add the multiple of the matrix with the vector p to the
   * result vector. */
  virtual void multiply(std::shared_ptr<const Vector> p,
                        std::shared_ptr<Vector> result,
                        double factor = 1.) const = 0;

  /** This methods add the multiple of the transpose of this matrix with the
   * vector p to the result vector. */
  virtual void multiply_transpose(std::shared_ptr<const Vector> p,
                                  std::shared_ptr<Vector> result,
                                  double factor = 1.) const = 0;

private:
  /** Copy Constructor */
  Matrix(const Matrix&);

  /** Overloaded Equals Operator */
  void operator=(const Matrix&);
};

// Forward definition
class SparseHbMatrix;

/** This is the class for SparseMatrix, it stores the sparse matrix
 *  data in triplet, in class member @values_. It contains the methods
 *  that can copy a Matrix, allocate data to the class member, and
 *  perform a matrix vector multiplication.
 */

class SparseTripletMatrix : public Matrix
{
public:
  /** constructor/destructor */
  //@{
  /** Constructor for an empty Sparse Matrix with N non-zero entries.
   *  It sets the property of the matrix, but */
  SparseTripletMatrix(int nnz, int num_rows, int num_columns,
                      bool is_symmetric = false, bool allocate = true);

  /**
   *@brief
   *
   */
  SparseTripletMatrix(const double* data, int num_rows, int num_columns,
                      bool row_oriented);

  /** Constructor from a sparse Harwell Boeing matrix. */
  SparseTripletMatrix(std::shared_ptr<const SparseHbMatrix> sp_hb_matrix);

  /** Destructor*/
  ~SparseTripletMatrix();
  //@}

  /**
   *@brief print the sparse matrix in triplet form
   */
  void print(const std::string& matrix_name,
             Ipopt::SmartPtr<Ipopt::Journalist> jnlst = nullptr,
             Ipopt::EJournalLevel level = Ipopt::J_ALL,
             Ipopt::EJournalCategory category = Ipopt::J_DBG) const;

  /**
   * @brief print the sparse matrix in the sense form
   */
  void print_dense(const char* name,
                   Ipopt::SmartPtr<Ipopt::Journalist> jnlst = nullptr,
                   Ipopt::EJournalLevel level = Ipopt::J_ALL,
                   Ipopt::EJournalCategory category = Ipopt::J_DBG) const;

  /** This methods add the multiple of the matrix with the vector p to the
   * result vector. */
  void multiply(std::shared_ptr<const Vector> p, std::shared_ptr<Vector> result,
                double factor = 1.) const override;

  /** This methods add the multiple of the transpose of this matrix with the
   * vector p to the result vector. */
  void multiply_transpose(std::shared_ptr<const Vector> p,
                          std::shared_ptr<Vector> result,
                          double factor = 1.) const override;

  //@}
  /**
   * @brief get the dense matrix corresponding to the matrix data stored in
   * class
   * members. The stored matrix will be column oriented.
   */
  void get_dense_matrix(double* dense_matrix, bool row_oriented = true) const;

  /**
   * @brief calculate the one norm of the matrix
   *
   * @return the calculated one-norm
   */
  bool is_compressed_row_format() const
  {
    return false;
  };

  bool is_symmetric() const
  {
    return is_symmetric_;
  }

  inline bool is_initialized() const
  {
    return is_initialized_;
  }

  /**
   *@brief make a deep copy of a matrix information
   */
  virtual void copy(std::shared_ptr<const SparseTripletMatrix> rhs,
                    bool deep_copy = true);

  /**@name Extract Matrix info*/
  //@{
  inline int get_num_columns() const
  {
    return num_columns_;
  }

  inline int get_num_rows() const
  {
    return num_rows_;
  }

  inline int get_num_entries() const
  {
    return num_entries_;
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
    return order_;
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
    return order_;
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
    return order_[i];
  }
  //@}

  inline void set_value_at_entry(int location, double value_to_assign)
  {
    values_[location] = value_to_assign;
  }

  inline void set_row_index_at_entry(int i, int value)
  {
    row_indices_[i] = value;
  }

  inline void set_column_index_at_entry(int i, int value)
  {
    column_indices_[i] = value;
  }

  ///////////////////////////////////////////////////////////
  //                     PRIVATE  METHODS                  //
  ///////////////////////////////////////////////////////////
private:
  /** Default constructor*/
  SparseTripletMatrix();

  /** Copy Constructor */
  SparseTripletMatrix(const SparseTripletMatrix&);

  /** Overloaded Equals Operator */
  void operator=(const SparseTripletMatrix&);

  /** free all memory*/
  void free_memory_();

  ///////////////////////////////////////////////////////////
  //                     PRIVATE  MEMBERS                  //
  ///////////////////////////////////////////////////////////

private:
  bool is_initialized_;
  bool is_allocated_;
  bool is_symmetric_;   /**< is the matrix symmetric, if yes, the non-diagonal
                                 *data will only be stored for once*/
  double* values_;      /**< the entry data of a matrix */
  int num_columns_;     /**< the number columns of a matrix */
  int num_entries_;     /**< number of non-zero entries in  matrix */
  int num_rows_;        /**< the number of rows of a matrix */
  int* column_indices_; /**< the column number of a matrix entry */
  int* row_indices_;    /**< the row number of a matrix entry */
  int* order_; /**< the corresponding original position of a matrix entry */
};
}

#endif
