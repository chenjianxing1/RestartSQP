/* Copyright (C) 2019
 * All Rights Reserved.
 *
 * Authors: Xinyi Luo
 * Date:2019-10
 */

#ifndef SQPHOTSTART_SPTRIPLETMAT_HPP_
#define SQPHOTSTART_SPTRIPLETMAT_HPP_

#include "sqphot/Matrix.hpp"
#include "sqphot/Utils.hpp"
#include "sqphot/Vector.hpp"

namespace SQPhotstart {

/** This is the class for SparseMatrix, it stores the sparse matrix
 *  data in triplet, in class member @values_. It contains the methods
 *  that can copy a Matrix, allocate data to the class member, and
 *  perform a matrix vector multiplication.
 */

class SpTripletMat : public Matrix
{
public:
  /** constructor/destructor */
  //@{
  /** Constructor for an empty Sparse Matrix with N non-zero entries.
   *  It sets the property of the matrix, but */
  SpTripletMat(int nnz, int num_rows, int num_columns,
               bool is_symmetric = false, bool allocate = true);

  /**
   *@brief
   *
   */
  SpTripletMat(const double* data, int num_rows, int num_columns,
               bool row_oriented);

  /** Default destructor*/
  ~SpTripletMat() override;
  //@}

  /**
   *@brief print the sparse matrix in triplet form
   */
  void print(const char* name,
             Ipopt::SmartPtr<Ipopt::Journalist> jnlst = nullptr,
             Ipopt::EJournalLevel level = Ipopt::J_ALL,
             Ipopt::EJournalCategory category = Ipopt::J_DBG) const override;

  /**
   * @brief print the sparse matrix in the sense form
   */
  void
  print_full(const char* name,
             Ipopt::SmartPtr<Ipopt::Journalist> jnlst = nullptr,
             Ipopt::EJournalLevel level = Ipopt::J_ALL,
             Ipopt::EJournalCategory category = Ipopt::J_DBG) const override;

  /**
   * @brief Multiplies the matrix with a vector p, the pointer to the
   * matrix-vector
   * product will be stored in the class member of another Vector class object
   * called "result"
   * */
  virtual void multiply(std::shared_ptr<const Vector> p,
                        std::shared_ptr<Vector> result) const;

  /**
   * @brief Multiplies the matrix transpose with a vector p, the pointer to the
   * matrix-vector
   * product will be stored in the class member of another Vector class object
   * called "result"
   * */
  virtual void multiply_transpose(std::shared_ptr<const Vector> p,
                                  std::shared_ptr<Vector> result) const;

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
#if 0
  double calc_one_norm() const;
#endif

  bool is_compressed_row_format() const override
  {
    return false;
  };

  bool is_symmetric() const override
  {
    return is_symmetric_;
  }

  inline bool is_initialized() const override
  {
    return is_initialized_;
  }

#if 0
  /**
   * @brief calculate the infinity norm of the matrix
   *
   * @return the calculated inf-norm
   */
  double calc_inf_norm();
#endif
  /**
   *@brief make a deep copy of a matrix information
   */
  virtual void copy(std::shared_ptr<const SpTripletMat> rhs,
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

  inline int get_num_entries() const override
  {
    return num_entries_;
  }

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
    return order_;
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
    return order_;
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
    return order_[i];
  }
#if 0
  inline const int get_row_indices(int i) const
  {
    return row_indices_[i];
  }

  inline const int get_column_indices(int i) const
  {
    return column_indices_[i];
  }

  inline const double get_values(int i) const
  {
    return values_[i];
  }

  inline const int order(int i) const
  {
    return order_[i];
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

  inline const int* order() const
  {
    return order_;
  }
#endif
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

  SpTripletMat();

  /** free all memory*/
  void free_memory_();

  /** Copy Constructor */
  SpTripletMat(const SpTripletMat&);

  /** Overloaded Equals Operator */
  void operator=(const SpTripletMat&);

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
