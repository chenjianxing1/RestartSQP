/* Copyright (C) 2019
 * All Rights Reserved.
 *
 * Authors: Xinyi Luo
 * Date:2019-07
 */
#ifndef SQPHOTSTART_MATRIX_HPP_
#define SQPHOTSTART_MATRIX_HPP_

#include "IpJournalist.hpp"
#include "sqphot/Vector.hpp"

#include <memory>

namespace SQPhotstart {

// This structure stores the positions of multiples of identity matrices in the
// Jacobian.
// TODO: Make this a claa
typedef struct
{
  int length;
  int* irow;
  int* jcol;
  int* size;
  double* value;
} IdentityMatrixPositionsStruct;

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
 *@brief
 * This is a virtual base class...
 *
 */
class Matrix
{

  ///////////////////////////////////////////////////////////
  //                     PUBLIC  METHODS                   //
  ///////////////////////////////////////////////////////////

public:
  /** Default constructor*/
  Matrix() = default;

  /** Default destructor*/
  virtual ~Matrix() = default;

  /**
   * @brief Print Matrix. If matrix is sparse, then print it in sparse form.
   */
  virtual void print(const char* name = nullptr,
                     Ipopt::SmartPtr<Ipopt::Journalist> jnlst = nullptr,
                     Ipopt::EJournalLevel level = Ipopt::J_ALL,
                     Ipopt::EJournalCategory category = Ipopt::J_DBG) const = 0;

  virtual void
  print_full(const char* name = nullptr,
             Ipopt::SmartPtr<Ipopt::Journalist> jnlst = nullptr,
             Ipopt::EJournalLevel level = Ipopt::J_ALL,
             Ipopt::EJournalCategory category = Ipopt::J_DBG) const = 0;

  virtual bool is_compressed_row_format() const = 0;

  virtual bool is_initialized() const = 0;

  virtual int get_num_entries() const = 0;

  virtual int get_num_rows() const = 0;

  virtual int get_num_columns() const = 0;

  virtual bool is_symmetric() const = 0;

  virtual double get_value_at_entry(int i) const = 0;

  virtual int get_column_index_at_entry(int i) const = 0;

  virtual int get_row_index_at_entry(int i) const = 0;

  virtual int get_order_at_entry(int i) const = 0;

  virtual const double* get_values() const = 0;

  virtual const int* get_column_indices() const = 0;

  virtual const int* get_row_indices() const = 0;

  virtual const int* get_order() const = 0;

  virtual double* get_nonconst_values() = 0;

  virtual int* get_nonconst_column_indices() = 0;

  virtual int* get_nonconst_row_indices() = 0;

  virtual int* get_nonconst_order() = 0;

  virtual void multiply(std::shared_ptr<const Vector> p,
                        std::shared_ptr<Vector> result) const = 0;

  virtual void multiply_transpose(std::shared_ptr<const Vector> p,
                                  std::shared_ptr<Vector> result) const = 0;

  ///////////////////////////////////////////////////////////
  //                     PRIVATE  METHODS                  //
  ///////////////////////////////////////////////////////////

private:
  /** Copy Constructor */
  Matrix(const Matrix&);

  /** Overloaded Equals Operator */
  void operator=(const Matrix&);
};
}

#endif // SQPHOTSTART_MATRIX_HPP_
