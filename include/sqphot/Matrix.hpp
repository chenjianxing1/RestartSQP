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
