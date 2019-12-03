/* Copyright (C) 2019
 * All Rights Reserved.
 *
 * Authors: Xinyi Luo
 * Date:2019-06
 */

#ifndef SQPHOTSTART_VECTOR_HPP
#define SQPHOTSTART_VECTOR_HPP

#include <cassert>
#include <cmath>
#include <memory>

#include "IpJournalist.hpp"
#include "sqphot/Types.hpp"

namespace SQPhotstart {

/**
 *  This class implements numerical vectors, with some basic
 *  functionality.
 */

class Vector
{
public:
  /** Constructor.
   *
   *  Initialize the size of the vector and allocate the memory.
   */
  Vector(int size);

  /** Constructor that initializes the size of the vector, allocates
   *  memory, and copies the values in vector_values to its internal
   *  array.
   */
  Vector(int size, const double* values);

  /** Copy Constructor */
  Vector(const Vector&);

  /** Destructor*/
  virtual ~Vector();

  /** Print the vector*/
  void print(const std::string name,
             Ipopt::SmartPtr<Ipopt::Journalist> jnlst = nullptr,
             Ipopt::EJournalLevel level = Ipopt::J_ALL,
             Ipopt::EJournalCategory category = Ipopt::J_DBG) const;

  /** add the element in a specific location by a number*/
  void add_number_to_element(int index, double increase_amount)
  {
    if (increase_amount != 0.) {
      values_[index] += increase_amount;
    }
  }

  /** Add a vector to this vector.
   *
   *  This vector += factor * rhs */
  void add_vector(double factor, std::shared_ptr<const Vector> rhs)
  {
    add_vector(factor, *rhs);
  }

  /** Add a vector to this vector.
   *
   *  This vector += factor * rhs */
  void add_vector(double factor, const Vector& rhs)
  {
    assert(size_ == rhs.size_);
    if (factor == 1.) {
      for (int i = 0; i < size_; i++) {
        values_[i] += rhs.values_[i];
      }
    } else if (factor == -1.) {
      for (int i = 0; i < size_; i++) {
        values_[i] -= rhs.values_[i];
      }
    } else {
      for (int i = 0; i < size_; i++) {
        values_[i] += factor * rhs.values_[i];
      }
    }
  }

  /** The values to the sum of two vectors.
   *
   *  this vector = fact1*vec1 + fact2*vec2 */
  void set_to_sum_of_vectors(double fact1, std::shared_ptr<const Vector> vec1,
                             double fact2, std::shared_ptr<const Vector> vec2)
  {
    set_to_sum_of_vectors(fact1, *vec1, fact2, *vec2);
  }

  /** The values to the sum of two vectors.
   *
   *  this vector = fact1*vec1 + fact2*vec2 */
  void set_to_sum_of_vectors(double fact1, const Vector& vec1, double fact2,
                             const Vector& vec2)
  {
    assert(size_ == vec1.size_);
    assert(size_ == vec2.size_);
    if (fact1 == 1. && fact2 == 1.) {
      for (int i = 0; i < size_; i++) {
        values_[i] = vec1.values_[i] + vec2.values_[i];
      }
    } else {
      for (int i = 0; i < size_; i++) {
        values_[i] = fact1 * vec1.values_[i] + fact2 * vec2.values_[i];
      }
    }
  }

  /** Copy all the entries from another vector */
  void copy_vector(std::shared_ptr<const Vector> rhs)
  {
    copy_vector(*rhs);
  }

  /** Copy all the entries from another vector */
  void copy_vector(const Vector& rhs)
  {
    assert(size_ == rhs.size_);
    for (int i = 0; i < size_; i++) {
      values_[i] = rhs.values_[i];
    }
  }

  /** Copy the values from a double* array */
  void copy_values(const double* values)
  {
    for (int i = 0; i < size_; i++) {
      values_[i] = values[i];
    }
  }

  /** set all entries to be 0*/
  void set_to_zero()
  {
    for (int i = 0; i < size_; i++) {
      values_[i] = 0;
    }
  }

  /** Calculate the one-norm of this vector. */
  double one_norm() const
  {
    return one_norm(0, size_);
  }

  /** calculate one norm of a subvector of this vector.  The subvector
      starts at element first_element and has a length of length */
  double one_norm(int first_element, int length) const
  {
    assert(first_element >= 0);
    assert(first_element + length <= size_);
    double oneNorm = 0.;
    for (int i = first_element; i < first_element + length; i++) {
      oneNorm += fabs(values_[i]);
    }
    return oneNorm;
  }

  /** calculate the infinity norm of the member _vector*/
  double inf_norm() const
  {
    double inf_norm = 0.;
    for (int i = 0; i < size_; i++) {
      inf_norm = std::max(inf_norm, fabs(values_[i]));
    }
    return inf_norm;
  }

  /** Compute the inner product of this vector with another */
  double inner_product(std::shared_ptr<Vector> rhs) const
  {
    double product = 0;
    for (int i = 0; i < size_; i++) {
      product += values_[i] * rhs->values()[i];
    }
    return product;
  }

  /** Scale vector by given scalar */
  void scale(double scaling_factor)
  {
    if (scaling_factor != 1.) {
      for (int i = 0; i < size_; i++) {
        values_[i] *= scaling_factor;
      }
    }
  }

  /** Return the dimension of the vector. */
  int dim() const
  {
    return size_;
  }

  /** Return the array with the vector elements (const version) */
  double* values() const
  {
    return values_;
  }

  /** Return the array with the vector elements (non-const version) */
  double* non_const_values()
  {
    return values_;
  }

  /** Return the value of the i-th vector element */
  double value(int i) const
  {
    return values_[i];
  }

  /** Set the value of the i-th element */
  void set_value(int location, double value)
  {
    values_[location] = value;
  }

  /** Print vector or write to a file. */
  void write_to_file(std::string name, Ipopt::SmartPtr<Ipopt::Journalist> jnlst,
                     Ipopt::EJournalLevel level,
                     Ipopt::EJournalCategory category, Solver qpsolver) const;
  // AW: Replace qpsolver argument
  // Also, why are the both a print and a write_to_file method

private:
  /** Default Constructor.  Make this private so that it is not
   *  automatically generated. */
  Vector();

  /** Overloaded Equals Operator.  Make this private so that it is
   *  not automatically generated. */
  void operator=(const Vector&);

  /** Dimension of the vector */
  int size_;

  /** Array with the elements of the vector. */
  double* values_;
};
}
#endif
