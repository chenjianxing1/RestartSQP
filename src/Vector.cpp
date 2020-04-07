/* Copyright (C) 2019
 * All Rights Reserved.
 *
 * Authors: Xinyi Luo
 * Date:2019-07
 */

#include "restartsqp/Vector.hpp"
#include <string>

using namespace Ipopt;
using namespace std;

namespace RestartSqp {

Vector::Vector(int size)
 : size_(size)
 , values_(NULL)
{
  // Reserve the memory
  values_ = new double[size_];
}

Vector::Vector(int size, const double* values)
 : size_(size)
 , values_(NULL)
{
  // Reserve the memory
  // AW: remove initialization
  values_ = new double[size_];

  // Initialize the values from provded values
  for (int i = 0; i < size_; i++) {
    values_[i] = values[i];
  }
}

Vector::Vector(const Vector& rhs)
 : size_(rhs.size_)
{
  values_ = new double[size_];
  for (int i = 0; i < size_; i++) {
    values_[i] = rhs.values_[i];
  }
}

/** Default destructor*/
Vector::~Vector()
{

  // Free the memory
  delete[] values_;
  values_ = NULL;
}

/** print the vector*/
void Vector::print(const string vector_name, SmartPtr<Journalist> jnlst,
                   EJournalLevel level, EJournalCategory category) const
{
  if (!IsValid(jnlst)) {
    write_to_file(stdout, vector_name);
  } else {
    jnlst->Printf(level, category, "Vector %s with %d elements:\n",
                  vector_name.c_str(), size_);
    for (int i = 0; i < size_; i++) {
      jnlst->Printf(level, category, "%5d %23.16e\n", i, values_[i]);
    }
  }
}

void Vector::write_to_file(FILE* file, const string& vector_name) const
{
  fprintf(file, "Vector %s with %d elements:\n", vector_name.c_str(), size_);
  for (int i = 0; i < size_; i++) {
    fprintf(file, "%5d %23.16e\n", i, values_[i]);
  }
}
}
