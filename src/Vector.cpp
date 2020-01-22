/* Copyright (C) 2019
 * All Rights Reserved.
 *
 * Authors: Xinyi Luo
 * Date:2019-07
 */

#include "sqphot/Vector.hpp"
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
void Vector::print(const string name, SmartPtr<Journalist> jnlst,
                   EJournalLevel level, EJournalCategory category) const
{
  if (!IsValid(jnlst)) {
    printf("%s = :\n", name.c_str());

    printf("{ ");
    for (int i = 0; i < size_; i++) {
      printf("%23.16e ", values_[i]);
    }
    printf("}\n\n");
  } else {
    jnlst->Printf(level, category, name.c_str());
    jnlst->Printf(level, category, " =: \n");
    jnlst->Printf(level, category, "{ ");
    for (int i = 0; i < size_; i++) {
      jnlst->Printf(level, category, "%23.16e ", values_[i]);
    }
    jnlst->Printf(level, category, "}\n\n");
  }
}

void Vector::write_to_file(FILE* file, const string& vector_name) const
{
  fprintf(file, "Vector %s with %d elements:\n", vector_name.c_str(), size_);
  for (int i = 0; i < size_; i++) {
    fprintf(file, "%23.16e\n", values_[i]);
  }
}

}
