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

namespace SQPhotstart {

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

void Vector::write_to_file(string name, SmartPtr<Journalist> jnlst,
                           EJournalLevel level, EJournalCategory category,
                           QpSolver qpsolver) const
{
// AW Instead of having DEBUG here, set the level accordingly
#ifdef DEBUG
#ifdef PRINT_QP_IN_CPP
  const char* var_type;
  var_type = (qpsolver == QPOASES) ? "real_t" : "double const";

  jnlst->Printf(level, category, "%s %s[%i] = {", var_type, name.c_str(),
                size_);
  for (int i = 0; i < size_; i++) {
    if (i % 10 == 0 && i > 1)
      jnlst->Printf(level, category, "\n");
    if (i == Dim() - 1)
      jnlst->Printf(level, category, "%23.16e};\n\n", values_[i]);
    else
      jnlst->Printf(level, category, "%23.16e, ", values_[i]);
  }
#else
  // print in file
  for (int i = 0; i < size_; i++) {
    jnlst->Printf(level, category, "%23.16e\n", values_[i]);
  }
#endif
#endif
}
}
