/* Copyright (C) 2019
* All Rights Reserved.
*
* Authors: Xinyi Luo
* Date:2019-05
*/

#include <iostream>
#include <cmath>
#include "sqphot/Utils.hpp"

namespace SQPhotstart {

void print_matrix(double *M, int length, int width) {
    for (int row = 0; row < length; row++) {
        for (int col = 0; col < width; col++) {
            std::cout << M[row * width + col] << " ";
        }
        std::cout << "\n" << std::endl;
    }
}



bool is_int_array_equal(const int* a, const int* b, int length)
{
  for (int i = 0; i <length; i++) {
    if(a[i] != b[i])
      return false;
  }
  return true;
}

bool is_double_array_equal(const double* a, const double* b, int length)
{
  for (int i = 0; i <length; i++) {
    if(fabs(a[i] - b[i])>1.0e-8)
      return false;
  }
  return true;
}


}
