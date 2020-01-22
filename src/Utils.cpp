/* Copyright (C) 2019
* All Rights Reserved.
*
* Authors: Xinyi Luo
* Date:2019-05
*/

#include "sqphot/Utils.hpp"
#include <chrono>
#include <cmath>
#include <iostream>

#if defined(__MACH__) || defined(__FreeBSD__)
#include <sys/time.h>
#endif
#if !defined(__MSVCRT__)
#include <sys/resource.h>
#endif

namespace RestartSqp {

void print_matrix(double* M, int length, int width)
{
  for (int row = 0; row < length; row++) {
    for (int col = 0; col < width; col++) {
      std::cout << M[row * width + col] << " ";
    }
    std::cout << "\n" << std::endl;
  }
}

bool is_int_array_equal(const int* a, const int* b, int length)
{
  for (int i = 0; i < length; i++) {
    if (a[i] != b[i])
      return false;
  }
  return true;
}

bool is_double_array_equal(const double* a, const double* b, int length)
{
  for (int i = 0; i < length; i++) {
    if (fabs(a[i] - b[i]) > 1.0e-8)
      return false;
  }
  return true;
}

double get_cpu_time_since_start()
{
  // This code is taken from Ipopt.  For windows, we will need to do something
  // different.
  double cpu_time;

  struct rusage usage;
  getrusage(RUSAGE_SELF, &usage);
  cpu_time = (double)usage.ru_utime.tv_sec;
  cpu_time += 1.0e-6 * ((double)usage.ru_utime.tv_usec);

  return cpu_time;
}

// The following might be an issue if there are several threads.  But it is
// probably OK since this static variable is overwritten only once at the very
// beginning, and if there are several changes at the same time, it does not
// matter which survives, since all time measurements should be relative anyway.
static const std::chrono::high_resolution_clock::time_point time_at_start_ =
    std::chrono::high_resolution_clock::now();

double get_wallclock_time_since_start()
{
  std::chrono::high_resolution_clock::time_point now =
      std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> time_span =
      std::chrono::duration_cast<std::chrono::duration<double>>(now -
                                                                time_at_start_);

  return time_span.count();
}
}
