/* Copyright (C) 2019
 * All Rights Reserved.
 *
 * Authors: Xinyi Luo
 * Date:2019-06
 */
#ifndef SQPHOTSTART_UTILS_HPP_
#define SQPHOTSTART_UTILS_HPP_

#if 0
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <iterator>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

#include "IpTNLP.hpp"
#include "sqphot/MessageHandling.hpp"
#include "sqphot/SQPDebug.hpp"
#include "sqphot/Types.hpp"
#endif

namespace SQPhotstart {

// AW the following should be options.
const double INF = 1.0e18;
const double m_eps = 1.0e-16;
const double sqrt_m_eps = 1.0e-8;

/** Function that returns true if all elements in the integer arrays a
 *  and b of length length are equal. */
bool is_int_array_equal(const int* a, const int* b, int length);

/** Function that returns true if all elements in the double arrays a
 *  and b of length length are equal. */
bool is_double_array_equal(const double* a, const double* b, int length);

/** Get CPU clock time since start. */
double get_cpu_time_since_start();

/** Get wallclock clock time since start. */
double get_wallclock_time_since_start();

} // namespace SQPhotstart
#endif
