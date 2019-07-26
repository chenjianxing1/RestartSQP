/* Copyright (C) 2019
 * All Rights Reserved.
 *
 * Authors: Xinyi Luo
 * Date:2019-06
 */
#ifndef SQPHOTSTART_UTILS_HPP_
#define SQPHOTSTART_UTILS_HPP_

#include <iostream>
#include <memory>
#include <cstdio>
#include <vector>
#include <stdlib.h>
#include <assert.h>
#include <iterator>
#include <string>
#include <algorithm>

#include <sqphot/SQPDebug.hpp>
#include <sqphot/Types.hpp>

using namespace std;
namespace SQPhotstart {
/*Print a matrix with given size*/
void print_matrix(Number* M, Index length, Index width);

const Number INF = 1.0e19;

/* check if x is finite*/
bool isFinite(Number* x, Index length);

/* Calculate the one norm of a n-dimension vector x*/
Number oneNorm(const Number* x, Index n);


/* Calculate the infinity norm of a n-dimension vector x*/
Number infNorm(const Number* x, Index n);

ConstraintType classify_single_constraint(Number lower_bound, Number upper_bound);



//debug tool print things out
template<typename T>
inline void print_(char* name, T* vec, Index length) {
    printf(" %s is \n", name);
    for (int i = 0; i < length; i++) {
        std::cout << vec[i] << std::endl;
    }
    printf("end of %s\n", name);
}


} // namespace SQPhotstart
#endif
