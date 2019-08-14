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
#define MAX(a,b) (a>b)? a :b
#define ABS(a) MAX(a,-a)

/*Print a matrix with given size*/
void print_matrix(double* M, int length, int width);

const double INF = 1.0e19;

/* check if x is finite*/
bool isFinite(double* x, int length);

/* Calculate the one norm of a n-dimension vector x*/
double oneNorm(const double* x, int n);


/* Calculate the infinity norm of a n-dimension vector x*/
double infNorm(const double* x, int n);

ConstraintType classify_single_constraint(double lower_bound, double upper_bound);



//debug tool print things out
template<typename T>
inline void print_(char* name, T* vec, int length) {
    printf(" %s is \n", name);
    for (int i = 0; i < length; i++) {
        std::cout << vec[i] << std::endl;
    }
    printf("end of %s\n", name);
}


} // namespace SQPhotstart
#endif
