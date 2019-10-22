/* Copyright (C) 2019
* All Rights Reserved.
*
* Authors: Xinyi Luo
* Date:2019-05
*/
#include <sqphot/Utils.hpp>
namespace SQPhotstart {

void print_matrix(double *M, int length, int width) {
    for (int row = 0; row < length; row++) {
        for (int col = 0; col < width; col++) {
            std::cout << M[row * width + col] << " ";
        }
        std::cout << "\n" << std::endl;
    }
}


bool isFinite(double *x, int length) {
    for (int i = 0; i < length; i++) {
        if (x[i] < INF && x[i] > -INF) {
            return true;
        }
    }
    return false;
}

ConstraintType classify_single_constraint(double lower_bound, double upper_bound) {
    if(lower_bound>-INF && upper_bound<INF) {
        if ((upper_bound-lower_bound)<1.0e-8) {
            return EQUAL;
        } else
            return BOUNDED;
    }
    else if(lower_bound>-INF&&upper_bound>INF) {
        return BOUNDED_BELOW;
    }
    else if(upper_bound<INF&&lower_bound<-INF) {
        return BOUNDED_ABOVE;
    }
    else {
        return UNBOUNDED;
    }
}


bool is_int_array_equal(const int* a, const int* b, int length) {
    for (int i = 0; i <length; i++) {
        if(a[i] != b[i])
            return false;
    }
    return true;
}

bool is_double_array_equal(const double* a, const double* b, int length) {
    for (int i = 0; i <length; i++) {
        if(fabs(a[i] - b[i])>1.0e-8)
            return false;
    }
    return true;
}


double oneNorm(const double* x, int n) {
    double sum = 0;
    for (int i = 0; i < n; i++) {
        if (x[i] < 0)sum -= x[i];
        else sum += x[i];
    }
    return sum;
}

double infNorm(const double* x, int n) {
    double infnorm = 0;
    for (int i = 0; i < n; i++) {
        double absxk;
        if (x[i] < 0) absxk = -x[i];
        else absxk = x[i];
        if (absxk > infnorm) infnorm = absxk;
    }
    return infnorm;
}

}
