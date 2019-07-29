/* Copyright (C) 2019
* All Rights Reserved.
*
* Authors: Xinyi Luo
* Date:2019-05
*/
#include <sqphot/Utils.hpp>
namespace SQPhotstart {

void print_matrix(Number *M, Index length, Index width) {
    for (int row = 0; row < length; row++) {
        for (int col = 0; col < width; col++) {
            std::cout << M[row * width + col] << " ";
        }
        std::cout << "\n" << std::endl;
    }
}


bool isFinite(Number *x, Index length) {
    for (int i = 0; i < length; i++) {
        if (x[i] < INF && x[i] > -INF) {
            return true;
        }
    }
    return false;
}

ConstraintType classify_single_constraint(Number lower_bound, Number upper_bound) {
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


Number oneNorm(const Number* x, Index n) {
    Number sum = 0;
    for (int i = 0; i < n; i++) {
        if (x[i] < 0)sum -= x[i];
        else sum += x[i];
    }
    return sum;
}

Number infNorm(const Number* x, Index n) {
    Number infnorm = 0;
    for (int i = 0; i < n; i++) {
        Number absxk;
        if (x[i] < 0) absxk = -x[i];
        else absxk = x[i];
        if (absxk > infnorm) infnorm = absxk;
    }
    return infnorm;
}

}
