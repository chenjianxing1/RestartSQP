/* Copyright (C) 2019
* All Rights Reserved.
*
* Authors: Xinyi Luo
* Date:2019-07
*/

#include <iostream>

#ifndef __SQPDEBUG_HPP__
#define __SQPDEBUG_HPP__

#define DEBUG false

#ifdef DEBUG
#define COMPARE_QP_SOLVER false
#define CHECK_TERMINATION false
#define CHECK_TR_ALG true //check trust region algorithm
#define CHECK_SOC false
#define CHECK_INFEA_CAL false
#define CHECK_QP_SOLVER false
#define CHECK_NLP_READER false
#define PRINT_QP_DATA false
#define CHECK_QP_INFEASIBILITY false
#define GET_QP_INTERFACE_MEMBERS false
#define PRINT_OUT_QP_WITH_ERROR false
#if PRINT_OUT_QP_WITH_ERROR
#define PRINT_QP_IN_CPP false
#endif
#endif

#endif /* __SQPDEBUG_HPP */