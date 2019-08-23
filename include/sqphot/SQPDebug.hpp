/* Copyright (C) 2019
* All Rights Reserved.
*
* Authors: Xinyi Luo
* Date:2019-07
*/

#include <iostream>

#ifndef __SQPDEBUG_HPP__
#define __SQPDEBUG_HPP__

#define DEBUG true

#define CHECK_LINEAR_ALGEBRA false


#ifdef DEBUG
#define CHECK_TERMINATION true
#define  CHECK_COMPLEMENTARITY true
#define CHECK_TR_ALG false //check trust region algorithm
#define CHECK_SOC false
#define CHECK_INFEA_CAL false
#define CHECK_QP_SOLVER false
#define CHECK_NLP_READER false
#define CHECK_QP_INFEASIBILITY false
#define PRINT_OUT_QP_WITH_ERROR false
#endif

#endif /* __SQPDEBUG_HPP */
