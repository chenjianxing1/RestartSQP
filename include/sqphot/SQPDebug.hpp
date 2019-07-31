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

#define CHECK_LINEAR_ALGEBRA false


#ifdef DEBUG
#define CHECK_TERMINATION false
#define CHECK_TR_ALG false //check trust region algorithm
#define CHECK_SOC false
#define CHECK_INFEA_CAL false
#define CHECK_QP_SOLVER false
#define CHECK_NLP_READER false
#define CHECK_QP_INFEASIBILITY true
#endif

#endif /* __SQPDEBUG_HPP */
