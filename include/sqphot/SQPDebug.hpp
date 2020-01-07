/* Copyright (C) 2019
* All Rights Reserved.
*
* Authors: Xinyi Luo
* Date:2019-07
*/

#include <iostream>

#ifndef __SQPDEBUG_HPP__
#define __SQPDEBUG_HPP__

#define DEBUG
#define NEW_FORMULATION

#ifdef DEBUG

// TODO: The following would be great to have
#undef COMPARE_QP_SOLVER
#undef CHECK_TERMINATION

#define CHECK_TR_ALG // check trust region algorithm
#undef CHECK_SOC
#undef CHECK_INFEA_CAL
#undef CHECK_QP_SOLVER
#undef CHECK_NLP_READER
#undef PRINT_QP_DATA
#define GET_QP_INTERFACE_MEMBERS

#define PRINT_OUT_QP_WITH_ERROR
#ifdef PRINT_OUT_QP_WITH_ERROR
#undef PRINT_DATA_FOR_QPOASES
#define PRINT_DATA_FOR_QORE
#endif
#endif

#endif /* __SQPDEBUG_HPP */
