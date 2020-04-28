# File:   FindUMFPACK.cmake
# Author: Jan Albersmeyer
# Author: Christian Hoffmann
# Date:   2007--2009
#
# This file is part of proprietary software of the
#   Simulation and Optimization Workgroup
#   Interdisciplinary Center for Scientific Computing (IWR)
#   University of Heidelberg, Germany.
# Copyright (C) 2007--2009 by the authors. All rights reserved.
#
###################################################################################################
#
# Find the UMFPACK includes and libraries
#
# This file fulfils the CMAKE modules guidelines:
#   http://www.cmake.org/cgi-bin/viewcvs.cgi/Modules/readme.txt?root=CMake&view=markup
# PLEASE READ THERE BEFORE CHANGING SOMETHING!
#
# Defines:
#   Variable: UMFPACK_FOUND        - TRUE, if the package has been completely found
#   Variable: UMFPACK_INCLUDE_DIRS - List of full paths to include directories required for using UMFPACK
#   Variable: UMFPACK_LIBRARIES    - List of full paths to libraries required for using UMFPACK
#   Function: USE_UMFPACK()      - Convenience function for preparing CMake for usage of UMFPACK.
#
###################################################################################################


####################################################################################################
#### SEARCH REQUIRED PACKAGE RESOURCES
####################################################################################################

SET( SUITESPARSE_DIR "${QORE_SRC_DIR}/suitesparse" CACHE PATH "Path to a SuiteSparse installation.")

MESSAGE( STATUS "Looking for UMFPACK: " )

SET( UMFPACK_DIR ${SUITESPARSE_DIR}/UMFPACK )

FIND_PACKAGE( AMD )

FIND_LIBRARY(CHOLMOD_LIB 
        NAMES  cholmod 
        NO_DEFAULT_PATH
	PATHS
        ${SUITESPARSE_DIR}/lib
		${UMFPACK_DIR}
		${DEFAULT_PACKAGE_DIRS}
		${DEFAULT_LIBRARY_DIRS}
	PATH_SUFFIXES UMFPACK/Lib UMFPACK/lib UMFPACK/LIB Lib lib  lib64 LIB
    )
FIND_PACKAGE_HANDLE_STANDARD_ARGS(cholmod_lib DEFAULT_MSG CHOLMOD_LIB)

MESSAGE( STATUS "Looking for UMFPACK: include directory" )
FIND_PATH( UMFPACK_INCLUDE_DIR umfpack.h
        NO_DEFAULT_PATH
	PATHS
        ${SUITESPARSE_DIR}/include 
		${UMFPACK_DIR}
		${DEFAULT_PACKAGE_DIRS}
		${DEFAULT_INCLUDE_DIRS}
	PATH_SUFFIXES UMFPACK/Include UMFPACK/INCLUDE UMFPACK/include Include INCLUDE include
)
MESSAGE( STATUS "Looking for UMFPACK: include directory (${UMFPACK_INCLUDE_DIR})" )

MESSAGE( STATUS "Looking for UMFPACK: library" )
FIND_LIBRARY( UMFPACK_LIBRARY
	NAMES umfpack
        NO_DEFAULT_PATH
	PATHS
        ${SUITESPARSE_DIR}/lib
		${UMFPACK_DIR}
		${DEFAULT_PACKAGE_DIRS}
		${DEFAULT_LIBRARY_DIRS}
	PATH_SUFFIXES UMFPACK/Lib UMFPACK/lib UMFPACK/LIB Lib lib  lib64 LIB
)
MESSAGE( STATUS "Looking for UMFPACK: library (${UMFPACK_LIBRARY})" )


####################################################################################################
#### EVALUATE SEARCH RESULTS
####################################################################################################
FIND_PACKAGE_HANDLE_STANDARD_ARGS( UMFPACK DEFAULT_MSG
	UMFPACK_LIBRARY
	UMFPACK_INCLUDE_DIR
)

IF( UMFPACK_FOUND )
	SET( UMFPACK_INCLUDE_DIRS
		"${AMD_INCLUDE_DIRS}"
		"${UMFPACK_INCLUDE_DIR}"
	)
	IF( NOT BLAS_FROM_IPOPT )
		SET( UMFPACK_LIBRARIES
			"${UMFPACK_LIBRARY}"
			"${AMD_LIBRARIES}"
			"${BLAS_LIBRARIES}"
	        "${CHOLMOD_LIB}"
		)
	ELSE()
		SET( UMFPACK_LIBRARIES
			"${UMFPACK_LIBRARY}"
			"${AMD_LIBRARIES}"
	        "${CHOLMOD_LIB}"
		)
	ENDIF()
	MARK_AS_ADVANCED(
		UMFPACK_DIR
		UMFPACK_INCLUDE_DIR
		UMFPACK_LIBRARIES
	)
ENDIF()
