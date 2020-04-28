SET ( COIN_DIR "$ENV{COINOR_ROOT_DIR}" CACHE STRING "Path to an COIN-OR installation.")

MESSAGE( "Looking for Ipopt in ${COIN_DIR}" )

FIND_PACKAGE( PkgConfig REQUIRED )
SET( PKG_CONFIG_USE_CMAKE_PREFIX_PATH 1 )
SET ( CMAKE_PREFIX_PATH ${COIN_DIR}/lib/pkgconfig )

### Look for the Ipopt libraries
pkg_check_modules(IPOPT ipopt)
IF( NOT IPOPT_FOUND )
  MESSAGE( "Ipopt package not found. Set the environment variable COINOR_ROOT_DIR to the base directory of an COIN-OR installation")
ENDIF()

### Now the HSL dependencies.  This is required because QORE also needs MA57
pkg_check_modules(COINHSL REQUIRED coinhsl)

### And if compiled, we probably also need metis and mumps
pkg_check_modules(COINMETIS coinmetis)
pkg_check_modules(COINMUMPS coinmumps)

SET( IPOPT_LIBRARIES ${IPOPT_LINK_LIBRARIES} ${COINHSL_LINK_LIBRARIES} ${COINMETIS_LINK_LIBRARIES} ${COINMUMPS_LINK_LIBRARIES} )

SET( IPOPT_INCLUDE_DIRS ${IPOPT_INCLUDE_DIRS} ${IPOPT_INCLUDE_DIRS}/.. ) 

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(IPOPT DEFAULT_MSG IPOPT_LIBRARIES IPOPT_INCLUDE_DIRS)

### We also want to link with the AMPL Solver Library for the AMPL solvers

### There are the AMPL related Ipopt classes
pkg_check_modules(IPOPTAMPLINTERFACE REQUIRED ipoptamplinterface)

### This is the ASL library provided by AMPL
pkg_check_modules(COINASL REQUIRED coinasl)

SET( AMPLINTERFACE_LIBRARIES ${IPOPTAMPLINTERFACE_LINK_LIBRARIES} ${COINASL_LINK_LIBRARIES} )

SET( AMPLINTERFACE_INCLUDE_DIRS ${IPOPTAMPLINTERFACE_INCLUDE_DIRS} ${COINASL_INCLUDE_DIRS}/.. ) 

find_package_handle_standard_args(AMPLINTERFACE DEFAULT_MSG AMPLINTERFACE_LIBRARIES AMPLINTERFACE_INCLUDE_DIRS)



if (true)

ELSE ()
	# Find IPOPT without relying on pkg-config

	find_path(IPOPT_INCLUDE_DIR
		NAMES IpNLP.hpp 
		HINTS ${PROJECT_SOURCE_DIR}/third_party/CoinIpopt/build/include/coin
		HINTS /usr/local/include/coin
		HINTS ${IPOPT_ROOT_DIR}/include/coin
		HINTS ${IPOPT_ROOT_DIR}/include
	)

	IF( APPLE )
		find_library(IPOPT_LIBRARY 
			libipopt.dylib
			HINTS /usr/local/lib
			HINTS ${PROJECT_SOURCE_DIR}/third_party/CoinIpopt/build/lib
			HINTS ${IPOPT_ROOT_DIR}/lib
		)

		find_library(IPOPT_LIBRARY2 
			libcoinasl.dylib
			HINTS /usr/local/lib
			HINTS ${PROJECT_SOURCE_DIR}/third_party/CoinIpopt/build/lib
			HINTS ${IPOPT_ROOT_DIR}/lib
		)

		find_library(IPOPT_LIBRARY3 
			libipoptamplinterface.dylib
			HINTS /usr/local/lib
			HINTS ${PROJECT_SOURCE_DIR}/third_party/CoinIpopt/build/lib
			HINTS ${IPOPT_ROOT_DIR}/lib
		)

		find_library(IPOPT_LIBRARY4 
			libcoinhsl.dylib
			HINTS /usr/local/lib
			HINTS ${PROJECT_SOURCE_DIR}/third_party/CoinIpopt/build/lib
			HINTS ${IPOPT_ROOT_DIR}/lib
		)

		# find_library(LAPACK_LIB 
		#         liblapack.dylib
		#         HINTS /usr/local/lib
		#         HINTS ${IPOPT_ROOT_DIR}/lib
		# )
		# find_library(BLAS_LIB 
		#         libblas.dylib
		#         HINTS /usr/local/lib
		#         HINTS ${IPOPT_ROOT_DIR}/lib
		# )

	ELSEIF( UNIX )

		find_library(IPOPT_LIBRARY 
			libipopt.so
			HINTS /usr/local/lib
			HINTS ${PROJECT_SOURCE_DIR}/third_party/CoinIpopt/build/lib
			HINTS ${IPOPT_ROOT_DIR}/lib
		)
		find_library(IPOPT_LIBRARY2 
			libcoinasl.so
			HINTS /usr/local/lib
			HINTS ${PROJECT_SOURCE_DIR}/third_party/CoinIpopt/build/lib
			HINTS ${IPOPT_ROOT_DIR}/lib
		)
		find_library(IPOPT_LIBRARY3 
			libipoptamplinterface.so
			HINTS /usr/local/lib
			HINTS ${PROJECT_SOURCE_DIR}/third_party/CoinIpopt/build/lib
			HINTS ${IPOPT_ROOT_DIR}/lib
		)
		find_library(LAPACK_LIB 
		        liblapack.so
		        HINTS /usr/lib
		        HINTS /usr/local/lib
		        HINTS ${IPOPT_ROOT_DIR}/lib
		)
		find_library(BLAS_LIB 
		        libblas.so
		        HINTS /usr/lib
		        HINTS /usr/local/lib
		        HINTS ${IPOPT_ROOT_DIR}/lib
		)

	ENDIF()

	include(FindPackageHandleStandardArgs)
	find_package_handle_standard_args(IPOPT DEFAULT_MSG IPOPT_LIBRARY IPOPT_INCLUDE_DIR)

	if(IPOPT_FOUND)
		message("—- Found Ipopt under ${IPOPT_INCLUDE_DIR}")
		message(${LAPACK_LIB} ${BLAS_LIB})
	    set(IPOPT_INCLUDE_DIRS ${IPOPT_INCLUDE_DIR})
	    set(IPOPT_LIBRARIES ${IPOPT_LIBRARY} ${IPOPT_LIBRARY2} ${IPOPT_LIBRARY3} ${IPOPT_LIBRARY4} ) # ${LAPACK_LIB} ${BLAS_LIB})
	    if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
	        set(IPOPT_LIBRARIES "${IPOPT_LIBRARIES};m;pthread")
	    endif(CMAKE_SYSTEM_NAME STREQUAL "Linux")
	endif(IPOPT_FOUND)

endif()

MARK_AS_ADVANCED( IPOPT_LDFLAGS IPOPT_LIBRARIES IPOPT_INCLUDE_DIRS )
