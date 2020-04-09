
SET( IPOPT_ROOT_DIR $ENV{IPOPT_ROOT_DIR} CACHE PATH "Path to an Ipopt installation." )
MESSAGE( "Looking for Ipopt in ${IPOPT_ROOT_DIR}" )
OPTION ( IPOPT_USE_PKGCONFIG "Find Ipopt using pkg-config" ON )

IF (IPOPT_USE_PKGCONFIG)
	find_package(PkgConfig REQUIRED)

	# The following does not work, it should automatically set PKG_CONFIG_PATH but it doesn't, so user has to set it manually for now.
	SET(PKG_CONFIG_USE_CMAKE_PREFIX_PATH 1)
	SET(CMAKE_PREFIX_PATH ${IPOPT_ROOT_DIR}/lib/pkgconfig)
	SET(ENV{PKG_CONFIG_PATH} "${CMAKE_PREFIX_PATH}:${PKG_CONFIG_PATH}")
	# message("CMAKE_PREFIX_PATH = ${CMAKE_PREFIX_PATH}")
	# message("*** RIGHT NOW, WE NEED TO SET THE ENVIRONMENT VARIABLE PGK_SEARCH_PATH TO THIS VALUE:")
	# message("Type: export PKG_CONFIG_PATH=${CMAKE_PREFIX_PATH}:$PKG_CONFIG_PATH")
	pkg_check_modules(IPOPT REQUIRED ipopt)
	# message("*** Make sure your environment variable LD_LIBRARY_PATH includes $IPOPT_ROOT_DIR/lib")

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
