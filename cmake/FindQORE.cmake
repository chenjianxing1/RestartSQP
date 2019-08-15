set(QORE_ROOT_DIR "$ENV{QORE_ROOT_DIR}" CACHE PATH "IPOPT root directory.")
message("Looking for QORE in ${QORE_ROOT_DIR}")


find_path(QORE_INCLUDE_DIR2
	NAMES qpsolver.h 
	HINTS /usr/local/include
	HINTS ${QORE_ROOT_DIR}
	)
MESSAGE(STATUS "QORE_INCLUDE_DIR2 = ${QORE_INCLUDE_DIR2}") 
find_path(QPPRESOLVER_INCLUDE_DIR
	NAMES qpPresolver.h 
	HINTS /usr/local/include
	HINTS ${QORE_ROOT_DIR}/qpPresolver/include
)

find_path(QORE_INCLUDE_DIR3
	NAMES qp_types.h 
	HINTS /usr/local/include
	HINTS ${QORE_ROOT_DIR}
)
 set(QORE_INCLUDE_DIR ${QORE_INCLUDE_DIR3} ${QORE_INCLUDE_DIR2} ${QPPRESOLVER_INCLUDE_DIR})

if(APPLE)
find_library(QORE_LIBRARY 
	libqore.dylib
	HINTS /usr/local/lib
	HINTS third_party/QORE/lib
	HINTS ${QORE_ROOT_DIR}/lib
)
find_library(QORE_LIBRARY2 
	libblas.dylib
	HINTS /usr/local/lib
	HINTS ${PROJECT_SOURCE_DIR}/third_party/CoinIpopt/build/lib
	HINTS ${QORE_ROOT_DIR}/lib
)

find_library(QORE_LIBRARY3 
	liblapack.dylib
	HINTS /usr/local/lib
	HINTS ${PROJECT_SOURCE_DIR}/third_party/CoinIpopt/build/lib
	HINTS ${QORE_ROOT_DIR}/lib
)

elseif(UNIX)
find_library(QORE_LIBRARY 
	libqore.so
	HINTS /usr/local/lib
	HINTS third_party/QORE/lib
	HINTS ${QORE_ROOT_DIR}/lib
)
find_library(QORE_LIBRARY2 
	libblas.so
	HINTS /usr/local/lib
	HINTS ${PROJECT_SOURCE_DIR}/third_party/CoinIpopt/build/lib
	HINTS ${QORE_ROOT_DIR}/lib
)

find_library(QORE_LIBRARY3 
	liblapack.so
	HINTS /usr/local/lib
	HINTS ${PROJECT_SOURCE_DIR}/third_party/CoinIpopt/build/lib
	HINTS ${QORE_ROOT_DIR}/lib
)


endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(QORE DEFAULT_MSG QORE_LIBRARY QORE_INCLUDE_DIR)

if(QORE_FOUND)
	message("—- Found QORE QORE under ${QORE_INCLUDE_DIR}")
    set(QORE_INCLUDE_DIRS ${QORE_INCLUDE_DIR})
    set(QORE_LIBRARIES ${QORE_LIBRARY}) 
    if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
        set(QORE_LIBRARIES "${QORE_LIBRARIES};m;pthread")
    endif(CMAKE_SYSTEM_NAME STREQUAL "Linux")
endif(QORE_FOUND)

mark_as_advanced(QORE_LIBRARIES QORE_INCLUDE_DIR)
