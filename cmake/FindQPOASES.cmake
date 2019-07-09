# Search supplied hint directories first if supplied.
set(QPOASES_ROOT_DIR "$ENV{QPOASES_ROOT_DIR}" CACHE PATH "QPOASES root directory.")
if(QPOASES_ROOT_DIR)
 message("Looking for QPOASES in ${QPOASES_ROOT_DIR}")
else(QPOASES_ROOT_DIR)
 message("QPOASES_ROOT_DIR not provided.")
 message("Looking for QPOASES in ${PROJECT_SOURCE_DIR}/ThirdParty/qpOASES-${QPOASES_VERSION}/include")
endif(QPOASES_ROOT_DIR)

find_path(QPOASES_INCLUDE_DIR
  NAMES qpOASES.hpp
  HINTS ${PROJECT_SOURCE_DIR}/ThirdParty/qpOASES-${QPOASES_VERSION}/include
  HINTS /usr/local/include
  HINTS /usr/include
  HINTS ${THIRDPARTY_INSTALL_PATH}/include
)


if(APPLE)
find_library(QPOASES_LIBRARY 
  libqpOASES.dylib
  HINTS /usr/local/lib
  HINTS ${PROJECT_SOURCE_DIR}/ThirdParty/qpOASES-${QPOASES_VERSION}/bin
  HINTS ${QPOASES_ROOT_DIR}/bin
)
elseif(UNIX)
find_library(QPOASES_LIBRARY 
  libqpOASES.so
  HINTS /usr/local/lib
  HINTS ${PROJECT_SOURCE_DIR}/ThirdParty/qpOASES-${QPOASES_VERSION}/bin
  HINTS ${QPOASES_ROOT_DIR}/bin
)
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(QPOASES DEFAULT_MSG QPOASES_LIBRARY QPOASES_INCLUDE_DIR)

if(QPOASES_FOUND)  
  message("—- Found QPOASES include under ${QPOASES_INCLUDE_DIR}")
    set(QPOASES_INCLUDE_DIRS ${QPOASES_INCLUDE_DIR})
    set(QPOASES_LIBRARIES ${QPOASES_LIBRARY})
    if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
        set(QPOASES_LIBRARIES "${QPOASES_LIBRARIES};m;pthread")
    endif(CMAKE_SYSTEM_NAME STREQUAL "Linux")
endif(QPOASES_FOUND)

mark_as_advanced(QPOASES_LIBRARIES QPOASES_INCLUDE_DIR)
