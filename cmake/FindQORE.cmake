set(QORE_SRC_DIR "$ENV{QORE_SRC_DIR}" CACHE PATH "QORE source directory.")
set(QORE_BIN_DIR "$ENV{QORE_BIN_DIR}" CACHE PATH "QORE binary directory.")
message("Looking for QORE in ${QORE_SRC_DIR} and ${QORE_BIN_DIR}")

###########################################################################
## Look for the QORE library itself
###########################################################################

find_path(QORE_INCLUDE_DIR2
  NAMES qpsolver.h 
  HINTS /usr/local/include
  HINTS ${QORE_SRC_DIR}/QPSOLVER/include
  )

find_path(QPPRESOLVER_INCLUDE_DIR
  NAMES qpPresolver.h 
  HINTS /usr/local/include
  HINTS ${QORE_SRC_DIR}/qpPresolver/include
)

find_path(QORE_INCLUDE_DIR3
  NAMES qp_types.h 
  HINTS /usr/local/include
  HINTS ${QORE_SRC_DIR}
)

set(QORE_INCLUDE_DIRS
  ${QORE_INCLUDE_DIR3}
  ${QORE_INCLUDE_DIR2}
  ${QPPRESOLVER_INCLUDE_DIR}
  )

find_library(QORE_LIBRARY 
  qore
  HINTS /usr/local/lib
  HINTS ${QORE_BIN_DIR}/lib
  )

set(QORE_LIBRARIES ${QORE_LIBRARY}) 

###########################################################################
## Now look for QORE dependency UMFPACK

FIND_PACKAGE( UMFPACK REQUIRED )

## Add all the libraries and include directories required by UMFPACK to QORE

SET( QORE_LIBRARIES ${QORE_LIBRARIES} ${UMFPACK_LIBRARIES} )

SET( QORE_INCLUDE_DIRS ${QORE_INCLUDE_DIRS} ${UMFPACK_INCLUDE_DIRS} )

message("QORE_LIBRARIES = ${QORE_LIBRARIES}")
message("QORE_INCLUDE_DIRS = ${QORE_INCLUDE_DIRS}")

INCLUDE( FindPackageHandleStandardArgs )
find_package_handle_standard_args(QORE DEFAULT_MSG QORE_LIBRARIES QORE_INCLUDE_DIRS)

#mark_as_advanced(QORE_LIBRARIES QORE_INCLUDE_DIR)
