SET ( COIN_DIR "$ENV{COINOR_ROOT_DIR}" CACHE STRING "Path to an COIN-OR installation.")

MESSAGE( "Looking for Ipopt in ${COIN_DIR}" )

set(ENV{PKG_CONFIG_PATH} "${COIN_DIR}/lib/pkgconfig")

FIND_PACKAGE( PkgConfig REQUIRED )
SET( PKG_CONFIG_USE_CMAKE_PREFIX_PATH 1 )
# THE FOLLOWING DOES NOT SEEM TO WORK!!!!
# pkg_config takes the environment variable
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

SET( AMPLINTERFACE_INCLUDE_DIRS ${IPOPTAMPLINTERFACE_INCLUDE_DIRS} ${COINASL_INCLUDE_DIRS} ) 

find_package_handle_standard_args(AMPLINTERFACE DEFAULT_MSG AMPLINTERFACE_LIBRARIES AMPLINTERFACE_INCLUDE_DIRS)


MARK_AS_ADVANCED( IPOPT_LDFLAGS IPOPT_LIBRARIES IPOPT_INCLUDE_DIRS )
