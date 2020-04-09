cmake_minimum_required(VERSION 3.2)

project(RESTARTSQP CXX)

# CK: This avoids the use of __cxx11 and __1 inline namespaces with STL data types
IF (APPLE)
  ADD_DEFINITIONS(-D_GLIBCXX_USE_CXX11_ABI=0)
ENDIF()

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake)
set(MACOSX_RPATH 1) 
# Debug mode
option(Debug "Turn on debug mode" OFF)

# Coverage
option(Cov "Turn on coverage" OFF)

# Find GUROBI (optional)
option(Gurobi "Link to GUROBI libraries" OFF)

# Find IPOPT (optional)
option(Ipopt "Link to IPOPT libraries" ON)

# Find CPLEX (optional)
option(Cplex "Link to CPLEX libraries" OFF)


option(Qore "Link to QORE libraries" ON)

# Find gurobi?
if (Gurobi)
    message(STATUS "Enable GUROBI")
    add_definitions(-DUSE_GUROBI)
    find_package(GUROBI)
    include_directories(${GUROBI_INCLUDE_DIRS})
    set(LIBS ${LIBS} ${GUROBI_LIBRARIES})
endif()

# Find CPLEX?
if (Cplex)
    message(STATUS "Enable CPLEX")
    add_definitions(-DUSE_CPLEX)
    find_package(CPLEX)
    include_directories(${CPLEX_INCLUDE_DIRS})
    set(LIBS ${LIBS} ${CPLEX_LIBRARIES})
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -DNDEBUG -DIL_STD")
endif()


set_directory_properties(PROPERTIES EP_BASE ${CMAKE_BINARY_DIR}/third_party)
get_directory_property(THIRDPARTY_BASE_PATH EP_BASE)

set(THIRDPARTY_INSTALL_PATH ${THIRDPARTY_BASE_PATH})


include(ExternalProject)
# Eigen
option(Eigen "Turn on Eigen" OFF)
if (Eigen)
  find_package(Eigen QUIET)
  if (EIGEN_FOUND)
  message("-- Found Eigen version ${EIGEN_VERSION}: ${EIGEN_INCLUDE_DIRS}")
  else(EIGEN_FOUND)
   set(EIGEN_VERSION 3.3.7)
    include(ExternalEigen)
    set(EIGEN_INCLUDE_DIR ${THIRDPARTY_INSTALL_PATH}/Install/eigen)
    set(ADD_EIGEN TRUE)
  endif (EIGEN_FOUND)
  include_directories(${EIGEN_INCLUDE_DIR})
endif()


#QORE
if(Qore)
	message(STATUS "Enable QORE")
	add_definitions(-DUSE_QORE)
	find_package(QORE)
	if(QORE_FOUND)
		message("-- Found QORE: ${QORE_INCLUDE_DIRS}") 
		set(LIBS ${LIBS} ${QORE_LIBRARIES})
		MESSAGE(STATUS "QORE_LIBRARIES = ${QORE_LIBRARIES}")
	endif (QORE_FOUND)
#INCLUDE ( FindUMFPACK )
FIND_PACKAGE ( UMFPACK REQUIRED )
if (UMFPACK_FOUND)
     set(LIBS ${LIBS} ${UMFPACK_LIBRARIES})
     MESSAGE(STATUS "UMFPACK_LIBRARIES = ${UMFPACK_LIBRARIES}")
endif (UMFPACK_FOUND)


endif(Qore)

# QPOASES
set(QPOASES_VERSION 3.2.1)
find_package(QPOASES)
if (QPOASES_FOUND)
 message("-- Found QPOASES version ${QPOASES_VERSION}: ${QPOASES_INCLUDE_DIRS}") 
 set(QPOASES_INCLUDE_DIR1 ${QPOASES_INCLUDE_DIR})
 set(QPOASES_INCLUDE_DIR2 ${QPOASES_INCLUDE_DIR}/qpOASES)
 set(LIBS ${LIBS} ${QPOASES_LIBRARIES})
else(QPOASES_FOUND) 
    message("-- Was not able to find QPOASES") 
 include(ExternalQPOASES)
 set(ADD_QPOASES TRUE)
 set(QPOASES_INCLUDE_DIR1 ${THIRDPARTY_INSTALL_PATH}/Install/QPOASES/include/qpOASES)
 set(QPOASES_INCLUDE_DIR2 ${THIRDPARTY_INSTALL_PATH}/Install/QPOASES/include)
endif (QPOASES_FOUND)


include_directories(${QPOASES_INCLUDE_DIR1})
include_directories(${QPOASES_INCLUDE_DIR2})
include_directories(${THIRDPARTY_INSTALL_PATH}/include)
include_directories(${QORE_INCLUDE_DIR})

# Find Ipopt.
if (Ipopt)
	message(STATUS "Enable IPOPT")
	add_definitions(-DUSE_IPOPT)
	find_package(IPOPT REQUIRED)
  include_directories(${IPOPT_INCLUDE_DIRS}/..)
  include_directories(${IPOPT_INCLUDE_DIRS})          # CK added this
	set(LIBS ${LIBS} ${IPOPT_LDFLAGS})
endif()

#set(LIBS ${LIBS} ${PROJECT_BINARY_DIR}/lib)
# Compiler options
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Ofast")

# Build options
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
#set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin/Release)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
set(PROJECT_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include)

#set(CMAKE_BUILD_TYPE "Release")
#if (Debug)
#set(CMAKE_BUILD_TYPE "Debug")
#set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin/Debug)
#endif()

if(Cov)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-arcs -ftest-coverage")
set(CMAKE_LINKER_FLAGS "${CMAKE_LINKER_FLAGS} -fprofile-arcs -lgcov")
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -g")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,-rpath,'\$ORIGIN/lib'")


if(LINUX)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")
endif()


add_subdirectory(src)
add_subdirectory(test)
