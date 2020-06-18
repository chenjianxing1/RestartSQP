cmake_minimum_required(VERSION 3.2)

project(RESTARTSQP CXX C)

# There are some C++-11 elements in the code
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)


# CK: This avoids the use of __cxx11 and __1 inline namespaces with STL data types
if(APPLE)
  ADD_DEFINITIONS(-D_GLIBCXX_USE_CXX11_ABI=0)
ENDIF()

LIST( APPEND CMAKE_MODULE_PATH
  ${PROJECT_SOURCE_DIR}/cmake
)

set(MACOSX_RPATH 1)

#######################################################################
## Options

option(Debug "Turn on debug mode" OFF)

option(Qore "Link to QORE libraries" ON)

option(qpOASES "Link to qpOASES libraries" ON)

########################################################################
## Get QORE
if(Qore)
  message(STATUS "Enable QORE")
  add_definitions(-DUSE_QORE)
  find_package(QORE REQUIRED)
  message("-- Found QORE: ${QORE_INCLUDE_DIRS}")
  include_directories(${QORE_INCLUDE_DIRS})
  set(DEPENDENCY_LIBS ${DEPENDENCY_LIBS} ${QORE_LIBRARIES})
  foreach (dir ${QORE_INCLUDE_DIRS})
    set(HEADER_FILES ${HEADER_FILES} ${dir}/*.h)
  endforeach(dir)
  MESSAGE(STATUS "QORE_LIBRARIES = ${QORE_LIBRARIES}")
endif(Qore)

########################################################################
## Get qpOASES
if (qpOASES)
  message(STATUS "Enable qpOASES")
  add_definitions(-DUSE_QPOASES)
  find_package(QPOASES REQUIRED)
  message("-- Found QPOASES: ${QPOASES_INCLUDE_DIR}")
  include_directories(${QPOASES_INCLUDE_DIR} ${QPOASES_INCLUDE_DIR}/qpOASES)
  set(DEPENDENCY_LIBS ${DEPENDENCY_LIBS} ${QPOASES_LIBRARY})
  set(HEADER_FILES ${HEADER_FILES} ${QPOASES_INCLUDE_DIR}/*.h* ${QPOASES_INCLUDE_DIR}/qpOASES/*.h*)
endif ()

#########################################################################
# Find Ipopt.
message(STATUS "Enable IPOPT")
find_package(IPOPT REQUIRED)
include_directories(${IPOPT_INCLUDE_DIRS})
foreach (dir ${IPOPT_INCLUDE_DIRS})
  set(HEADER_FILES ${HEADER_FILES} ${dir}/*.h*)
endforeach(dir)
set(DEPENDENCY_LIBS ${DEPENDENCY_LIBS} ${IPOPT_LIBRARIES})

#########################################################################
# Find LAPACK and BLAS
find_package(LAPACK REQUIRED)
find_package(BLAS REQUIRED)

set(DEPENDENCY_LIBS ${DEPENDENCY_LIBS} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} -lstdc++)

#########################################################################
## Set include directories
set(PROJECT_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include)
include_directories(${PROJECT_INCLUDE_DIR})

## Get the list of all header files from dependencies
file(GLOB DEPEND_HEADERS ${HEADER_FILES})

message(STATUS "DEPEND_HEADERS ${DEPEND_HEADERS}")




## Build options
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
##set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin/Release)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

#set(CMAKE_BUILD_TYPE "Release")
#if (Debug)
#set(CMAKE_BUILD_TYPE "Debug")
#set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin/Debug)
#endif()

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -g")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,-rpath,'\$ORIGIN/lib'")


if(LINUX)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")
endif()

add_subdirectory(src)
add_subdirectory(interfaces/AmplInterface)
add_subdirectory(test)
