include_directories(${PROJECT_INCLUDE_DIR})

file(GLOB HEADERS *.hpp ${PROJECT_INCLUDE_DIR}/sqphot/*.hpp)
file(GLOB SOURCES *.cpp ${PROJECT_SRC_DIR}/src/*.cpp)

if (NOT Gurobi)
  list(REMOVE_ITEM SOURCES ${PROJECT_SOURCE_DIR}/src/GurobiInterface.cpp)
endif()

if (NOT Cplex)
    list(REMOVE_ITEM SOURCES ${PROJECT_SOURCE_DIR}/src/CplexInterface.cpp)
endif()


add_library(restartsqp STATIC ${SOURCES} ${HEADERS})

if(ADD_EIGEN)
 add_dependencies(restartsqp eigen)
endif(ADD_EIGEN)

if(ADD_QPOASES)
 add_dependencies(restartsqp QPOASES)
endif(ADD_QPOASES)

target_link_libraries(restartsqp ${CPLEX_LIBRARIES} ${GUROBI_LIBRARIES} ${QPOASES_LIBRARIES} ${QORE_LIBRARIES}) # ${IPOPT_LDFLAGS})
#set(CMAKE_EXE_LINKER_FLAGS "${IPOPT_LDFLAGS}")

if(Cov)
  set_target_properties(restartsqp
  PROPERTIES
    COMPILE_FLAGS "-ftest-coverage -fprofile-arcs"
    LINK_FLAGS    "-lgcov --coverage -fprofile-arcs"
  )
endif()

SET(CMAKE_SKIP_BUILD_RPATH  FALSE)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE) 
SET(CMAKE_INSTALL_RPATH "${PROJ_LIB_DIR}")
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)


file(GLOB SQPHOTSTART_SRC ${PROJECT_SOURCE_DIR}/src/*.cpp)

