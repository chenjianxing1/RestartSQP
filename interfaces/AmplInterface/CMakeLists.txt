include_directories(${PROJECT_INCLUDE_DIR})

file(GLOB AMPL_SOURCES ${PROJECT_SOURCE_DIR}/interfaces/AmplInterface/*.cpp)
file(GLOB AMPL_HEADERS ${PROJECT_SOURCE_DIR}/interfaces/AmplInterface/*.hpp)

message(STATUS "AMPLINTERFACE_INCLUDE_DIRS = ${AMPLINTERFACE_INCLUDE_DIRS}")

add_executable(crossoversqp_ampl ${AMPL_SOURCES} ${AMPL_HEADERS})

include_directories(${QPOASES_INCLUDE_DIRS})
include_directories(${QORE_INCLUDE_DIR})
include_directories(${AMPLINTERFACE_INCLUDE_DIRS})

target_link_libraries(crossoversqp_ampl restartsqp ${AMPLINTERFACE_LIBRARIES} ${DEPENDENCY_LIBS})


