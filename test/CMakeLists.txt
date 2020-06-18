include_directories(${PROJECT_INCLUDE_DIR})

file(GLOB TEST_SOURCES ${PROJECT_SOURCE_DIR}/test/hs071_lazy/hs071_main.cpp)
file(GLOB TEST_SOURCES ${PROJECT_SOURCE_DIR}/test/hs071_lazy/hs071_nlp.cpp)
file(GLOB TEST_SOURCES ${PROJECT_SOURCE_DIR}/test/sqp_ampl.cpp)
file(GLOB TEST_SOURCES ${PROJECT_SOURCE_DIR}/test/lazy_sqp_ampl.cpp)
file(GLOB TEST_SOURCES ${PROJECT_SOURCE_DIR}/test/hs071_c/hs071_c.c)

add_executable(hs071_lazy ${PROJECT_SOURCE_DIR}/test/hs071_lazy/hs071_main.cpp ${PROJECT_SOURCE_DIR}/test/hs071_lazy/hs071_nlp.cpp)
add_executable(QPsolvers_test ${PROJECT_SOURCE_DIR}/test/QPsolvers_testers.cpp)
add_executable(unitTest_SpHbMat ${PROJECT_SOURCE_DIR}/test/unitTest/test_SpHbMat.cpp)
add_executable(unitTest_SpTripletMat ${PROJECT_SOURCE_DIR}/test/unitTest/test_SpTripletMat.cpp)
add_executable(sqp_ampl ${PROJECT_SOURCE_DIR}/test/sqp_ampl.cpp)
add_executable(lazy_sqp_ampl ${PROJECT_SOURCE_DIR}/test/lazy_sqp_ampl.cpp)
add_executable(hs071_c ${PROJECT_SOURCE_DIR}/test/hs071_c/hs071_c.c)


include_directories(${PROJECT_SOURCE_DIR}/test/unitTest)
include_directories(${PROJECT_SOURCE_DIR}/test/hs071_lazy)
include_directories(${QPOASES_INCLUDE_DIRS})
include_directories(${QORE_INCLUDE_DIR})


message( "DEPENDENCY_LIBS = ${DEPENDENCY_LIBS}") 
target_link_libraries(sqp_ampl restartsqp ${AMPLINTERFACE_LIBRARIES} ${DEPENDENCY_LIBS})
target_link_libraries(lazy_sqp_ampl restartsqp ${AMPLINTERFACE_LIBRARIES} ${DEPENDENCY_LIBS})
target_link_libraries(hs071_lazy restartsqp ${DEPENDENCY_LIBS})
target_link_libraries(QPsolvers_test restartsqp ${DEPENDENCY_LIBS})
target_link_libraries(unitTest_SpHbMat restartsqp ${DEPENDENCY_LIBS})
target_link_libraries(unitTest_SpTripletMat restartsqp ${DEPENDENCY_LIBS})
target_link_libraries(hs071_c restartsqp ${DEPENDENCY_LIBS})


