include_directories(${PROJECT_INCLUDE_DIR})



file(GLOB TEST_SOURCES ${PROJECT_SOURCE_DIR}/test/*.cpp)

add_executable(simple_test ${TEST_SOURCES})
target_link_libraries(simple_test sqphotstart)
