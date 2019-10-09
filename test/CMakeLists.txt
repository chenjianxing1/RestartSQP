include_directories(${PROJECT_INCLUDE_DIR})



file(GLOB TEST_SOURCES ${PROJECT_SOURCE_DIR}/test/*.cpp)

add_executable(simple_test ${TEST_SOURCES})
add_executable(unitTest_SpHbMatrix ${PROJECT_SOURCE_DIR}/test/unitTest/test_SpHbMat.cpp)

set(IPOPTDIR $ENV{IPOPTDIR})

find_library(AMPLINTERFACE_LIBRARY 
        ipoptamplinterface
        HINTS ${IPOPTDIR}/lib
)

find_library(ASL_LIBRARY 
        coinasl
        HINTS ${IPOPTDIR}/lib
)

#string(REPLACE ";" " " MY_LD_FLAGS "${IPOPT_LDFLAGS}")
#message(${LD_FLAGS})
##set(CMAKE_EXE_LINKER_FLAGS "${LD_FLAGS}")

#set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${IPOPT_LINK_FLAGS}")
#message("LIBS = ${LIBS}")
#message("IPOPT_LDFLAGS = ${IPOPT_LDFLAGS}")
target_link_libraries(simple_test sqphotstart ${AMPLINTERFACE_LIBRARY} ${ASL_LIBRARY} ${LIBS})
target_link_libraries(unitTest_SpHbMatrix sqphotstart ${AMPLINTERFACE_LIBRARY} ${ASL_LIBRARY} ${LIBS})


