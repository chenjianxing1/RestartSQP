include_directories(${PROJECT_INCLUDE_DIR})



file(GLOB TEST_SOURCES ${PROJECT_SOURCE_DIR}/test/*.cpp)

add_executable(simple_test ${TEST_SOURCES})

set(IPOPTDIR $ENV{IPOPTDIR})

find_library(AMPLINTERFACE_LIBRARY 
        ipoptamplinterface
        HINTS ${IPOPTDIR}/lib
)

find_library(ASL_LIBRARY 
        coinasl
        HINTS ${IPOPTDIR}/lib
)

#set(AMPLSOLVER_DIR $ENV{IPOPTDIR}/Ipopt/src/Apps/AmplSolver)
#file(GLOB AMPLTNLP_OBJ ${AMPLSOLVER_DIR}/AmplTNLP.o ${AMPLSOLVER_DIR}/.libs/AmplTNLP.o)
#message("AMPLTNLP_OBJ = ${AMPLTNLP_OBJ}")

string(REPLACE ";" " " MY_LD_FLAGS "${IPOPT_LDFLAGS}")
#message(${LD_FLAGS})
##set(CMAKE_EXE_LINKER_FLAGS "${LD_FLAGS}")

#target_link_libraries(simple_test sqphotstart ${IPOPT_LDFLAGS} ${AMPLTNLP_OBJ})
target_link_libraries(simple_test sqphotstart ${AMPLINTERFACE_LIBRARY} ${ASL_LIBRARY} ${MY_LD_FLAGS})

