include_directories(${PROJECT_INCLUDE_DIR})



# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
#SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)


# the RPATH to be used when installing, but only if it's not a system directory
#SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib:'$ORIGIN/lib'")

file(GLOB HEADERS ${PROJECT_SOURCE_DIR}/src/*.hpp)

file(GLOB SQPHOTSTART_SRC ${PROJECT_SOURCE_DIR}/src/*.cpp)

add_executable(SQPHOTSTART ${SQPHOTSTART_SRC} ${HEADERS})
