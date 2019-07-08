# Create download URL derived from version number.
set(QPOASES_HOME https://github.com/startcode/qp-oases/archive)
set(QPOASES_DOWNLOAD_URL ${QPOASES_HOME}/v${QPOASES_VERSION}.zip)
unset(QPOASES_HOME)

# Download and build the Eigen library and add its properties to the third party arguments.
set(QPOASES_ROOT ${THIRDPARTY_INSTALL_PATH} CACHE INTERNAL "")
ExternalProject_Add(QPOASES
    DOWNLOAD_DIR ${THIRDPARTY_INSTALL_PATH}
    DOWNLOAD_COMMAND export HTTPS_PROXY=$ENV{HTTPS_PROXY} && curl -k -L ${QPOASES_DOWNLOAD_URL} -o QPOASES.zip && unzip QPOASES.zip && mv qp-oases-${QPOASES_VERSION} QPOASES && rm -fr ./Install/QPOASES && mv QPOASES ./Install && cd ./Install/QPOASES && mkdir bin && make src
    URL ${QPOASES_DOWNLOAD_URL}
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${QPOASES_ROOT}
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
)

list(APPEND GLOBAL_THIRDPARTY_LIB_ARGS "-DQPOASES_ROOT:PATH=${QPOASES_ROOT}")
set(QPOASES_LIBRARY ${THIRDPARTY_INSTALL_PATH}/Install/QPOASES/bin)
mark_as_advanced(QPOASES_LIBRARY)
unset(QPOASES_DOWNLOAD_URL)
unset(QPOASES_ROOT)
