# Create download URL derived from version number.
set(QPOASES_HOME https://www.qpoases.org/go/release)
set(QPOASES_DOWNLOAD_URL ${QPOASES_HOME})
unset(QPOASES_HOME)

# Download and build the QPOASES library
set(QPOASES_ROOT ${THIRDPARTY_INSTALL_PATH} CACHE INTERNAL "")
ExternalProject_Add(QPOASES
    DOWNLOAD_DIR ${THIRDPARTY_INSTALL_PATH}
    DOWNLOAD_COMMAND export HTTPS_PROXY=$ENV{HTTPS_PROXY} && curl -k -L ${QPOASES_DOWNLOAD_URL} -o QPOASES.zip && unzip QPOASES.zip && mv qpOASES-${QPOASES_VERSION} QPOASES && rm -fr ./Install/QPOASES && mv QPOASES ./Install && cd ./Install/QPOASES && make src && cd bin && ranlib libqpOASES.a
    URL ${QPOASES_DOWNLOAD_URL}
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${QPOASES_ROOT}
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
)

#list(APPEND QPOASES_LIBRARIES "${QPOASES_ROOT}/Install/QPOASES/bin/libqpOASES.a")

#find_library(QPOASES_LIBRARY libqpOASES.a HINTS ${THIRDPARTY_INSTALL_PATH}/Install/QPOASES/bin)
#mark_as_advanced(QPOASES_LIBRARY)
#set(QPOASES_LIBRARIES ${QPOASES_LIBRARY})
#message("-- Found QPOASES libraries ${QPOASES_LIBRARIES}")
mark_as_advanced(QPOASES_ROOT)
# set(LIBS ${LIBS} ${QPOASES_LIBRARIES})
# unset(QPOASES_DOWNLOAD_URL)
# unset(QPOASES_ROOT)
