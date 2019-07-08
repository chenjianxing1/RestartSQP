# Search supplied hint directories first if supplied.
set(QPOASES_ROOT_DIR "$ENV{QPOASES_ROOT_DIR}" CACHE PATH "QPOASES root directory.")
if(QPOASES_ROOT_DIR)
 message("Looking for QPOASES in ${QPOASES_ROOT_DIR}")
else(QPOASES_ROOT_DIR)
 message("QPOASES_ROOT_DIR not provided.")
endif(QPOASES_ROOT_DIR)

find_path(QPOASES_INCLUDE_DIR
  NAMES QPOASES/Core
  HINTS /usr/local/include
  HINTS /usr/include
  HINTS ${THIRDPARTY_INSTALL_PATH}/include)

# Extract QPOASES version from QPOASES/src/Core/util/Macros.h
if (QPOASES_INCLUDE_DIR)
  set(QPOASES_VERSION_FILE ${QPOASES_INCLUDE_DIR}/QPOASES/src/Core/util/Macros.h)
  if (NOT EXISTS ${QPOASES_VERSION_FILE})
    QPOASES_report_not_found(
      "Could not find file: ${QPOASES_VERSION_FILE} "
      "containing version information in QPOASES install located at: "
      "${QPOASES_INCLUDE_DIR}.")
  else (NOT EXISTS ${QPOASES_VERSION_FILE})
    file(READ ${QPOASES_VERSION_FILE} QPOASES_VERSION_FILE_CONTENTS)

    string(REGEX MATCH "#define QPOASES_WORLD_VERSION [0-9]+"
      QPOASES_WORLD_VERSION "${QPOASES_VERSION_FILE_CONTENTS}")
    string(REGEX REPLACE "#define QPOASES_WORLD_VERSION ([0-9]+)" "\\1"
      QPOASES_WORLD_VERSION "${QPOASES_WORLD_VERSION}")

    string(REGEX MATCH "#define QPOASES_MAJOR_VERSION [0-9]+"
      QPOASES_MAJOR_VERSION "${QPOASES_VERSION_FILE_CONTENTS}")
    string(REGEX REPLACE "#define QPOASES_MAJOR_VERSION ([0-9]+)" "\\1"
      QPOASES_MAJOR_VERSION "${QPOASES_MAJOR_VERSION}")

    string(REGEX MATCH "#define QPOASES_MINOR_VERSION [0-9]+"
      QPOASES_MINOR_VERSION "${QPOASES_VERSION_FILE_CONTENTS}")
    string(REGEX REPLACE "#define QPOASES_MINOR_VERSION ([0-9]+)" "\\1"
      QPOASES_MINOR_VERSION "${QPOASES_MINOR_VERSION}")

    # This is on a single line s/t CMake does not interpret it as a list of
    # elements and insert ';' separators which would result in 3.;2.;0 nonsense.
    set(QPOASES_VERSION "${QPOASES_WORLD_VERSION}.${QPOASES_MAJOR_VERSION}.${QPOASES_MINOR_VERSION}")
  endif (NOT EXISTS ${QPOASES_VERSION_FILE})
endif (QPOASES_INCLUDE_DIR)

# Set standard CMake FindPackage variables if found.
if (QPOASES_FOUND)
  set(QPOASES_INCLUDE_DIRS ${QPOASES_INCLUDE_DIR})
else (QPOASES_FOUND)
 message("Cannot find QPOASES, will try pulling it from github.")
endif (QPOASES_FOUND)
