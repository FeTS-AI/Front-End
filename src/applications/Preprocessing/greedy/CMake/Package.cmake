# ----------------------------------------------------------------
# INSTALLATION AND PACKAGING with CPack
# ----------------------------------------------------------------

# On Win32, we must include the redistributable
IF(MSVC_VERSION GREATER 1399)
  FIND_PROGRAM(VCREDIST_X86 vcredist_x86.exe)
  IF(VCREDIST_X86)
    INSTALL(FILES ${VCREDIST_X86} DESTINATION bin)
    SET(CPACK_NSIS_EXTRA_INSTALL_COMMANDS 
      "ExecWait '\\\"$INSTDIR\\\\bin\\\\vcredist_x86.exe\\\" /passive'")
  ENDIF(VCREDIST_X86)
ENDIF(MSVC_VERSION GREATER 1399)

# Allow package generation
SET(CPACK_PACKAGE_NAME "greedy")
SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Greedy medical image registration tool")
SET(CPACK_PACKAGE_VENDOR "itksnap.org")
SET(CPACK_PACKAGE_VERSION_MAJOR "${GREEDY_VERSION_MAJOR}")
SET(CPACK_PACKAGE_VERSION_MINOR "${GREEDY_VERSION_MINOR}")
SET(CPACK_PACKAGE_VERSION_PATCH "${GREEDY_VERSION_PATCH}")
SET(CPACK_NSIS_MODIFY_PATH ON)


# Shamelessly stolen from ParaView_
SET(CPACK_SOURCE_PACKAGE_FILE_NAME "greedy-${GREEDY_VERSION_FULL}")
IF (CMAKE_SYSTEM_PROCESSOR MATCHES "unknown")
  EXEC_PROGRAM(uname ARGS "-m" OUTPUT_VARIABLE CMAKE_SYSTEM_PROCESSOR)
ENDIF (CMAKE_SYSTEM_PROCESSOR MATCHES "unknown")
IF(NOT DEFINED CPACK_SYSTEM_NAME)
  SET(CPACK_SYSTEM_NAME ${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR})
ENDIF(NOT DEFINED CPACK_SYSTEM_NAME)
IF(${CPACK_SYSTEM_NAME} MATCHES Windows)
  IF(CMAKE_CL_64)
    SET(CPACK_SYSTEM_NAME win64-${CMAKE_SYSTEM_PROCESSOR})
  ELSE(CMAKE_CL_64)
    SET(CPACK_SYSTEM_NAME win32-${CMAKE_SYSTEM_PROCESSOR})
  ENDIF(CMAKE_CL_64)
ENDIF(${CPACK_SYSTEM_NAME} MATCHES Windows)

# For Apple, we need to base the filename on the architecture
IF(CMAKE_SYSTEM_NAME MATCHES Darwin)
  IF(NOT DEFINED CMAKE_OSX_ARCHITECTURES)
    MESSAGE(ERROR "CMAKE_OSX_ARCHITECTURES must be defined")
  ENDIF(NOT DEFINED CMAKE_OSX_ARCHITECTURES)
  STRING(REPLACE ";" "-" ARCH "${CMAKE_OSX_ARCHITECTURES}")
  SET(CPACK_SYSTEM_NAME "MacOS-${ARCH}")
  MESSAGE(STATUS "   ARCH ${ARCH}")
  MESSAGE(STATUS "   CMAKE_OSX_ARCHITECTURES ${CMAKE_OSX_ARCHITECTURES}")
ENDIF(CMAKE_SYSTEM_NAME MATCHES Darwin)

MESSAGE(STATUS "   CPACK_SYSTEM_NAME ${CPACK_SYSTEM_NAME}")

SET(CPACK_PACKAGE_FILE_NAME "${CPACK_SOURCE_PACKAGE_FILE_NAME}-${CPACK_SYSTEM_NAME}")

MESSAGE(STATUS "   CPACK_PACKAGE_FILE_NAME ${CPACK_PACKAGE_FILE_NAME}")

# Show GPL license
SET(CPACK_RESOURCE_FILE_LICENSE "${GREEDY_TOOL_SOURCE_DIR}/COPYING.txt")

IF(WIN32 AND NOT UNIX)

  SET(CPACK_GENERATOR "NSIS")
  SET(CPACK_EXTENSION "exe")

ELSE(WIN32 AND NOT UNIX)

  # Set the generator to either STGZ or Apple
  IF(NOT APPLE)
    SET(CPACK_DEBIAN_PACKAGE_MAINTAINER "pauly2@mail.med.upenn.edu")
    SET(CPACK_GENERATOR "TGZ")
    SET(CPACK_EXTENSION "tar.gz")
  ELSE(NOT APPLE)
    SET(CPACK_GENERATOR "DragNDrop")
    SET(CPACK_EXTENSION "dmg")
  ENDIF(NOT APPLE)

ENDIF(WIN32 AND NOT UNIX)

#--------------------------------------------------------------------------------
# Construct the name of the package
SET(CPACK_PACKAGE_FILE_NAME_WEXT "${CPACK_PACKAGE_FILE_NAME}.${CPACK_EXTENSION}")


INCLUDE(CPack)

#--------------------------------------------------------------------------------
# Uploading code to SourceForge
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# Configure SCP

FIND_PROGRAM(SCP_PROGRAM NAMES scp DOC "Location of the scp program (optional)")
MARK_AS_ADVANCED(SCP_PROGRAM)

SET(SCP_ARGUMENTS "-v" CACHE STRING "Optional arguments to the scp command for uploads to SourceForge")
MARK_AS_ADVANCED(SCP_ARGUMENTS)

SET(SCP_USERNAME "" CACHE STRING "SourceForge.net account id for uploads")
MARK_AS_ADVANCED(SCP_USERNAME)

SET(NIGHTLY_TARGET "greedy-nightly-${CPACK_SYSTEM_NAME}.${CPACK_EXTENSION}")

SET(SCP_ROOT "github.com/pyushkevich/greedy")

#--------------------------------------------------------------------------------
# Create targets

ADD_CUSTOM_TARGET(upload_nightly 
  VERBATIM COMMAND "${SCP_PROGRAM}" ${SCP_ARGUMENTS}
  ${CPACK_PACKAGE_FILE_NAME_WEXT} ${SCP_USERNAME},c3d@${SCP_ROOT}/Nightly/${NIGHTLY_TARGET}
  DEPENDS ${CPACK_TARGET}
  WORKING_DIRECTORY ${SNAP_BINARY_DIR}
  COMMENT "Uploading package ${CPACK_PACKAGE_FILE_NAME_WEXT} to SourceForge.net as ${NIGHTLY_TARGET}")

ADD_CUSTOM_TARGET(upload_experimental 
	VERBATIM COMMAND "${SCP_PROGRAM}" ${SCP_ARGUMENTS} ${CPACK_PACKAGE_FILE_NAME_WEXT} ${SCP_USERNAME},c3d@${SCP_ROOT}/Experimental
  DEPENDS ${CPACK_TARGET}
  WORKING_DIRECTORY ${SNAP_BINARY_DIR}
  COMMENT "Uploading package ${CPACK_PACKAGE_FILE_NAME_WEXT} to SourceForge.net to Experimental directory")


