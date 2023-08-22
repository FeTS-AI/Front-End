# This CMakeLists file will add the source files in this directory for compilation

# common_includes will have all the functions and classes from CBICA_Toolkit_ITK branch and more

# each header file, "$application.h" : main entry point of an application whose input will be either a single "itk::Image" or a "std::vector< itk::Image >"

# each application-specific include folder, "/$application_includes" : dependencies which are not taken care by the default dependencies (i.e., ITK and OpenCV), for example libSVM for recurrence

INCLUDE(CaPTk_macros)

FILE( GLOB
  COMMON_DEPENDENCIES
  "${CMAKE_CURRENT_SOURCE_DIR}/common_includes/*.*"

)
SET( SOURCES "" )
SET( DEPENDENCIES ${COMMON_DEPENDENCIES} )
SET( INCLUDES "")
SET( CAPTK_MACOSX_BUNDLE_FILES "")

SET( CAPTK_APP_LIST "" )

SET( CAPTK_BRAIN_APPS "" )
SET( CAPTK_BREAST_APPS "" )
SET( CAPTK_LUNG_APPS "" )
SET( CAPTK_MISC_APPS "" )

ADD_SUBDIRECTORY( Preprocessing )
ADD_SUBDIRECTORY( Utilities )

INCLUDE_DIRECTORIES( 
  ${CMAKE_CURRENT_SOURCE_DIR}/Preprocessing/src/depends
  ${CMAKE_CURRENT_SOURCE_DIR}/Utilities/src/depends
)

SET( INCLUDES ${INCLUDES} ${GEODESIC_TRAINING_INCLUDE_DIRS})

SET(PREPROCESSING_ALGOS "DCM2NIfTI Denoise-SUSAN Register BiasCorrect-N3 HistogramMatching SkullStripping" CACHE STRING "All Preprocessing algorithms available")

## this is for all stand-alone applications which do not integrate at source code level with CaPTk
SET( CAPTK_APP_LIST_PY_GUI "" )

## download and extract the prepared binaries of external apps

SET( DOWNLOADED_APPS_DIR ${PROJECT_BINARY_DIR}/externalApps )
SET( SOURCE_APPLICATIONS_PATH ${CMAKE_CURRENT_SOURCE_DIR}/individualApps )
SET( PLATFORM_STRING "" )

IF(WIN32)
  SET( PLATFORM_STRING "windows" )
ELSEIF(APPLE)
  SET( PLATFORM_STRING "macos" )
ELSE()
  SET( PLATFORM_STRING "linux" )
ENDIF()

SET( FILENAME_TO_EXTRACT "binaries_${PLATFORM_STRING}")
# SET( DOWNLOAD_LINK "ftp://www.nitrc.org/home/groups/captk/downloads/precompiledApps/${PLATFORM_STRING}.zip" )
SET( DOWNLOAD_LINK "https://github.com/CBICA/CaPTk/raw/master/binaries/precompiledApps/${PLATFORM_STRING}.zip" )
SET( LFS_FILE_TO_CHECK "${PROJECT_SOURCE_DIR}/binaries/precompiledApps/${PLATFORM_STRING}.zip" )

SET( FILE_TO_EXTRACT ${PROJECT_BINARY_DIR}/${FILENAME_TO_EXTRACT}.zip)

# # putting this condition so that the extraction doesn't happen every single time
# IF( NOT EXISTS "${DOWNLOADED_APPS_DIR}" )

#   IF( NOT EXISTS "${FILE_TO_EXTRACT}" )
  
#     # copy from LFS folder
#     IF( EXISTS ${LFS_FILE_TO_CHECK} )
#       CONFIGURE_FILE( ${LFS_FILE_TO_CHECK} ${FILE_TO_EXTRACT} )
#     ENDIF()

#     # do not re-download if the LFS fetch worked
#     IF(NOT EXISTS ${FILE_TO_EXTRACT})
#       # download exe from url
#       MESSAGE( STATUS "Downloading pre-compiled external applications" )
#       FILE(DOWNLOAD "${DOWNLOAD_LINK}" "${FILE_TO_EXTRACT}" TIMEOUT 1000000 STATUS STATUS_CODE SHOW_PROGRESS)
#       IF(NOT STATUS_CODE EQUAL 0)
#         MESSAGE(FATAL_ERROR "Failed to download Precompiled packages. Status=${STATUS_CODE}")
#       ENDIF()
#     ENDIF()
#   ENDIF()

#   # FILE(MAKE_DIRECTORY ${DOWNLOADED_APPS_DIR})

#   MESSAGE( STATUS "Extracting pre-compiled external applications" )
#   IF( EXISTS "${FILE_TO_EXTRACT}" )

#     EXECUTE_PROCESS( COMMAND ${CMAKE_COMMAND} -E tar xfz ${FILE_TO_EXTRACT}
#       WORKING_DIRECTORY ${DOWNLOADED_APPS_DIR}
#       RESULT_VARIABLE RESULT_CODE
#     )

#     IF(NOT RESULT_CODE EQUAL 0)
#       MESSAGE( WARNING "Extracting the pre-compiled applications failed" )
#     ENDIF()

#   ENDIF()
# ENDIF()

SUBDIRLIST(EXTERNAL_APPS ${DOWNLOADED_APPS_DIR}/${FILENAME_TO_EXTRACT}/)

# copy the prepared app into the appropriate folder
# FOREACH(external_app ${EXTERNAL_APPS})

#   IF( NOT EXISTS "${SOURCE_APPLICATIONS_PATH}/${external_app}" )

#     FILE(MAKE_DIRECTORY ${SOURCE_APPLICATIONS_PATH}/${external_app})
    
#     EXECUTE_PROCESS( COMMAND ${CMAKE_COMMAND} -E copy_directory
#       ${DOWNLOADED_APPS_DIR}/${FILENAME_TO_EXTRACT}/${external_app} ${SOURCE_APPLICATIONS_PATH}/${external_app}
#     )
  
#   ENDIF()

# ENDFOREACH()

SUBDIRLIST(SUBDIRECTORIES ${SOURCE_APPLICATIONS_PATH})
INCLUDE(${CMAKE_ROOT}/Modules/ExternalProject.cmake)

SET(BUILD_DEEPMEDIC false CACHE BOOL "disable deepmedic")

IF (BUILD_DEEPMEDIC)

  SET( CAPTK_APP_LIST_PY_GUI "${CAPTK_APP_LIST_PY_GUI} deepmedic" CACHE STRING "Available stand-alone apps" FORCE)
  ### this should be removed and put in with application decoupling
  ## this is a TEMPORARY work-around until https://github.com/CBICA/CaPTk/pull/1043 is merged
  SET( DOWNLOAD_LINK "ftp://www.nitrc.org/home/groups/captk/downloads/deepMedicInference_${PLATFORM_STRING}.zip" )
  SET( FILENAME_TO_EXTRACT "deepMedicInference_${PLATFORM_STRING}")
  SET( FILE_TO_EXTRACT ${PROJECT_BINARY_DIR}/${FILENAME_TO_EXTRACT}.zip)

  SET( DOWNLOADED_APPS_DIR ${PROJECT_BINARY_DIR}/deepMedicInference )

  IF( NOT EXISTS "${DOWNLOADED_APPS_DIR}" )

    IF( NOT EXISTS "${FILE_TO_EXTRACT}" )

      # do not re-download if the LFS fetch worked
      IF(NOT EXISTS ${FILE_TO_EXTRACT})
        # download exe from url
        MESSAGE( STATUS "Downloading pre-compiled deepMedicInference" )
        FILE(DOWNLOAD "${DOWNLOAD_LINK}" "${FILE_TO_EXTRACT}" TIMEOUT 1000000 STATUS STATUS_CODE SHOW_PROGRESS)
        IF(NOT STATUS_CODE EQUAL 0)
          MESSAGE(FATAL_ERROR "Failed to download Precompiled deepMedicInference. Status=${STATUS_CODE}")
        ENDIF()
      ENDIF()
    ENDIF()

    FILE(MAKE_DIRECTORY ${DOWNLOADED_APPS_DIR})

    MESSAGE( STATUS "Extracting pre-compiled deepMedicInference" )
    IF( EXISTS "${FILE_TO_EXTRACT}" )

      EXECUTE_PROCESS( COMMAND ${CMAKE_COMMAND} -E tar xfz ${FILE_TO_EXTRACT}
        WORKING_DIRECTORY ${DOWNLOADED_APPS_DIR}
        RESULT_VARIABLE RESULT_CODE
      )

      IF(NOT RESULT_CODE EQUAL 0)
        MESSAGE( WARNING "Extracting the pre-compiled deepMedicInference failed" )
      ENDIF()

    ENDIF()
  ENDIF()

  FILE(GLOB DMInferenceEXE "${DOWNLOADED_APPS_DIR}/*")

  IF(APPLE)
    INSTALL(FILES ${DMInferenceEXE}
      DESTINATION ${EXE_NAME}.app/Contents/Resources/bin
      PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
    )
  ELSE()
    INSTALL(FILES ${DMInferenceEXE}
      DESTINATION bin
      PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
    )
  ENDIF()
ENDIF()
### this should be removed and put in with application decoupling

FOREACH(subdir ${SUBDIRECTORIES})

  SET(SUBDIRPATH "${SOURCE_APPLICATIONS_PATH}/${subdir}")
  
  IF (BUILD_DEEPMEDIC)

    SET( CAPTK_APP_LIST_PY_GUI "${CAPTK_APP_LIST_PY_GUI} deepmedic" CACHE STRING "Available stand-alone apps" FORCE)
    ### this should be removed and put in with application decoupling
    ## this is a TEMPORARY work-around until https://github.com/CBICA/CaPTk/pull/1043 is merged
    SET( DOWNLOAD_LINK "ftp://www.nitrc.org/home/groups/captk/downloads/deepMedicInference_${PLATFORM_STRING}.zip" )
    SET( FILENAME_TO_EXTRACT "deepMedicInference_${PLATFORM_STRING}")
    SET( FILE_TO_EXTRACT ${PROJECT_BINARY_DIR}/${FILENAME_TO_EXTRACT}.zip)
    
    SET( DOWNLOADED_APPS_DIR ${PROJECT_BINARY_DIR}/deepMedicInference )
    
    IF( NOT EXISTS "${DOWNLOADED_APPS_DIR}" )
    
      IF( NOT EXISTS "${FILE_TO_EXTRACT}" )
    
        # do not re-download if the LFS fetch worked
        IF(NOT EXISTS ${FILE_TO_EXTRACT})
          # download exe from url
          MESSAGE( STATUS "Downloading pre-compiled deepMedicInference" )
          FILE(DOWNLOAD "${DOWNLOAD_LINK}" "${FILE_TO_EXTRACT}" TIMEOUT 1000000 STATUS STATUS_CODE SHOW_PROGRESS)
          IF(NOT STATUS_CODE EQUAL 0)
            MESSAGE(FATAL_ERROR "Failed to download Precompiled deepMedicInference. Status=${STATUS_CODE}")
          ENDIF()
        ENDIF()
      ENDIF()
    
      FILE(MAKE_DIRECTORY ${DOWNLOADED_APPS_DIR})
    
      MESSAGE( STATUS "Extracting pre-compiled deepMedicInference" )
      IF( EXISTS "${FILE_TO_EXTRACT}" )
    
        EXECUTE_PROCESS( COMMAND ${CMAKE_COMMAND} -E tar xfz ${FILE_TO_EXTRACT}
          WORKING_DIRECTORY ${DOWNLOADED_APPS_DIR}
          RESULT_VARIABLE RESULT_CODE
        )
    
        IF(NOT RESULT_CODE EQUAL 0)
          MESSAGE( WARNING "Extracting the pre-compiled deepMedicInference failed" )
        ENDIF()
    
      ENDIF()
    ENDIF()
    
    FILE(GLOB DMInferenceEXE "${DOWNLOADED_APPS_DIR}/*")
    
    IF(APPLE)
    INSTALL(FILES ${DMInferenceEXE}
      DESTINATION ${EXE_NAME}.app/Contents/Resources/bin
      PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
    )
    ELSE()
    INSTALL(FILES ${DMInferenceEXE}
      DESTINATION bin
      PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
    )
    ENDIF()
    ### this should be removed and put in with application decoupling

    # IF(WIN32)
    #   INSTALL(FILES ${SUBDIRPATH}/deepMedicRun.exe
    #     DESTINATION bin
    #     PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
    #   )
    # ELSEIF(APPLE)
    #   INSTALL(FILES ${SUBDIRPATH}/deepMedicRun
    #     DESTINATION ${EXE_NAME}.app/Contents/Resources/bin
    #     PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
    #   )
    # ELSE()
    #   INSTALL(FILES ${SUBDIRPATH}/deepMedicRun
    #     DESTINATION bin
    #     PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
    #   )
    # ENDIF()

  ENDIF()

ENDFOREACH()

#dcm2niix: Start 
#downloads to: [captk bin folder]
#extracts to: [captk bin folder]/dcm2niix
#copies to: [captk src folder]/src/applications/individualsApps/dcm2niix
#installs to: [bin folder]

#make dir
IF( NOT EXISTS "${SOURCE_APPLICATIONS_PATH}/dcm2niix" )
  FILE(MAKE_DIRECTORY ${SOURCE_APPLICATIONS_PATH}/dcm2niix)
ENDIF()

#specify download link
SET( DOWNLOAD_LINK "ftp://www.nitrc.org/home/groups/captk/downloads/precompiledApps/dcm2niix_${PLATFORM_STRING}.zip" )
SET( FILENAME_TO_EXTRACT "dcm2niix_${PLATFORM_STRING}")
SET( FILE_TO_EXTRACT ${PROJECT_BINARY_DIR}/${FILENAME_TO_EXTRACT}.zip)

SET( DOWNLOADED_APPS_DIR ${PROJECT_BINARY_DIR}/dcm2niix )

IF( NOT EXISTS "${DOWNLOADED_APPS_DIR}" )

  IF( NOT EXISTS "${FILE_TO_EXTRACT}" )

	# do not re-download if file exists
	IF(NOT EXISTS ${FILE_TO_EXTRACT})
#Download from url
	  MESSAGE( STATUS "Downloading pre-compiled dcm2niix" )
	  FILE(DOWNLOAD "${DOWNLOAD_LINK}" "${FILE_TO_EXTRACT}" TIMEOUT 1000000 STATUS STATUS_CODE SHOW_PROGRESS)
	  IF(NOT STATUS_CODE EQUAL 0)
		MESSAGE(FATAL_ERROR "Failed to download Precompiled dcm2niix. Status=${STATUS_CODE}")
	  ENDIF()
	ENDIF()
  ENDIF()

  FILE(MAKE_DIRECTORY ${DOWNLOADED_APPS_DIR})

#Extract
  MESSAGE( STATUS "Extracting pre-compiled dcm2niix" )
  IF( EXISTS "${FILE_TO_EXTRACT}" )
	EXECUTE_PROCESS( COMMAND ${CMAKE_COMMAND} -E tar xfz ${FILE_TO_EXTRACT}
	  WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/dcm2niix
	  RESULT_VARIABLE RESULT_CODE)

	IF(NOT RESULT_CODE EQUAL 0) #if extraction fails
	  MESSAGE( WARNING "Extracting the pre-compiled dcm2niix failed" )
	ELSE() # if extraction succeeds
	  MESSAGE( STATUS "Extraction succeeded. Now copying to ${SOURCE_APPLICATIONS_PATH}/dcm2niix")
	  FILE(GLOB allFiles ${PROJECT_BINARY_DIR}/dcm2niix/dcm2niix_${PLATFORM_STRING}/* )
#Copy
	  FILE( COPY ${allFiles} 
		DESTINATION ${SOURCE_APPLICATIONS_PATH}/dcm2niix )
	ENDIF()
  ENDIF()
ENDIF()

#Install
SET( DCM2NIIX_EXE "" )
IF (WIN32)
  FILE(GLOB allFiles ${SOURCE_APPLICATIONS_PATH}/dcm2niix/*.*)
  SET( DCM2NIIX_EXE ${allFiles} )
ELSE ()
  SET( DCM2NIIX_EXE "${SOURCE_APPLICATIONS_PATH}/dcm2niix/dcm2niix" )
ENDIF()

IF (APPLE)
  SET(CAPTK_MACOSX_BUNDLE_FILES ${CAPTK_MACOSX_BUNDLE_FILES} ${DCM2NIIX_EXE})
ELSE()
  INSTALL(FILES ${DCM2NIIX_EXE}
	DESTINATION bin
	PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
  )
ENDIF()
#dcm2niix: End

SET( CAPTK_APP_LIST ${CAPTK_APP_LIST} CACHE STRING "Available Code-level Applications" FORCE)
SET( APPLICATION_SOURCES  "${SOURCES}" CACHE STRING "Code-level Application sources" FORCE)
SET( APPLICATION_INCLUDES "${INCLUDES}" CACHE STRING "Code-level Application dependency folders" FORCE)
SET( APPLICATION_DEPENDS "${DEPENDENCIES}" CACHE STRING "Code-level Application dependencies" FORCE)
SET( CAPTK_MACOSX_BUNDLE_FILES_INSTALL_TARGETS "${CAPTK_MACOSX_BUNDLE_FILES}" CACHE STRING "macOS specific stuff" FORCE)
