# This is a common CMake include file for creating version variables across
# different projects. It just defines some common Macros for versioning
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/CMake/")
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/CMake/rpavlik/")

# Get today's date (see http://cmake.3232098.n2.nabble.com/How-to-get-the-current-date-td5776870.html)
MACRO (TODAY RESULT)
    IF (WIN32)
        EXECUTE_PROCESS(COMMAND "cmd" " /C date /T" OUTPUT_VARIABLE ${RESULT})
        string(REGEX REPLACE "(..)/(..)/(....).*" "\\1/\\2/\\3" ${RESULT} ${${RESULT}})
    ELSEIF(UNIX)
        EXECUTE_PROCESS(COMMAND "date" "+%b %d, %Y" OUTPUT_VARIABLE ${RESULT})
        string(REGEX REPLACE "(...) (..), (....).*" "\\1 \\2, \\3" ${RESULT} ${${RESULT}})
    ELSE (WIN32)
        MESSAGE(SEND_ERROR "date not implemented")
        SET(${RESULT} 000000)
    ENDIF (WIN32)
	string(REPLACE "\n" "" ${RESULT} ${${RESULT}})
  # string(REPLACE " " "" ${RESULT} ${${RESULT}})
ENDMACRO (TODAY)

# Set the semantic version variables in a consistent manner for my projects
MACRO(VERSION_VARS MAJOR MINOR PATCH QUALIFIER DATE FMT_DATE)

  # Set the tool's semantic version
  SET(${PROJECT_NAME}_VERSION_MAJOR ${MAJOR})
  SET(${PROJECT_NAME}_VERSION_MINOR ${MINOR})
  SET(${PROJECT_NAME}_VERSION_PATCH ${PATCH})
  SET(${PROJECT_NAME}_VERSION_QUALIFIER "${QUALIFIER}")

  # These fields should also be modified each time
  SET(${PROJECT_NAME}_VERSION_RELEASE_DATE ${DATE})
  SET(${PROJECT_NAME}_VERSION_RELEASE_DATE_FORMATTED ${FMT_DATE})

  # Get today's date (see http://cmake.3232098.n2.nabble.com/How-to-get-the-current-date-td5776870.html)
  TODAY(${PROJECT_NAME}_VERSION_COMPILE_DATE)

  # This should not need to change
  SET(${PROJECT_NAME}_VERSION_FULL 
    "${${PROJECT_NAME}_VERSION_MAJOR}.${${PROJECT_NAME}_VERSION_MINOR}.${${PROJECT_NAME}_VERSION_PATCH}${${PROJECT_NAME}_VERSION_QUALIFIER}")

  # Get the current git hash
  include(GetGitRevisionDescription)
  get_git_head_revision(GIT_REFSPEC ${PROJECT_NAME}_VERSION_GIT_SHA1)

  # Get the current git branch
  include(GitBranch)
  get_git_branch(${PROJECT_NAME}_VERSION_GIT_BRANCH)
  get_git_commit_date(${${PROJECT_NAME}_VERSION_GIT_SHA1} ${PROJECT_NAME}_VERSION_GIT_TIMESTAMP)

  # Print the Git information
  MESSAGE(STATUS "${PROJECT_NAME} Version: ${${PROJECT_NAME}_VERSION_FULL} Released ${${PROJECT_NAME}_VERSION_RELEASE_DATE_FORMATTED}")
  MESSAGE(STATUS "GIT Info:")
  MESSAGE(STATUS "  Branch : ${${PROJECT_NAME}_VERSION_GIT_BRANCH}")
  MESSAGE(STATUS "  SHA    : ${${PROJECT_NAME}_VERSION_GIT_SHA1}")
  MESSAGE(STATUS "  Date   : ${${PROJECT_NAME}_VERSION_GIT_TIMESTAMP}")

ENDMACRO(VERSION_VARS)

