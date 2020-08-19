PROJECT(GREEDY)

# Required CMake version
CMAKE_MINIMUM_REQUIRED(VERSION 3.3)

#--------------------------------------------------------------------------------
# VERSION
#--------------------------------------------------------------------------------
# INCLUDE(CMake/PYVersion.cmake)

# Set the semantic version and version date
# VERSION_VARS(1 0 1 "" "20190321" "Mar 21, 2019")

#--------------------------------------------------------------------------------
# Build configuration
#--------------------------------------------------------------------------------
INCLUDE(CMakeDependentOption)

# Subproject build options
OPTION(GREEDY_BUILD_AS_SUBPROJECT "Build greedy inside of another CMake build" OFF)
MARK_AS_ADVANCED(GREEDY_BUILD_AS_SUBPROJECT)
IF(GREEDY_BUILD_AS_SUBPROJECT)
  OPTION(GREEDY_SUBPROJECT_BUILD_CLI "Build the executables" OFF)
  OPTION(GREEDY_SUBPROJECT_INSTALL_CLI "Install the executables" OFF)
  SET(BUILD_CLI ${GREEDY_SUBPROJECT_BUILD_CLI})
  SET(INSTALL_CLI ${GREEDY_SUBPROJECT_INSTALL_CLI})
ELSE()
  SET(BUILD_CLI ON)
  SET(INSTALL_CLI ON)
ENDIF()

# Do we want to build LDDMM?
OPTION(GREEDY_BUILD_LDDMM "Build experimental LDDMM implementation" OFF)

# Do we want to use FFTW?
CMAKE_DEPENDENT_OPTION(GREEDY_USE_FFTW 
  "Use features provided by the FFTW library, only for experimental LDDMM code" OFF
  "GREEDY_BUILD_LDDMM" OFF)

# Do we want to build MACF tools?
OPTION(GREEDY_BUILD_MACF "Build experimental multi-atlas correspondence fusion code" OFF)

# Do we want to build the stack registration binary
OPTION(GREEDY_BUILD_STACK_REG "Build experimental slice stack registration code" OFF)

#--------------------------------------------------------------------------------
# Dependent packages
#--------------------------------------------------------------------------------

# ITK
FIND_PACKAGE(ITK 4.12.2 REQUIRED)
INCLUDE(${ITK_USE_FILE})

# Deal with FFTW - only used by experimental LDDMM code
IF(GREEDY_USE_FFTW)
  SET(ITK_USE_FFTWF ON)
  SET(ITK_USE_FFTWD ON)
  INCLUDE(CMake/FindFFTW.cmake)
  ADD_DEFINITIONS(-D_LDDMM_FFT_)
ENDIF()

# Include the header directories
INCLUDE_DIRECTORIES(
  ${GREEDY_SOURCE_DIR}/src
  ${GREEDY_SOURCE_DIR}/src/ITKFilters/include
)

#--------------------------------------------------------------------------------
# Hard-coded documentation
#--------------------------------------------------------------------------------

# Markdown documentation compiled into the C code
# modified from: https://github.com/starseeker/tinyscheme-cmake/blob/master/CMakeLists.txt
# # Rather than load the init.scm file at run time,
# # with the uncertainties as to where exactly the file
# # resides, use William Ahern's hexdump to generate
# # an embeddable version. Build our own copy of hexdump
# # to ensure consistent behavior and portability.
# # See http://25thandclement.com/~william/projects/hexdump.c.html
#ADD_EXECUTABLE(text_to_hex utilities/hexdump.c)
#SET_PROPERTY(TARGET text_to_hex APPEND PROPERTY COMPILE_DEFINITIONS "HEXDUMP_MAIN")

# List of docs that we want to generate
#SET(HARDCODED_DOCS
#  stackg_usage_main stackg_usage_init stackg_usage_recon
#  stackg_usage_volmatch stackg_usage_voliter stackg_usage_splat)
#
#UNSET(HARDCODED_DOCS_HEADERS)
#FOREACH(doc ${HARDCODED_DOCS})
#
#  ADD_CUSTOM_COMMAND(
#    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${doc}.h
#    COMMAND text_to_hex -i ${CMAKE_CURRENT_SOURCE_DIR}/docs/usage/${doc}.txt > #${CMAKE_CURRENT_BINARY_DIR}/${doc}.h
#    DEPENDS text_to_hex ${CMAKE_CURRENT_SOURCE_DIR}/docs/usage/${doc}.txt)
#
#  SET(HARDCODED_DOCS_HEADERS ${HARDCODED_DOCS_HEADERS} ${CMAKE_CURRENT_BINARY_DIR}/${doc}.h)
#
#ENDFOREACH()

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

#ADD_CUSTOM_TARGET(docs_to_hex ALL DEPENDS ${HARDCODED_DOCS_HEADERS})

#--------------------------------------------------------------------------------
# Source file specification
#--------------------------------------------------------------------------------

# One of the files needs to be configured (to insert version info)
CONFIGURE_FILE(
  ${CMAKE_CURRENT_SOURCE_DIR}/src/GreedyVersion.cxx.in
  ${CMAKE_CURRENT_BINARY_DIR}/GreedyVersion.cxx @ONLY IMMEDIATE)

# Define header files
SET(HEADERS
  src/ITKFilters/include/FastLinearInterpolator.h
  src/ITKFilters/include/FastWarpCompositeImageFilter.h
  src/ITKFilters/include/FastWarpCompositeImageFilter.txx
  src/ITKFilters/include/JacobianDeterminantImageFilter.h
  src/ITKFilters/include/JacobianDeterminantImageFilter.txx
  src/ITKFilters/include/MultiComponentImageMetricBase.h
  src/ITKFilters/include/MultiComponentImageMetricBase.txx
  src/ITKFilters/include/MultiComponentMetricReport.h
  src/ITKFilters/include/MultiComponentMutualInfoImageMetric.h
  src/ITKFilters/include/MultiComponentMutualInfoImageMetric.txx
  src/ITKFilters/include/MultiComponentNCCImageMetric.h
  src/ITKFilters/include/MultiComponentNCCImageMetric.txx
  src/ITKFilters/include/MultiImageAffineMSDMetricFilter.h
  src/ITKFilters/include/MultiImageAffineMSDMetricFilter.txx
  src/ITKFilters/include/MultiImageOpticalFlowImageFilter.h
  src/ITKFilters/include/MultiImageOpticalFlowImageFilter.txx
  src/ITKFilters/include/OneDimensionalInPlaceAccumulateFilter.h
  src/ITKFilters/include/OneDimensionalInPlaceAccumulateFilter.txx
  src/ITKFilters/include/SimpleWarpImageFilter.h
  src/ITKFilters/include/SimpleWarpImageFilter.txx
  src/ITKFilters/include/itkGaussianInterpolateImageFunction.h
  src/ITKFilters/include/itkOptVectorLinearInterpolateImageFunction.h
  src/ITKFilters/include/itkOptVectorLinearInterpolateImageFunction.txx
  src/lddmm_common.h
  src/lddmm_data.h
  src/AffineCostFunctions.h
  src/AffineTransformUtilities.h
  src/GreedyAPI.h
  src/GreedyException.h
  src/GreedyParameters.h
  src/MultiImageRegistrationHelper.h
  src/CommandLineHelper.h
)

# Define greedy library files
SET(GREEDY_LIB_SRC
  src/lddmm_data.cxx
  src/GreedyAPI.cxx
  src/GreedyParameters.cxx
  src/MultiImageRegistrationHelper.cxx
  src/AffineCostFunctions.cxx
)

SET(LDDMM_SRC src/lddmm_main.cxx)
SET(GREEDY_SRC 
  ${CMAKE_CURRENT_BINARY_DIR}/GreedyVersion.cxx
  src/greedy_main.cxx)

SET(MACF_SRC src/macf_optimization.cxx)

SET(STACK_GREEDY_SRC
  ${CMAKE_CURRENT_BINARY_DIR}/GreedyVersion.cxx
  src/dijkstra/ShortestPath.cxx
  src/greedy_stack_reg.cxx)

ADD_LIBRARY(greedyapi ${GREEDY_LIB_SRC} ${HEADERS})
TARGET_INCLUDE_DIRECTORIES(greedyapi PUBLIC ${GREEDY_INCLUDE_DIRS})

# The executables are only compiled when the software is built as its own project
IF(BUILD_CLI)
  IF(GREEDY_BUILD_LDDMM)
    ADD_EXECUTABLE(lddmm ${LDDMM_SRC})
    TARGET_LINK_LIBRARIES(lddmm greedyapi
      ${ITK_LIBRARIES} ${FFTWF_LIB}
      ${FFTWF_THREADS_LIB})
  ENDIF()

  IF(GREEDY_BUILD_MACF)
    ADD_EXECUTABLE(macf ${MACF_SRC})
    TARGET_LINK_LIBRARIES(macf greedyapi ${ITK_LIBRARIES})
  ENDIF()

  IF(GREEDY_BUILD_STACK_REG)
    ADD_EXECUTABLE(stack_greedy ${STACK_GREEDY_SRC})
    INCLUDE_DIRECTORIES(${GREEDY_SOURCE_DIR}/src/dijkstra)
    TARGET_LINK_LIBRARIES(stack_greedy greedyapi ${ITK_LIBRARIES} ${FFTWF_LIB} ${FFTWF_THREADS_LIB})
    ADD_DEPENDENCIES(stack_greedy docs_to_hex)
  ENDIF()

  ADD_EXECUTABLE(greedy ${GREEDY_SRC})
  TARGET_LINK_LIBRARIES(greedy greedyapi ${ITK_LIBRARIES} ${FFTWF_LIB} ${FFTWF_THREADS_LIB})
  #ADD_EXECUTABLE(test_accum testing/src/TestOneDimensionalInPlaceAccumulateFilter.cxx)
  #TARGET_LINK_LIBRARIES(test_accum ${ITK_LIBRARIES})
ENDIF(BUILD_CLI)

# Install command-line executables
IF(INSTALL_CLI)
  IF(GREEDY_SUBPROJECT_CLI_INSTALL_PATH)
    SET(CLI_INSTALL_PATH ${GREEDY_SUBPROJECT_CLI_INSTALL_PATH})
  ELSE()
    SET(CLI_INSTALL_PATH bin)
  ENDIF()

  INSTALL(TARGETS greedy DESTINATION ${CLI_INSTALL_PATH} COMPONENT Runtime)
ENDIF(INSTALL_CLI)

IF(NOT GREEDY_BUILD_AS_SUBPROJECT)

  # ENABLE_TESTING()
  # INCLUDE(CTest)

ENDIF(NOT GREEDY_BUILD_AS_SUBPROJECT)
