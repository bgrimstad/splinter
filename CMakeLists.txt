# This file is part of the SPLINTER library.
# Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

project(SPLINTER CXX)
cmake_minimum_required(VERSION 2.8)
string(TOLOWER ${PROJECT_NAME} PROJECT_NAME_LOWER)

# Read the SPLINTER version from version file
file(STRINGS "version" VERSION)

# GCC 4.* required for some applications
# NOTE for Python users: if the library is compiled with GCC > 4.8 the libgcc library must be installed (conda install libgcc)
#set(CMAKE_CXX_COMPILER /usr/bin/g++-4.8)

# Default configuration values
set(DEFAULT_BUILD_TYPE "Release")
set(DEFAULT_ARCH "x86-64")
set(DEFAULT_HEADER_INSTALL_DIRECTORY "include")
set(DEFAULT_LIBRARY_INSTALL_DIRECTORY "lib")
set(DEFAULT_EIGEN_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/Eigen)

# Detect platform
string(COMPARE EQUAL "Linux"   ${CMAKE_SYSTEM_NAME} LINUX)
string(COMPARE EQUAL "Windows" ${CMAKE_SYSTEM_NAME} WINDOWS)
string(COMPARE EQUAL "Darwin"  ${CMAKE_SYSTEM_NAME} OS_X)
if(LINUX)
    set(OS_STRING "linux")
elseif(WINDOWS)
    set(OS_STRING "windows")
elseif(OS_X)
    set(OS_STRING "osx")
else()
    set(OS_STRING "Unknown")
endif()

# Detect compiler
string(COMPARE EQUAL "Clang"      ${CMAKE_CXX_COMPILER_ID} CLANG)
string(COMPARE EQUAL "AppleClang" ${CMAKE_CXX_COMPILER_ID} APPLECLANG)
string(COMPARE EQUAL "GNU"        ${CMAKE_CXX_COMPILER_ID} GCC)
string(COMPARE EQUAL "Intel"      ${CMAKE_CXX_COMPILER_ID} INTEL)
string(COMPARE EQUAL "MSVC"       ${CMAKE_CXX_COMPILER_ID} MSVC)

if (CLANG OR APPLECLANG)
    set(CLANG 1)
endif()

# Set build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE ${DEFAULT_BUILD_TYPE})
endif()
string(TOUPPER ${CMAKE_BUILD_TYPE} CMAKE_BUILD_TYPE)

# Set architecture to default if not set
if(NOT ARCH)
    set(ARCH ${DEFAULT_ARCH})
endif()

# Set bitness
if(${ARCH} STREQUAL "x86-64")
    set(BITNESS "64")
    set(BITNESS_FLAG "-m${BITNESS}")
elseif(${ARCH} STREQUAL "x86")
    set(BITNESS "32")
    set(BITNESS_FLAG "-m${BITNESS}")
elseif(${ARCH} STREQUAL "armv8")
    set(BITNESS "64")
    set(BITNESS_FLAG "")
else()
    set(BITNESS "32")
    set(ARCH "x86")
    set(BITNESS_FLAG "-m${BITNESS}")
    message("Unknown architecture selected, defaulting to x86")
endif()

if(WINDOWS AND GCC)
    # MinGW (GCC on Windows) version 4.8.1 has a bug where _hypot is not found
    # when using any level of optimization. Defining __NO_INLINE__ fixes the bug.
    # The bug probably exists on other versions of MinGW as well, but that has not been verified by us.
    if(CMAKE_CXX_COMPILER_VERSION VERSION_EQUAL "4.8.1")
        set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -D__NO_INLINE__")
    endif()
    if(CMAKE_CXX_COMPILER_VERSION VERSION_EQUAL "4.9.2")
        set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O2")
        # From the GCC man pages:
        # If you use multiple -O options, with or without level numbers, the last such option is the one that is effective.
        set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -finline-functions")
        set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -funswitch-loops")
        set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -fpredictive-commoning")
        set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -fgcse-after-reload")
        set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -ftree-loop-vectorize")
        set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -ftree-loop-distribute-patterns")
        set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -ftree-slp-vectorize")
        set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -fvect-cost-model")
        set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -ftree-partial-pre")
        # -fipa-cp-clone causes a SegFault in Eigen when enabled with MinGW i686 4.9.2 with dwarf exception model
        #set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -fipa-cp-clone")
    endif()

    # To avoid an external dependency (libgcc_s_seh11.dll or libgcc_s_dw2-1.dll)
    # we statically link the libraries required by MinGW runtimes
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -static")
endif()

if(MSVC)
    # User cannot specify bitness with MSVC, so set it to whatever the generator is.
    string(TOLOWER ${CMAKE_GENERATOR} GENERATOR)
    if(GENERATOR MATCHES ".*win64.*")
        set(BITNESS "64")
    else()
        set(BITNESS "32")
    endif()
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Ox")

elseif(GCC OR CLANG)
    message("Identified compiler as GCC or Clang")
    # Treat warning return-type as error to avoid undefined behaviour
    # when a non-void function does not return a value.
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${BITNESS_FLAG} -std=c++11 -Werror=return-type")
    message("CMAKE_CXX_FLAGS: ${CMAKE_CXX_FLAGS}")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Wno-long-long")

elseif(INTEL)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${BITNESS_FLAG} -std=c++11")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall")
endif()

if(NOT EIGEN_DIRECTORY)
    set(EIGEN_DIRECTORY ${DEFAULT_EIGEN_DIRECTORY})
endif()
# CACHE STRING = Display the option with help text in CMakeCache.txt
set(EIGEN_DIRECTORY ${EIGEN_DIRECTORY} CACHE STRING "Directory where the Eigen library is located.")

if(NOT HEADER_INSTALL_DIRECTORY)
    set(HEADER_INSTALL_DIRECTORY ${DEFAULT_HEADER_INSTALL_DIRECTORY})
endif()
# CACHE STRING = Display the option with help text in CMakeCache.txt
set(HEADER_INSTALL_DIRECTORY ${HEADER_INSTALL_DIRECTORY} CACHE STRING "Absolute path, or, if relative, relative to CMAKE_INSTALL_PREFIX to install the header files.")

if(NOT LIBRARY_INSTALL_DIRECTORY)
    set(LIBRARY_INSTALL_DIRECTORY ${DEFAULT_LIBRARY_INSTALL_DIRECTORY})
endif()
# CACHE STRING = Display the option with help text in CMakeCache.txt
set(LIBRARY_INSTALL_DIRECTORY ${LIBRARY_INSTALL_DIRECTORY} CACHE STRING "Absolute path, or, if relative, relative to CMAKE_INSTALL_PREFIX to install the library file.")

# Header directories
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
include_directories(${EIGEN_DIRECTORY})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/test)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/thirdparty)

set(CINTERFACE_SRC_LIST
    include/cinterface/cinterface.h
    include/cinterface/utilities.h
    src/cinterface/bspline.cpp
    src/cinterface/bspline_builders.cpp
    src/cinterface/cinterface.cpp
    src/cinterface/data_table.cpp
    src/cinterface/utilities.cpp
)
# These are the sources we need for compilation of the library
set(SRC_LIST
    ${CINTERFACE_SRC_LIST}
    include/bspline.h
    include/bspline_basis.h
    include/bspline_basis_1d.h
    include/knot_builders.h
    include/knot_vector.h
    include/bspline_builders.h
    include/bspline_utils.h
    include/data_point.h
    include/data_table.h
    include/function.h
    include/definitions.h
    include/linear_solvers.h
    include/kronecker_product.h
    include/utilities.h
    include/json_parser.h
    src/bspline.cpp
    src/bspline_basis.cpp
    src/bspline_basis_1d.cpp
    src/knot_builders.cpp
    src/knot_vector.cpp
    src/bspline_builders.cpp
    src/bspline_utils.cpp
    src/data_point.cpp
    src/data_table.cpp
    src/function.cpp
    src/kronecker_product.cpp
    src/utilities.cpp
    src/json_parser.cpp)

set(TEST_SRC_LIST
    ${SRC_LIST}
    test/main.cpp
    test/approximation/bspline.cpp
    test/approximation/bspline_wls.cpp
    test/approximation/pspline.cpp
    test/general/bspline.cpp
    test/general/bspline_builder.cpp
    test/unit/data_table.cpp
    test/general/utilities.cpp
    test/examples/approximation.cpp
    test/utils/op_overloads.h
    test/utils/op_overloads.cpp
    test/utils/test_function.h
    test/utils/test_function.cpp
    test/utils/test_function_collection.h
    test/utils/test_function_collection.cpp
    test/utils/test_function_utils.h
    test/utils/test_function_utils.cpp
    test/utils/test_utils.h
    test/utils/test_utils.cpp
    test/utils/bspline_test_utils.h
    test/utils/bspline_test_utils.cpp
    test/unit/bspline_basis_1d.cpp
    test/unit/bspline_hessian.cpp
    test/unit/bspline_utils.cpp
    test/unit/knots.cpp
    test/unit/knot_vector.cpp
    test/unit/json_parser.cpp
    test/utils/timer.h
    test/utils/bspline_collection.cpp
    test/utils/bspline_collection.h
    test/utils/data_table_collection.cpp
    test/utils/data_table_collection.h
    test/examples/fit_random_1d_spline.cpp)

set(SHARED_LIBRARY ${PROJECT_NAME_LOWER}-${VERSION})
set(STATIC_LIBRARY ${PROJECT_NAME_LOWER}-static-${VERSION})

set(SHARED_LIBRARY_TEST "${PROJECT_NAME_LOWER}-shared-test")
set(STATIC_LIBRARY_TEST "${PROJECT_NAME_LOWER}-static-test")
set(TEST "${PROJECT_NAME_LOWER}-test")

# Add output library: add_library(libname [SHARED | STATIC] sourcelist)
add_library(${SHARED_LIBRARY} SHARED ${SRC_LIST})
add_library(${STATIC_LIBRARY} STATIC ${SRC_LIST})

# Create test executable
add_executable(${TEST} ${TEST_SRC_LIST})
target_link_libraries(${TEST} ${STATIC_LIBRARY})

# Copy test resources
file(COPY test/resources/ DESTINATION ${CMAKE_BINARY_DIR}/test-resources)

# Make the directory structure expected by the Python interface
# Python files (.py)
install(
        DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python/
        DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/splinter-python
        PATTERN "*.pyc" EXCLUDE
        PATTERN "__pycache__*" EXCLUDE
)
install(
        FILES ${CMAKE_CURRENT_SOURCE_DIR}/LICENSE
        DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/splinter-python/
)
# The Python interface reads the library version from this file
file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/splinter-python/splinterpy/version ${VERSION})
install(
        TARGETS ${SHARED_LIBRARY}
        DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/splinter-python/splinterpy/lib/${OS_STRING}/${ARCH}
)

# Install the header files (including Eigen) to the header directory
install(
    DIRECTORY
    ${CMAKE_CURRENT_SOURCE_DIR}/include/
    ${EIGEN_DIRECTORY}/Eigen
    ${EIGEN_DIRECTORY}/unsupported
    DESTINATION ${HEADER_INSTALL_DIRECTORY}/SPLINTER
)
# Install the shared library file
install(
    TARGETS ${SHARED_LIBRARY}
    DESTINATION ${LIBRARY_INSTALL_DIRECTORY}
)
# Install the static library file
install(
    TARGETS ${STATIC_LIBRARY}
    DESTINATION ${LIBRARY_INSTALL_DIRECTORY}
)

# Doxygen documentation
find_package(Doxygen)

if(DOXYGEN_FOUND)
  configure_file(
    ${PROJECT_SOURCE_DIR}/docs/Doxyfile.in
    ${PROJECT_BINARY_DIR}/Doxyfile
    @ONLY)

  add_custom_target(
    doc
#      ALL
    ${DOXYGEN_EXECUTABLE}
    ${PROJECT_BINARY_DIR}/Doxyfile
    COMMENT "Generating API documentation with Doxygen"
    VERBATIM)

endif()

if(CMAKE_BUILD_TYPE STREQUAL DEBUG)
    message("Detected compiler and platform:")
    message("Clang:   ${CLANG}")
    message("GCC:     ${GCC}")
    message("Intel:   ${INTEL}")
    message("MSVC:    ${MSVC}")
    message("Linux:   ${LINUX}")
    message("Windows: ${WINDOWS}")
    message("OS X:    ${OS_X}")
endif()

message("Configuring ${PROJECT_NAME} version ${VERSION} in ${CMAKE_BUILD_TYPE} mode for ${ARCH} (${BITNESS} bit)")
message("Compiler flags: ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE}}")
