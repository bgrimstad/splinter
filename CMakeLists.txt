# This file is part of the SPLINTER library.
# Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

project(SPLINTER)
cmake_minimum_required(VERSION 2.8)
string(TOLOWER ${PROJECT_NAME} PROJECT_NAME_LOWER)

set(VERSION_MAJOR 1)
set(VERSION_MINOR 3)
set(VERSION ${VERSION_MAJOR}-${VERSION_MINOR})

# Default configuration values
set(DEFAULT_BUILD_TYPE "Release")
set(DEFAULT_BITNESS "32") # Can be overriden by user or generator selection
set(DEFAULT_HEADER_INSTALL_DIRECTORY "include")
set(DEFAULT_LIBRARY_INSTALL_DIRECTORY "lib")
set(DEFAULT_EIGEN_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/Eigen)

# Detect platform
string(COMPARE EQUAL "Linux"   ${CMAKE_SYSTEM_NAME} LINUX)
string(COMPARE EQUAL "Windows" ${CMAKE_SYSTEM_NAME} WINDOWS)
string(COMPARE EQUAL "Darwin"  ${CMAKE_SYSTEM_NAME} OS_X)
if(LINUX)
    set(OS_STRING "Linux")
elseif(WINDOWS)
    set(OS_STRING "Windows")
elseif(OS_X)
    set(OS_STRING "OSX")
else()
    set(OS_STRING "Unknown")
endif()

# Detect compiler
string(COMPARE EQUAL "Clang" ${CMAKE_CXX_COMPILER_ID} CLANG)
string(COMPARE EQUAL "GNU"   ${CMAKE_CXX_COMPILER_ID} GCC)
string(COMPARE EQUAL "Intel" ${CMAKE_CXX_COMPILER_ID} INTEL)
string(COMPARE EQUAL "MSVC"  ${CMAKE_CXX_COMPILER_ID} MSVC)

# Set build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE ${DEFAULT_BUILD_TYPE})
endif()
string(TOUPPER ${CMAKE_BUILD_TYPE} CMAKE_BUILD_TYPE)

# Set bitness
if(BITNESS EQUAL "64")
	set(BITNESS "64")
else()
	set(BITNESS ${DEFAULT_BITNESS})
endif()

# GCC on Windows generates wrong code with -O3
if(WINDOWS AND GCC)
    set(CMAKE_CXX_FLAGS_RELEASE "-O2")
    # -ftree-vectorize is turned on by -O3 but not -O2.
    # Investigate if the issue of generating wrong code still exists before uncommenting the following line.
    # set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -ftree-vectorize")
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
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m${BITNESS} -std=c++11")
	set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Wno-long-long")
#   TODO: Intel flags
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

# These are the sources we need for compilation of the library
set(SRC_LIST
    src/bspline.cpp
    src/bsplinebasis.cpp
    src/bsplinebasis1d.cpp
    src/datasample.cpp
    src/datatable.cpp
    src/mykroneckerproduct.cpp
    src/pspline.cpp
    src/radialbasisfunction.cpp
)
set(MATLAB_SRC_LIST ${SRC_LIST} src/matlab.cpp)

set(SHARED_LIBRARY ${PROJECT_NAME_LOWER}-${VERSION})
set(STATIC_LIBRARY ${PROJECT_NAME_LOWER}-static-${VERSION})
set(MATLAB_LIBRARY ${PROJECT_NAME_LOWER}-matlab-${VERSION})

set(SHARED_LIBRARY_TEST "${PROJECT_NAME_LOWER}-shared-test")
set(STATIC_LIBRARY_TEST "${PROJECT_NAME_LOWER}-static-test")
set(SERIALIZE_TEST "${PROJECT_NAME_LOWER}-serialize-test")
set(TESTING_UTILITIES ${CMAKE_CURRENT_SOURCE_DIR}/test/testingutilities.cpp)

# Add output library: add_library(libname [SHARED | STATIC] sourcelist)
add_library(${SHARED_LIBRARY} SHARED ${SRC_LIST})
add_library(${STATIC_LIBRARY} STATIC ${SRC_LIST})
add_library(${MATLAB_LIBRARY} SHARED ${MATLAB_SRC_LIST})

# Testing executables
add_executable(${SHARED_LIBRARY_TEST} ${CMAKE_CURRENT_SOURCE_DIR}/test/main.cpp ${TESTING_UTILITIES})
target_link_libraries(${SHARED_LIBRARY_TEST} ${SHARED_LIBRARY})
add_executable(${STATIC_LIBRARY_TEST} ${CMAKE_CURRENT_SOURCE_DIR}/test/main.cpp ${TESTING_UTILITIES})
target_link_libraries(${STATIC_LIBRARY_TEST} ${STATIC_LIBRARY})

# Serialize testing executable
add_executable(${SERIALIZE_TEST} ${CMAKE_CURRENT_SOURCE_DIR}/test/serialize.cpp ${TESTING_UTILITIES})
target_link_libraries(${SERIALIZE_TEST} ${STATIC_LIBRARY})

# Make the directory structure expected by the MatLab interface
install(
    DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/matlab
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/splinter-matlab
)
install(
    FILES ${CMAKE_CURRENT_SOURCE_DIR}/include/matlab.h
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/splinter-matlab/include
)
install(
    TARGETS ${MATLAB_LIBRARY}
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/splinter-matlab/lib/${OS_STRING}/${BITNESS}bit
)

# Install the header files
install(
    DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/
    DESTINATION ${HEADER_INSTALL_DIRECTORY}
)
# Install Eigen to the header directory
install(
	DIRECTORY ${EIGEN_DIRECTORY}
	DESTINATION ${HEADER_INSTALL_DIRECTORY}
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

message("Configuring ${PROJECT_NAME} version ${VERSION} in ${CMAKE_BUILD_TYPE} mode (${BITNESS} bit)")
message("Compiler flags: ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE}}")