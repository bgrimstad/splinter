# This file is part of the Multivariate Splines library.
# Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com)
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.


project(multivariate-splines)
cmake_minimum_required(VERSION 2.8)
set(CMAKE_BUILD_TYPE Release)

set(MULTIVARIATESPLINES_VERSION_MAJOR 1)
set(MULTIVARIATESPLINES_VERSION_MINOR 1)
set(MULTIVARIATESPLINES_VERSION -${MULTIVARIATESPLINES_VERSION_MAJOR}-${MULTIVARIATESPLINES_VERSION_MINOR})

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")

include("${CMAKE_SOURCE_DIR}/cmake/c11Support.cmake")

if(UNIX)
    set(DEFAULT_EIGEN_DIRECTORY "/usr/local/include/eigen3")
elseif(WIN32)
    # TODO: This surely can't be the most logical place to look for Eigen
    set(DEFAULT_EIGEN_DIRECTORY "C:/Program Files/Eigen")
endif(UNIX)

set(DEFAULT_HEADER_DIRECTORY "include")
set(DEFAULT_LIBRARY_DIRECTORY "lib")

if(NOT EIGEN_DIRECTORY)
    find_package(Eigen3 REQUIRED) # Attempt to automatically locate Eigen
    if(EIGEN3_FOUND)
        set(EIGEN_DIRECTORY ${EIGEN3_INCLUDE_DIR})
    else(EIGEN3_FOUND)
        set(EIGEN_DIRECTORY ${DEFAULT_EIGEN_DIRECTORY})
    endif(EIGEN3_FOUND)
endif(NOT EIGEN_DIRECTORY)
# CACHE STRING = Display the option with help text in CMakeCache.txt
set(EIGEN_DIRECTORY ${EIGEN_DIRECTORY} CACHE STRING "Directory where the Eigen library is located.")

if(NOT HEADER_DIRECTORY)
    set(HEADER_DIRECTORY ${DEFAULT_HEADER_DIRECTORY})
endif(NOT HEADER_DIRECTORY)
# CACHE STRING = Display the option with help text in CMakeCache.txt
set(HEADER_DIRECTORY ${HEADER_DIRECTORY} CACHE STRING "Absolute path, or, if relative, relative to CMAKE_INSTALL_PREFIX to install the header files.")

if(NOT LIBRARY_DIRECTORY)
    set(LIBRARY_DIRECTORY ${DEFAULT_LIBRARY_DIRECTORY})
endif(NOT LIBRARY_DIRECTORY)
# CACHE STRING = Display the option with help text in CMakeCache.txt
set(LIBRARY_DIRECTORY ${LIBRARY_DIRECTORY} CACHE STRING "Absolute path, or, if relative, relative to CMAKE_INSTALL_PREFIX to install the library file.")

# These are the headers we need for compilation
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
include_directories(${EIGEN_DIRECTORY})

# These are the sources we need for compilation
#aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR}/src SRC_LIST)
set(SRC_LIST
    include/spline.h
    include/bspline.h
    include/bsplinebasis.h
    include/bsplinebasis1d.h
    include/pspline.h
    include/rbfspline.h
    include/datasample.h
    include/datatable.h
    include/generaldefinitions.h
    include/linearsolvers.h
    include/mykroneckerproduct.h
    src/bspline.cpp
    src/bsplinebasis.cpp
    src/bsplinebasis1d.cpp
    src/pspline.cpp
    src/rbfspline.cpp
    src/datasample.cpp
    src/datatable.cpp
    src/mykroneckerproduct.cpp
)

# Add output library: add_library(libname [SHARED | STATIC] sourcelist)
add_library(multivariate-splines${MULTIVARIATESPLINES_VERSION} SHARED ${SRC_LIST})
add_library(multivariate-splines-static${MULTIVARIATESPLINES_VERSION} STATIC ${SRC_LIST})

# Testing executable
add_executable("multivariate-splines-test" ${CMAKE_CURRENT_SOURCE_DIR}/test/main.cpp)
target_link_libraries("multivariate-splines-test" multivariate-splines-static${MULTIVARIATESPLINES_VERSION})

if(WIN32)
    set(CMAKE_EXE_LINKER_FLAGS "-static-libgcc -static-libstdc++") # Necessary for test executable to work with MinGW
endif()


# Install the header files
install(
    DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/
    DESTINATION ${HEADER_DIRECTORY}
)

# Install the shared library file
install(
    TARGETS multivariate-splines${MULTIVARIATESPLINES_VERSION}
    DESTINATION ${LIBRARY_DIRECTORY}
)
# Install the static library file
install(
    TARGETS multivariate-splines-static${MULTIVARIATESPLINES_VERSION}
    DESTINATION ${LIBRARY_DIRECTORY}
)

# Used to debug warnings
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-long-long -pedantic")
