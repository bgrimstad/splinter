# This file is part of the Multivariate Splines library.
# Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com)
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

project(multivariate-splines)
cmake_minimum_required(VERSION 2.8)
set(CMAKE_BUILD_TYPE Release)

set(MULTIVARIATESPLINES_VERSION_MAJOR 0)
set(MULTIVARIATESPLINES_VERSION_MINOR 1)
set(MULTIVARIATESPLINES_VERSION ${MULTIVARIATESPLINES_VERSION_MAJOR}.${MULTIVARIATESPLINES_VERSION_MINOR})

set(DEFAULT_EIGEN_DIRECTORY "/usr/local/include/eigen3")
set(DEFAULT_HEADER_DIRECTORY "include")
set(DEFAULT_LIBRARY_DIRECTORY "lib")

if(NOT EIGEN_DIRECTORY)
    set(EIGEN_DIRECTORY ${DEFAULT_EIGEN_DIRECTORY})
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
    include/basis.h
    include/basis1d.h
    include/bspline.h
    include/datasample.h
    include/datatable.h
    include/generaldefinitions.h
    include/linearsolvers.h
    include/mykroneckerproduct.h
    include/pspline.h
    include/rbfspline.h
    include/spline.h
    src/basis.cpp
    src/basis1d.cpp
    src/bspline.cpp
    src/datasample.cpp
    src/datatable.cpp
    src/mykroneckerproduct.cpp
    src/pspline.cpp
    src/rbfspline.cpp
)

# Add output library: add_library(libname [SHARED | STATIC] sourcelist)
add_library(multivariate-splines SHARED ${SRC_LIST})
add_library(multivariate-splines-static STATIC ${SRC_LIST})

# Testing executable
add_executable("multivariate-splines-test" ${CMAKE_CURRENT_SOURCE_DIR}/test/main.cpp)
target_link_libraries("multivariate-splines-test" multivariate-splines)

# Compiler flags
add_definitions(
    -std=c++11
)

# Install the header files
install(
    DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/
    DESTINATION ${HEADER_DIRECTORY}
)

# Install the shared library file
install(
    TARGETS multivariate-splines
    DESTINATION ${LIBRARY_DIRECTORY}
)
# Install the static library file
install(
    TARGETS multivariate-splines-static
    DESTINATION ${LIBRARY_DIRECTORY}
)
