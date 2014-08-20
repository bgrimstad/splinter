project(InterpolationLibrary)
cmake_minimum_required(VERSION 2.8)
set(CMAKE_BUILD_TYPE Release)

set(BSPLINE_VERSION_MAJOR 0)
set(BSPLINE_VERSION_MINOR 1)
set(BSPLINE_VERSION ${BSPLINE_VERSION_MAJOR}.${BSPLINE_VERSION_MINOR})


set(DEFAULT_EIGEN_DIRECTORY "/usr/local/include/eigen3")
set(DEFAULT_HEADER_DIRECTORY "include/bspline")
set(DEFAULT_LIBRARY_DIRECTORY "lib/bspline")


if(NOT EIGEN_DIRECTORY)
    set(EIGEN_DIRECTORY ${DEFAULT_EIGEN_DIRECTORY})
endif(NOT EIGEN_DIRECTORY)
# Display the option with help text in CMakeCache.txt
set(EIGEN_DIRECTORY ${EIGEN_DIRECTORY} CACHE STRING "Directory where the Eigen library is located.")

if(NOT HEADER_DIRECTORY)
    set(HEADER_DIRECTORY ${DEFAULT_HEADER_DIRECTORY})
endif(NOT HEADER_DIRECTORY)
# Display the option with help text in CMakeCache.txt
set(HEADER_DIRECTORY ${HEADER_DIRECTORY} CACHE STRING "Absolute path, or, if relative, relative to CMAKE_INSTALL_PREFIX to install the header files.")

if(NOT LIBRARY_DIRECTORY)
    set(LIBRARY_DIRECTORY ${DEFAULT_LIBRARY_DIRECTORY})
endif(NOT LIBRARY_DIRECTORY)
# Display the option with help text in CMakeCache.txt
set(LIBRARY_DIRECTORY ${LIBRARY_DIRECTORY} CACHE STRING "Absolute path, or, if relative, relative to CMAKE_INSTALL_PREFIX to install the library file.")


# These are the headers we need for compilation
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
include_directories(${EIGEN_DIRECTORY})

# These are the sources we need for compilation
aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR}/src SRC_LIST)

# Remove main.cpp from library sources (it is only used for testing)
set(LIB_SRC_LIST ${SRC_LIST})
list(REMOVE_ITEM LIB_SRC_LIST ${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp)

# Add output library: add_library(libname [SHARED | STATIC] sourcelist)
add_library(BSpline SHARED ${LIB_SRC_LIST})
add_library(BSplineStatic STATIC ${LIB_SRC_LIST})

add_executable("bspline_test" ${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp)
target_link_libraries("bspline_test" BSpline)

# Compiler flags
add_definitions(
    -std=c++11
)

# Install the header files the user will see
install(
    FILES bspline.h sorteddatatable.h radialbasisfunction.h
    DESTINATION ${HEADER_DIRECTORY}
)

# Install the header files that are used internally
install(
    DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/
    DESTINATION ${HEADER_DIRECTORY}/include
)

# Install the shared library file
install(
    TARGETS BSpline
    DESTINATION ${LIBRARY_DIRECTORY}
)
# Install the static library file
install(
    TARGETS BSplineStatic
    DESTINATION ${LIBRARY_DIRECTORY}
)
