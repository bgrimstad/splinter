project(InterpolationLibrary)
cmake_minimum_required(VERSION 2.8)
aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR} SRC_LIST)
#aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR}/Eigen SRC_LIST)
aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR}/src SRC_LIST)

# Found here: http://stackoverflow.com/questions/17511496/create-a-shared-library-with-cmake
set(CMAKE_BUILD_TYPE Release)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
include_directories(/usr/local/include/eigen3)
add_library(BSplineInterpolation SHARED ${SRC_LIST})

add_definitions(
    -std=c++11
)

# Install the header files
install(
    DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/
    DESTINATION ${CMAKE_INSTALL_PREFIX}/include/bsplineinterpolation
)

# Install the .so file
install(
    TARGETS BSplineInterpolation
    DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/bsplineinterpolation
)
# Uncomment this line to generate executable
#add_executable(${PROJECT_NAME} ${SRC_LIST})
