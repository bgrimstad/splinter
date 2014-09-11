# Enable C++11
if(CMAKE_COMPILER_IS_GNUCXX)
  add_definitions(-std=c++0x)
  add_definitions(-std=c++11)
endif()
