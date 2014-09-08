##Installation

You can compile the library yourself by following the guide below, or you can download binaries and headers [here](https://github.com/bgrimstad/multivariate-splines/releases).

There are two separate guides for compilation:
* [UNIX](#compile-on-unix)
* [Windows](#compile-on-windows)

Note that for both guides the same [options](#options) to CMake apply.

##Compile on UNIX
####Requirements
* [CMake](http://www.cmake.org/)
* [Git](http://git-scm.com/)
* [GCC](https://gcc.gnu.org/) or an equivalent C++11 compiler


1. Download and install [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page)
2. `sudo apt-get update && sudo apt-get install git cmake build-essential`
3. `git clone https://github.com/bgrimstad/multivariate-splines.git`
4. `cd multivariate-splines`
5. `mkdir build && cd build`
6. `cmake ..` Also, see [options](#options)
7. `make`
8. `make install`

####Troubleshooting
`fatal error: Eigen/Dense: No such file or directory`: The compiler could not find Eigen. You either need to install Eigen, and then run step #6 again (with `-DEIGEN_DIRECTORY="/path/to/eigen"` if Eigen did not install to the default directory (/usr/local/include/eigen3)).

`make: *** [install] Error 1`: You probably need elevated rights to install the library because you are trying to write to a directory you don't have permission to write to. Either change the install paths via the options, or run step #8 again like this: `sudo make install`.

---

##Compile on Windows

1. Clone https://github.com/bgrimstad/multivariate-bsplines
2. Download Eigen: http://eigen.tuxfamily.org/index.php?title=Main_Page
  1. Extract the zip-file into a new folder, and write down the location of that folder
3. Download and install CMake: http://www.cmake.org/
4. Download and install Qt Creator: http://qt-project.org/downloads
  1. Make sure that MinGW is marked for installation
5. Run Qt Creator, select `Open project`
  1. Navigate into the multivariate-bsplines folder, select `CMakeLists.txt`
  2. In the arguments field, write: `-DEIGEN_DIRECTORY=C:/path/to/eigen/from/step/2.1`. Also, see [options](#options)
  3. Run CMake
6. Now you can build the library with Qt Creator, and the library files will be output to your build directory.

* You may have to add -static as a flag to your linker if you are compiling with MinGW.
* C++11 must be enabled.
* If you get asked to specify where CMake is located, then you can typically find the file cmake.exe in C:/Program Files (x86)/CMake/bin.

####Troubleshooting
`fatal error: Eigen/Dense: No such file or directory`: The compiler could not find Eigen. You either need to install Eigen, and then run step #6 again (with `-DEIGEN_DIRECTORY="/path/to/eigen"` if Eigen did not install to the default directory (/usr/local/include/eigen3)).

---

##Options (both platforms)

| Variable name     | Default value                 | Description                               |
| ----------------- | ----------------------------- | ----------------------------------------- |
| EIGEN_DIRECTORY   | /usr/local/include/eigen3     | Path to the Eigen lib.                    |
| HEADER_DIRECTORY  | include                       | Where the headers should be installed.    |
| LIBRARY_DIRECTORY | lib                           | Where to install the library file.        |

These options go along with UNIX step #6 or Windows step #5.2, and are used like this:

    -DEIGEN_DIRECTORY=/home/me/eigen
    
    -DEIGEN_DIRECTORY=/path/to/eigen -DHEADER_DIRECTORY=/home/me/c++/multivariate-bsplines/includes

The syntax is: `-D<VARIABLE_NAME>=<VARIABLE_VALUE>`. If you have any spaces in your value you must surround it with double quotes (").

Note for the header and library paths:
If the path is relative (the first character is not / on UNIX or C:/ (or equivalent) on Windows), then the actual path used will be relative to [CMAKE_INSTALL_PREFIX](http://www.cmake.org/cmake/help/v2.8.12/cmake.html#variable:CMAKE_INSTALL_PREFIX).
