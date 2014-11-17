##Installation

You can compile the library yourself by following the guide below, or you can download binaries and headers [here](https://github.com/bgrimstad/multivariate-splines/releases).

There are two separate guides for compilation:
* [UNIX](#compile-on-unix)
* [Windows](#compile-on-windows)

Note that for both guides the same [options](#options-both-platforms) to CMake apply.

CMake will try to automatically find where Eigen is installed, but if this fails you can specify it yourself with the option `EIGEN_DIRECTORY`. See [options](#options-both-platforms) for further information on how to use it.

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
6. `cmake ..` Also, see [options](#options-both-platforms)
7. `make`
8. `make install`

####Troubleshooting
`fatal error: Eigen/Dense: No such file or directory`: The compiler could not find Eigen. You either need to install Eigen, and then run step #6 again (with `-DEIGEN_DIRECTORY="/path/to/eigen"` if Eigen did not install to the default directory (/usr/local/include/eigen3)).

`make: *** [install] Error 1`: You probably need elevated rights to install the library because you are trying to write to a directory you don't have permission to write to. Either change the install paths via the options, or run step #8 again like this: `sudo make install`.

---

##Compile on Windows

1. Clone https://github.com/bgrimstad/multivariate-bsplines
2. Download [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page).
  1. Extract the zip-file into a new folder, and write down the location of that folder
3. Download and install [CMake](http://www.cmake.org/download/).
4. Download and install [MinGW](http://sourceforge.net/projects/mingw/files/latest/download?source=files).
5. Run MinGW (it should start automatically after you install it).
  1. Mark mingw32-base, mingw32-gcc-g++ and msys-base for installation.
  2. Select Installation -> Apply Changes.
6. Add "C:\MinGW\bin;C:\MinGW\msys\1.0\bin" to your PATH variable (Google is your friend).
7. Run CMake
  1. In "Where is the source code", select the multivariate-splines folder.
  2. In "Where to build the binaries", make a new folder somewhere and select that one (from now on called your "build" folder).
  3. Click Configure.
  4. Select MinGW Makefiles and Use default native compilers and click Finish.
  5. If you get errors, just click OK then repeat the past two steps and it should work.
  6. If you get errors saying CMake could not find Eigen, then you need to specify the path from step 2.1 manually:
    1. Click Add Entry.
    2. In "Name", write "EIGEN_DIRECTORY"
    3. Select Type: String
    4. In "Value", write the path from step 2.1.
    5. Click OK.
    6. Re-run Configure.
  7. If you wish to add more arguments you can do this the same way as in step 6.
  3. Click Generate.
8. Open a command window (Start -> type cmd -> hit enter)
  1. Navigate into your build folder
  2. Type "mingw32-make" and wait for the library to compile.

The output files should be in your build folder. If you wish to have them installed you can type "mingw32-make install".

####Troubleshooting
`fatal error: Eigen/Dense: No such file or directory`: The compiler could not find Eigen. You either need to install Eigen, and then run step #6 again (with `-DEIGEN_DIRECTORY="/path/to/eigen"` if Eigen did not install to the default directory (C:/Program Files/eigen3)).

`cc1plus.exe:-1: error: out of memory allocating (...)`: There is a bug with MinGW that causes this error. The newest version we've found that correctly compiles the library can be found [here](http://sourceforge.net/projects/mingwbuilds/files/host-windows/releases/4.7.0/64-bit/threads-posix/sjlj/). If that doesn't work, you can try [this one](https://code.google.com/p/mingw-builds/downloads/detail?name=x86_64-mingw32-gcc-4.7.0-release-c%2Cc%2B%2B%2Cfortran-sjlj.zip&can=2&q=). To install those alternative MinGW versions you can just delete everything in C:\mingw and replace it with what you find in the archives you downloaded.
Note: Somehow the 32-bit version of 4.7.0 will compile the library, but crashes at runtime. If you know a workaround for this issue, please contact us or send us a pull request with the workaround!

---

##Options (both platforms)

| Variable name     | Default value                    | Description                               |
| ----------------- | -------------------------------- | ----------------------------------------- |
| EIGEN_DIRECTORY   | UNIX: /usr/local/include/eigen3  | Path to the Eigen lib.                    |
|                   | Windows: C:/Program Files/eigen3 | Path to the Eigen lib.                    |
| HEADER_DIRECTORY  | include                          | Where the headers should be installed.    |
| LIBRARY_DIRECTORY | lib                              | Where to install the library file.        |

These options go along with UNIX step #6 or Windows step #5.2, and are used like this:

    UNIX:
    cmake .. -DEIGEN_DIRECTORY=/path/to/eigen -DHEADER_DIRECTORY=/home/me/multivariate-bsplines/includes
    
    
    Windows (in the arguments field):
    -DEIGEN_DIRECTORY=/path/to/eigen -DHEADER_DIRECTORY=/home/me/multivariate-bsplines/includes

The syntax is: `-D<VARIABLE_NAME>=<VARIABLE_VALUE>`. If you have any spaces in your value you must surround it with double quotes (").

Note for the header and library paths:
If the path is relative (the first character is not / on UNIX or C:/ (or equivalent) on Windows), then the actual path used will be relative to [CMAKE_INSTALL_PREFIX](http://www.cmake.org/cmake/help/v2.8.12/cmake.html#variable:CMAKE_INSTALL_PREFIX).
