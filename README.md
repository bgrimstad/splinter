##BSpline
BSpline is a curve fitting library utilizing a speedy implementation of [BSplines](http://en.wikipedia.org/wiki/B-spline) with [Eigens](http://eigen.tuxfamily.org/index.php?title=Main_Page) SparseMatrix.

###Requirements: 
* [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page)
* [CMake](http://www.cmake.org/)
* [Git](http://git-scm.com/)
* [GCC](https://gcc.gnu.org/) or an equivalent C++11 compiler

If you are using Debian / Ubuntu or any derivatives of those then you can probably just open a terminal and write
```shell
sudo apt-get update && sudo apt-get install git cmake build-essential
```
to install Git, CMake and GCC.

Install instructions for Eigen can be found at their website.

###Steps to install (on UNIX, steps are probably similar on Windows)
1. `git clone`
2. `cd bspline`
3. `cmake .`
4. `make`
5. `make install`

####Options:
These options go along with step #3, and are used like this:

*     cmake . -DEIGEN_DIRECTORY=/home/me/eigen

*     cmake . -DEIGEN_DIRECTORY=/path/to/eigen -DHEADER_DIRECTORY=/home/me/c++/bspline/includes

The syntax is: `-D<VARIABLE_NAME>=<VARIABLE_VALUE>`. If you have any spaces in your value you must surround it with ".

| Variable name     | Default value             | Description                               |
| ----------------- | ------------------------- | ----------------------------------------- |
| EIGEN_DIRECTORY   | /usr/local/include/eigen3 | Path to the Eigen lib.                    |
| HEADER_DIRECTORY  | include/bspline           | Where the headers should be installed.    |
| LIBRARY_DIRECTORY | lib/bspline               | Where to install the library file.        |

Note for the header and library paths:
If the path is relative (the first character is not / on UNIX or C:/ (or equivalent) on Windows), then the actual path used will be relative to [CMAKE_INSTALL_PREFIX](http://www.cmake.org/cmake/help/v2.8.12/cmake.html#variable:CMAKE_INSTALL_PREFIX).

####Troubleshooting
`fatal error: Eigen/Dense: No such file or directory`: The compiler could not find Eigen. You either need to install Eigen, and then run step #3 again (with `-DEIGEN_DIRECTORY="/path/to/eigen"` if Eigen did not install to the default directory (/usr/local/include/eigen3)).

`make: *** [install] Error 1`: You probably need elevated rights to install the library because you are trying to write to a directory you don't have permission to write to. Either change the install paths via the options, or run step #5 again like this: `sudo make install`.

###Usage
This is a simple example demonstrating the use of BSpline. Note that there is no restrictions to the dimension of x or y (except that they have to be >= 1, of course), nor is there any requirement that their dimensions should be equal.

```c++
#include <iostream>
#include <sorteddatatable.h>
#include <bspline.h>
using std::cout;
using std::endl;

// Note: DenseVector is a typedef of Eigen::VectorXd (from generaldefinitions.h)
// typedef Eigen::VectorXd DenseVector;

// The function we are sampling.
// In reality this would probably be a set of sensors or other, more useful, data sources.
DenseVector func(DenseVector &x)
{
    DenseVector y(2);
    y(0) = x(0) / (1.3 * x(1) + 0.01) + 2 * x(1);
    y(1) = 10.2 * x(0) - 3 * x(1);
    return y;
}

int main()
{
    // The BSpline class needs it's data sorted by x0, x1 ... xn.
    // This class makes sure that it is. It also allows / disallows duplicates
    // and checks that the grid is complete.
    SortedDataTable table;

    DenseVector x(2);
    for(double i = 0.0; i <= 2.6; i += 0.11)
    {
        for(double j = 1.4; j <= 5.3; j += 0.05)
        {
            x(0) = i; x(1) = j;
            // addSample has several signatures, see the header file.
            table.addSample(x, func(x));
        }
    }

    Bspline bSpline(table, 3);

    // Make sure we evaluate the BSpline inside the domain of x
    DenseVector y, b;
    for(double i = 0.0; i <= 2.6; i += 0.07)
    {
        for(double j = 1.4; j <= 5.3; j += 0.09)
        {
            x(0) = i; x(1) = j;
            y = func(x);
            b = bSpline.evaluate(x);
            cout << y(0) << ", " << y(1) << " =? ";
            cout << b(0) << ", " << b(1) << endl;
        }
    }

    return 0;
}
```

###Grid
The "grid", meaning all the values of x where you have sampled the function, must be complete. This means that if you have sampled the function in `x = [0 0]`, `x = [1 0]` and `x = [2 1]`, you must also sample the function in `x = [1 1]`, `x = [0 1]` and `x = [2 0]`. You must have sampled the function in all permutations of x within the possible values of x0, x1 ... xn. The number of samples will then (disregarding duplicates) be num(x1) * num(x2) * ... * num(xn) where num(x) is the number of distinct values of x the the function has been sampled in. You can check if the grid is complete by calling `isGridComplete()` on your SortedDataTable.


This is **not** a complete grid:

| y0   | y1   | x0    | x1   |
| ---- | ---- | ----- | ---- |
| 1    | 2    | 2.1   | 1    |
| 0    | 20   | 2.3   | 3    |
| 19   | -1   | 2.1   | 3    |


This is a complete grid:

| y0   | y1   | x0    | x1   |
| ---- | ---- | ----- | ---- |
| 1    | 2    | 2.1   | 1    |
| 0    | 20   | 2.3   | 3    |
| 19   | -1   | 2.1   | 3    |
| -9   | 2.1  | 2.3   | 1    |

Please note that whether the grid is complete or not only depends on the values of x, not those of y.
