
##Basic usage

The workflow to construct an approximation is simple: sample a function and construct an approximation. As the following figure illustrates, this process can be run iteratively until a satisfactory approximation has been built. To assess the accuracy of the approximation one can use existing samples for cross-validation or perform additional sampling. Note that the current version of SPLINTER only facilitates sampling and model construction. 

![Possbile workflow with SPLINTER.](../assets/workflow.png)
Figure: A possible workflow for building approximations with SPLINTER.

The header files and classes intended for the end user of this library are:
[DataTable](../include/datatable.h), [BSpline](../include/bspline.h), [BSplineType](../include/bspline.h), [PSpline](../include/pspline.h), [RadialBasisFunction](../include/radialbasisfunction.h), [RadialBasisFunctionType](../include/radialbasisfunctionterm.h) and [PolynomialRegression](../include/polynomialregression.h).

This is a simple example demonstrating the use of Splinter. Note that there is no restrictions to the dimension of x (except that it has to be >= 1, of course).

Remember to compile with a c++11 compatible compiler! That means you probably have to add a flag when compiling.

```c++
#include <iostream>
#include "datatable.h"
#include "bspline.h"
#include "pspline.h"
#include "radialbasisfunction.h"
#include "polynomialregression.h"

using std::cout;
using std::endl;

using namespace SPLINTER;

// Six-hump camelback function
double f(DenseVector x)
{
    assert(x.rows() == 2);
    return (4 - 2.1*x(0)*x(0)
           + (1/3.)*x(0)*x(0)*x(0)*x(0))*x(0)*x(0)
           + x(0)*x(1)
           + (-4 + 4*x(1)*x(1))*x(1)*x(1);
}

int main(int argc, char *argv[])
{
    // Create new DataTable to manage samples
    DataTable samples;

    // Sample the function
    DenseVector x(2);
    double y;
    for(int i = 0; i < 20; i++)
    {
        for(int j = 0; j < 20; j++)
        {
            // Sample function at x
            x(0) = i*0.1;
            x(1) = j*0.1;
            y = f(x);

            // Store sample
            samples.addSample(x,y);
        }
    }

    // Build B-splines that interpolate the samples
    BSpline bspline1(samples, BSplineType::LINEAR);
    BSpline bspline3(samples, BSplineType::CUBIC);

    // Build penalized B-spline (P-spline) that smooths the samples
    PSpline pspline(samples, 0.03);

    // Build radial basis function spline that interpolate the samples
    RadialBasisFunction rbfspline(samples, RadialBasisFunctionType::THIN_PLATE_SPLINE);

    /* The six-hump camelback is a function of degree 6 in x(0) and degree 4 in x(1),
     * therefore a polynomial of degree 6 in the first variable and 4 in the second
     * will be sufficient to interpolate the function
     */
    auto degrees = std::vector<unsigned int>(2);
    degrees.at(0) = 6;
    degrees.at(1) = 4;
    PolynomialRegression polyfit(samples, degrees);

    /* Evaluate the approximants at x = (1,1)
     * Note that the error will be 0 at that point (except for the P-spline, which may introduce an error
     * in favor of a smooth approximation) because it is a point we sampled at.
     */
    x(0) = 1; x(1) = 1;
    cout << "-------------------------------------------------"        << endl;
    cout << "Function at x: \t\t\t"            << f(x)                 << endl;
    cout << "Linear B-spline at x: \t\t"       << bspline1.eval(x)     << endl;
    cout << "Cubic B-spline at x: \t\t"        << bspline3.eval(x)     << endl;
    cout << "P-spline at x: \t\t\t"            << pspline.eval(x)      << endl;
    cout << "Thin-plate spline at x:\t\t"      << rbfspline.eval(x)    << endl;
    cout << "Polynomial of degree 6 at x:\t"   << polyfit.eval(x)      << endl;
    cout << "-------------------------------------------------"        << endl;

    return 0;
}
```

###Sampling with DataTable
Function samples are managed by and stored in the DataTable data structure. The following code snippet shows how DataTable is used to manage samples. 
```c++
// Create new data structure
DataTable samples; 

// Add some samples (x,y), where y = f(x)
samples.addSample(1,0);
samples.addSample(2,5);
samples.addSample(3,10);
samples.addSample(4,15);

// The order in which the samples are added does not matter
// since DataTable keeps the samples sorted internally.
```
The grid, meaning all the values of x where you have sampled the function, must be complete. This means that if you have sampled the function at `x = [0 0]`, `x = [1 0]` and `x = [2 1]`, you must also sample the function at `x = [1 1]`, `x = [0 1]` and `x = [2 0]`. You must have sampled the function in all permutations of x within the possible values of x<sub>0</sub>, x<sub>1</sub> ... x<sub>n</sub>. The number of samples will then (disregarding duplicates) be num(x<sub>0</sub>) * num(x<sub>1</sub>) * ... * num(x<sub>n</sub>), where num(x) is the number of distinct sample points x. You can check if the grid is complete by calling `isGridComplete()` on your DataTable.


This is an **incomplete** grid:

| x<sub>0</sub>   | x<sub>1</sub>   | y   |
| --------------- | --------------- | --- |
| 2.1             | 1               | -7  |
| 2.3             | 3               | 10  |
| 2.1             | 3               | 9.3 |


This is a **complete** grid:

| x<sub>0</sub>   | x<sub>1</sub>   | y   |
| --------------- | --------------- | --- |
| 2.1             | 1               | -7  |
| 2.3             | 3               | 10  |
| 2.1             | 3               | 9.3 |
| 2.3             | 1               | 0   |

Please note that whether the grid is complete or not only depends on the values of x, not those of y.
