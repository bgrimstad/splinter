
##Basic usage

The workflow to construct an approximation is simple: sample a function and construct an approximation. As the following figure illustrates, this process can be run iteratively until a satisfactory approximation has been built. To assess the accuracy of the approximation one can use existing samples for cross-validation or perform additional sampling. Note that the current version of SPLINTER only facilitates sampling and model construction. 

![Possible workflow with SPLINTER.](../assets/workflow.png)
Figure: A possible workflow for building approximations with SPLINTER.

The header files and classes intended for the end user of this library are:
[DataTable](../include/datatable.h), [BSpline](../include/bspline.h), [BSplineBuilder](../include/bsplinebuilder.h).

This is a simple example demonstrating the use of SPLINTER.

Remember to compile with a C++11 compatible compiler! That means you probably have to add a flag when compiling.

```c++
#include <iostream>
#include "datatable.h"

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
    BSpline bspline1 = BSpline::Builder(samples).degree(1).build();
    BSpline bspline3 = BSpline::Builder(samples).degree(3).build();

    // Build penalized B-spline (P-spline) that smooths the samples
    BSpline pspline = BSpline::Builder(samples)
            .degree(3)
            .smoothing(BSpline::Smoothing::PSPLINE)
            .alpha(0.03)
            .build();

    /* Evaluate the approximants at x = (1,1)
     * Note that the error will be 0 at that point (except for the P-spline, which may introduce an error
     * in favor of a smooth approximation) because it is a point we sampled at.
     */
    x(0) = 1; x(1) = 1;
    cout << "-----------------------------------------------------" << endl;
    cout << "Function at x:                 " << f(x)               << endl;
    cout << "Linear B-spline at x:          " << bspline1.eval(x)   << endl;
    cout << "Cubic B-spline at x:           " << bspline3.eval(x)   << endl;
    cout << "P-spline at x:                 " << pspline.eval(x)    << endl;
    cout << "-----------------------------------------------------" << endl;

    return 0;
}
```
