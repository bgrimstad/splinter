
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
    BSpline bspline1 = BSpline::Builder(samples).degree(BSpline::Degree::LINEAR).build();
    BSpline bspline3 = BSpline::Builder(samples).degree(BSpline::Degree::CUBIC).build();

    // Build penalized B-spline (P-spline) that smooths the samples
    BSpline pspline = BSpline::Builder(samples)
            .degree(BSpline::Degree::CUBIC)
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

##Restrictions
###B-splines
For the current implementation of B-splines we require that the samples you provide form a complete [grid](https://en.wikipedia.org/wiki/Regular_grid). This means that if the function you are sampling is two-dimensional with variables `x0` and `x1`, then all combinations of `x0` and `x1` must be present in the samples. This means that if you choose to sample `x1` in a new value, say 1, then you must sample `[x0 1]` for all previous values of `x0` used so far. In 2D you can visualize this as [graphing paper](https://en.wikipedia.org/wiki/Graph_paper#/media/File:Log_paper.svg), where all lines intersect. If a sample were missing, one of the intersections would be missing, and the grid would be incomplete. You can check if the grid is complete by calling `isGridComplete()` on your DataTable. This restriction will be removed in a later implementation.


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


###Polynomial regression
For polynomial regression it is important that the number of unique samples is larger or equal to the number of variables of the function to approximate. If not, you get an [underdetermined system](https://en.wikipedia.org/wiki/Underdetermined_system) with bad approximations.
