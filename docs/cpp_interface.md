## C++ interface
The C++ interface is the native interface to SPLINTER. At any given time, the C++ interface will be the most comprehensive among the interfaces, exposing the most of SPLINTER's features.

Below is a simple example demonstrating the use of SPLINTER. Remember to compile with a C++11 compatible compiler!

More examples can be found in the directory `/test/examples`.

```c++
#include <iostream>
#include "bspline_builders.h"

using std::cout;
using std::endl;

using namespace SPLINTER;

int main(int argc, char *argv[])
{
    // We want to approximate this function, known as the six-hump camelback function
    auto f = [](std::vector<double> x) {
        assert(x.size() == 2);
        auto x0 = x.at(0);
        auto x1 = x.at(1);
        return (4 - 2.1*x0*x0 + (1/3.)*x0*x0*x0*x0)*x0*x0 + x0*x1 + (-4 + 4*x1*x1)*x1*x1;
    };

    // Create new DataTable to manage samples
    DataTable samples;

    // Sample the function
    for (auto i = 0; i < 20; i++)
    {
        for(auto j = 0; j < 20; j++)
        {
            // Sample function at x
            std::vector<double> x = {i*0.1, j*0.1};
            auto y = f(x);

            // Store sample
            samples.add_sample(x, y);
        }
    }

    // Build B-splines that interpolate the samples
    BSpline bspline1 = bspline_interpolator(samples, 1);
    BSpline bspline3 = bspline_interpolator(samples, 3);

    // Build degree 3 penalized B-spline (P-spline) that smooths the samples
    // The smoothing/regularization parameter is set to 0.03
    BSpline pspline = bspline_smoother(samples, 3, BSpline::Smoothing::PSPLINE, 0.03);

    /*
     * Evaluate the splines at x = (1,1)
     * NOTE1: the error will be 0 at that point (except for the P-spline, which may introduce an error
     * in favor of a smooth approximation) because it is a point we sampled at.
     * NOTE2: The BSpline::eval function returns an output vector (in this case of size 1)
     */
    std::vector<double> x = {1, 1};
    auto y_f = f(x);
    auto y_bs1 = bspline1.eval(x).at(0);
    auto y_bs3 = bspline3.eval(x).at(0);
    auto y_ps = pspline.eval(x).at(0);

    // Print results
    cout << "-----------------------------------------------------" << endl;
    cout << "Function at x:                 " << y_f                << endl;
    cout << "Linear B-spline at x:          " << y_bs1              << endl;
    cout << "Cubic B-spline at x:           " << y_bs3              << endl;
    cout << "P-spline at x:                 " << y_ps               << endl;
    cout << "-----------------------------------------------------" << endl;

    return 0;
}
```

### Sampling with DataTable
To simplify sampling in C++, SPLINTER comes with a DataTable data structure for managing and storing sample points. The following code snippet shows how DataTable can be used to manage samples. 
```c++
// Create new data structure
DataTable samples; 

// Add some samples (x, y), where y = f(x).
// Note that the order in which the samples are added does not matter.
samples.add_sample(1, 0);
samples.add_sample(2, 5);
samples.add_sample(3, 10);
samples.add_sample(4, 15);
```
