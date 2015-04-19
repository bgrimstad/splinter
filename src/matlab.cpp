#include "matlab.h"
#include <bspline.h>
#include <datatable.h>

using namespace Splinter;

	double f(double x)
	{
		return (4 - 2.1*x*x
			+ (1 / 3.)*x*x*x*x)*x*x
			+ x;
	}

#ifdef __cplusplus
	extern "C"
	{
#endif

		int *init() {
			DataTable samples;
			double y;
			double x;
			for (int i = 0; i < 20; i++)
			{
				// Sample function at x
				x = i*0.1;
				y = f(x);

				// Store sample
				samples.addSample(x, y);
			}

			return (int *) new BSpline(samples, BSplineType::CUBIC_FREE);
		}
		double eval(int *bspline, double x) {

			return ((BSpline *)bspline)->eval(x);
		}
#ifdef __cplusplus
	}
#endif
