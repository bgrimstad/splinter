#include "matlab.h"
#include <bspline.h>
#include <datatable.h>

using namespace Splinter;

	/*double f(double x)
	{
		return (4 - 2.1*x*x
			+ (1 / 3.)*x*x*x*x)*x*x
			+ x;
	}*/

// TODO: #ifdef probably isn't needed here because we will always compile this file
// with a C++11 compiler. Make sure it is correct, though.
#ifdef __cplusplus
	extern "C"
	{
#endif
		/* DataTable interface */
		/* Cast the obj_ptr to a DataTable * */
		static DataTable *get_datatable(obj_ptr datatable_ptr) {
			return (DataTable *)datatable_ptr;
		}

		/* Constructor */
		obj_ptr datatable_init() {
			return (obj_ptr) new DataTable();
		}

		void datatable_add_sample(obj_ptr datatable_ptr, double x, double y) {
			DataTable *dataTable = get_datatable(datatable_ptr);

			dataTable->addSample(x, y);
		}

		void datatable_delete(obj_ptr datatable_ptr) {
			DataTable *dataTable = get_datatable(datatable_ptr);

			delete dataTable;
		}

		/* BSpline interface */
		/* Cast the obj_ptr to a BSpline * */
		static BSpline *get_bspline(obj_ptr datatable_ptr) {
			return (BSpline *)datatable_ptr;
		}

		// TODO: Make sure MatLab provides the correct type for the BSplineType enum
		/* Constructor */
		obj_ptr bspline_init(obj_ptr datatable_ptr, int type) {
			DataTable *table = get_datatable(datatable_ptr);

			BSplineType bsplineType;
			switch (type) {
			case 1: {
				bsplineType = BSplineType::LINEAR;
				break;
			}
			case 3: {
				bsplineType = BSplineType::CUBIC_FREE;
				break;
			}
			case 4: {
				bsplineType = BSplineType::QUADRATIC_FREE;
				break;
			}
			default: {
				bsplineType = BSplineType::CUBIC_FREE;
				break;
			}
			}

			return (obj_ptr) new BSpline(*table, bsplineType);
		}

		double bspline_eval(obj_ptr bspline_ptr, double x) {
			BSpline *bspline = get_bspline(bspline_ptr);

			return bspline->eval(x);
		}

		void bspline_delete(obj_ptr bspline_ptr) {
			BSpline *bspline = get_bspline(bspline_ptr);

			delete bspline;
		}

		/*int *init() {
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
		}*/
#ifdef __cplusplus
	}
#endif
