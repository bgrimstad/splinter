#ifndef SPLINTER_MATLAB_H
#define SPLINTER_MATLAB_H

#define obj_ptr int *

#ifndef API
# define API __declspec(dllexport)
#endif

#ifdef __cplusplus
	extern "C"
	{
#endif
		/* DataTable interface */
		/* Constructor */
		API obj_ptr datatable_init();

		API void datatable_add_sample(obj_ptr datatable_ptr, double x, double y);

		API void datatable_delete(obj_ptr datatable_ptr);

		/* BSpline interface */
		// TODO: Make sure MatLab provides the correct type for the BSplineType enum
		/* Constructor */
		API obj_ptr bspline_init(obj_ptr datatable_ptr, int type);

		API double bspline_eval(obj_ptr bspline_ptr, double x);

		API void bspline_delete(obj_ptr bspline_ptr);

		/*int API *init();
		double API eval(int *spline, double x);*/
#ifdef __cplusplus
	}
#endif


#endif // SPLINTER_MATLAB_H