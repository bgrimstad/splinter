#ifndef SPLINTER_MATLAB_H
#define SPLINTER_MATLAB_H

/*#define obj_ptr int * */
typedef int *obj_ptr;

#ifndef API
# ifdef _MSC_VER
#  define API __declspec(dllexport)
# else
#  define API
# endif
#endif

#ifdef __cplusplus
	extern "C"
	{
#endif
		/* 1 if the previous function call resulted in an error, 0 otherwise. */
		API int get_error();


		API obj_ptr datatable_init();

		API void datatable_add_sample(obj_ptr datatable_ptr, double *x, int x_dim, double y);
		
		API void datatable_add_samples(obj_ptr datatable_ptr, double *x, int n_samples, int x_dim);

		API void datatable_delete(obj_ptr datatable_ptr);


		API obj_ptr bspline_init(obj_ptr datatable_ptr, int type);

		API double bspline_eval(obj_ptr bspline_ptr, double *x, int n);

		API double *bspline_eval_jacobian(obj_ptr bspline_ptr, double *x, int n);

		API double *bspline_eval_hessian(obj_ptr bspline_ptr, double *x, int n);

		API void bspline_delete(obj_ptr bspline_ptr);


		API obj_ptr pspline_init(obj_ptr pspline_ptr, double lambda);

		API double pspline_eval(obj_ptr pspline_ptr, double *x, int x_dim);

		API double *pspline_eval_jacobian(obj_ptr pspline_ptr, double *x, int x_dim);

		API double *pspline_eval_hessian(obj_ptr pspline_ptr, double *x, int x_dim);

		API void pspline_delete(obj_ptr pspline_ptr);


		API obj_ptr rbf_init(obj_ptr rbf_ptr, int type_index, int normalized);

		API double rbf_eval(obj_ptr rbf_ptr, double *x, int x_dim);

		API double *rbf_eval_jacobian(obj_ptr rbf_ptr, double *x, int x_dim);

		API double *rbf_eval_hessian(obj_ptr rbf_ptr, double *x, int x_dim);

		API void rbf_delete(obj_ptr rbf_ptr);
#ifdef __cplusplus
	}
#endif


#endif // SPLINTER_MATLAB_H