#include "matlab.h"
#include <bspline.h>
#include <pspline.h>
#include <datatable.h>
#include <generaldefinitions.h>
#include <iostream>
#include <fstream>

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

		void datatable_add_sample(obj_ptr datatable_ptr, double *x, int x_dim, double y) {
			DataTable *dataTable = get_datatable(datatable_ptr);
			DenseVector vec(x_dim);
			for (int i = 0; i < x_dim; i++) {
				vec(i) = x[i];
			}

			dataTable->addSample(vec, y);
		}

		API void datatable_add_samples(obj_ptr datatable_ptr, double *x, int n_samples, int x_dim) {
			DataTable *dataTable = get_datatable(datatable_ptr);
			DenseVector vec(x_dim);

			int k = 0;
			for (int i = 0; i < n_samples; ++i) {
				for (int j = 0; j < x_dim; ++j) {
					vec(j) = x[i + j * n_samples];
				}
				dataTable->addSample(vec, x[i + x_dim * n_samples]);
			}
		}

		void datatable_delete(obj_ptr datatable_ptr) {
			DataTable *dataTable = get_datatable(datatable_ptr);

			delete dataTable;
		}

		/* BSpline interface */
		/* Cast the obj_ptr to a BSpline * */
		BSpline *get_bspline(obj_ptr bspline_ptr) {
			return (BSpline *) bspline_ptr;
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
			case 2: {
				bsplineType = BSplineType::QUADRATIC;
				break;
			}
			case 3: {
				bsplineType = BSplineType::CUBIC;
				break;
			}
			case 4: {
				bsplineType = BSplineType::QUARTIC;
				break;
			}
			default: {
				return 0;
			}
			}

			return (obj_ptr) new BSpline(*table, bsplineType);
		}

		double bspline_eval(obj_ptr bspline_ptr, double *x, int x_dim) {
			BSpline *bspline = get_bspline(bspline_ptr);
			DenseVector vec(x_dim);
			for (int i = 0; i < x_dim; i++) {
				vec(i) = x[i];
			}

			return bspline->eval(vec);
		}
		
		/*
			Return double ptr so MatLab identifies the libpointer type as
			'doublePtr'. Then we can do reshape(ptr, x_dim, y_dim) in MatLab
			to get a MatLab matrix.
		*/
		double *bspline_eval_jacobian(obj_ptr bspline_ptr, double *x, int x_dim) {
			BSpline *bspline = get_bspline(bspline_ptr);
			DenseVector vec(x_dim);
			for (int i = 0; i < x_dim; i++) {
				vec(i) = x[i];
			}

			DenseMatrix jacobian = bspline->evalJacobian(vec);

			int numCoefficients = jacobian.cols() * jacobian.rows();
			double *res = (double *) malloc(sizeof(double) * numCoefficients);
			memcpy(res, jacobian.data(), sizeof(double) * numCoefficients);
			return res;
		}

		double *bspline_eval_hessian(obj_ptr bspline_ptr, double *x, int x_dim) {
			BSpline *bspline = get_bspline(bspline_ptr);
			DenseVector vec(x_dim);
			for (int i = 0; i < x_dim; i++) {
				vec(i) = x[i];
			}

			DenseMatrix hessian = bspline->evalHessian(vec);

			// The DenseMatrix lives on the stack, so we need to copy
			// the data into a newly allocated memory area.
			int numCoefficients = hessian.cols() * hessian.rows();
			double *res = (double *) malloc(sizeof(double) * numCoefficients);
			memcpy(res, hessian.data(), sizeof(double) * numCoefficients);
			return res;
		}

		void bspline_delete(obj_ptr bspline_ptr) {
			BSpline *bspline = get_bspline(bspline_ptr);

			delete bspline;
		}

		/* PSpline interface */
		/* Cast the obj_ptr to a PSpline * */
		PSpline *get_pspline(obj_ptr pspline_ptr) {
			return (PSpline *) pspline_ptr;
		}

		/* Constructor */
		obj_ptr pspline_init(obj_ptr pspline_ptr, double lambda) {
			DataTable *table = get_datatable(pspline_ptr);

			return (obj_ptr) new PSpline(*table, lambda);
		}

		double pspline_eval(obj_ptr pspline_ptr, double *x, int x_dim) {
			PSpline *pspline = get_pspline(pspline_ptr);
			DenseVector vec(x_dim);
			for (int i = 0; i < x_dim; i++) {
				vec(i) = x[i];
			}

			return pspline->eval(vec);
		}

		/*
		Return double ptr so MatLab identifies the libpointer type as
		'doublePtr'. Then we can do reshape(ptr, x_dim, y_dim) in MatLab
		to get a MatLab matrix.
		*/
		double *pspline_eval_jacobian(obj_ptr pspline_ptr, double *x, int x_dim) {
			PSpline *pspline = get_pspline(pspline_ptr);
			DenseVector vec(x_dim);
			for (int i = 0; i < x_dim; i++) {
				vec(i) = x[i];
			}

			DenseMatrix jacobian = pspline->evalJacobian(vec);

			int numCoefficients = jacobian.cols() * jacobian.rows();
			double *res = (double *)malloc(sizeof(double) * numCoefficients);
			memcpy(res, jacobian.data(), sizeof(double) * numCoefficients);
			return res;
		}

		double *pspline_eval_hessian(obj_ptr pspline_ptr, double *x, int x_dim) {
			PSpline *pspline = get_pspline(pspline_ptr);
			DenseVector vec(x_dim);
			for (int i = 0; i < x_dim; i++) {
				vec(i) = x[i];
			}

			DenseMatrix hessian = pspline->evalHessian(vec);

			// The DenseMatrix lives on the stack, so we need to copy
			// the data into a newly allocated memory area.
			int numCoefficients = hessian.cols() * hessian.rows();
			double *res = (double *)malloc(sizeof(double) * numCoefficients);
			memcpy(res, hessian.data(), sizeof(double) * numCoefficients);
			return res;
		}

		void pspline_delete(obj_ptr pspline_ptr) {
			PSpline *pspline = get_pspline(pspline_ptr);

			delete pspline;
		}
#ifdef __cplusplus
	}
#endif
