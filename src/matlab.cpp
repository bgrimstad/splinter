#include "matlab.h"
#include <bspline.h>
#include <pspline.h>
#include <radialbasisfunction.h>
#include <datatable.h>
#include <generaldefinitions.h>
#include <iostream>
#include <fstream>
#include <set>

using namespace SPLINTER;

// Keep a list of objects so we avoid performing operations on objects that don't exist
std::set<obj_ptr> objects = std::set<obj_ptr>();

/* Cast the obj_ptr to a DataTable * */
static DataTable *get_datatable(obj_ptr datatable_ptr) {
	if (objects.count(datatable_ptr) > 0) {
		return (DataTable *)datatable_ptr;
	}

	return nullptr;
}

/* Cast the obj_ptr to a BSpline * */
static BSpline *get_bspline(obj_ptr bspline_ptr) {
	if (objects.count(bspline_ptr) > 0) {
		return (BSpline *)bspline_ptr;
	}

	return nullptr;
}

/* Cast the obj_ptr to a PSpline * */
static PSpline *get_pspline(obj_ptr pspline_ptr) {
	if (objects.count(pspline_ptr) > 0) {
		return (PSpline *)pspline_ptr;
	}

	return nullptr;
}

/* Cast the obj_ptr to a RadialBasisFunction * */
static RadialBasisFunction *get_rbf(obj_ptr rbf_ptr) {
	if (objects.count(rbf_ptr) > 0) {
		return (RadialBasisFunction *)rbf_ptr;
	}

	return nullptr;
}

extern "C"
{
	/* 1 if the last call to the library resulted in an error,
	*  0 otherwise. This is used to avoid throwing exceptions across library boundaries,
	*  and we expect the caller to manually check the value of this flag.
	*/
	int lastFuncCallError = 0;
	int get_error() {
		return lastFuncCallError;
	}

	/* DataTable interface */
	/* Constructor */
	obj_ptr datatable_init() {
		lastFuncCallError = 0;
		obj_ptr dataTable = (obj_ptr) new DataTable();
		if (dataTable == nullptr) {
			lastFuncCallError = 1;
			return nullptr;
		}

		objects.insert(dataTable);
		return dataTable;
	}

	void datatable_add_sample(obj_ptr datatable_ptr, double *x, int x_dim, double y) {
		lastFuncCallError = 0;
		DataTable *dataTable = get_datatable(datatable_ptr);
		if (dataTable == nullptr) {
			lastFuncCallError = 1;
			return;
		}

		DenseVector vec(x_dim);
		for (int i = 0; i < x_dim; i++) {
			vec(i) = x[i];
		}

		dataTable->addSample(vec, y);
	}

	void datatable_add_samples(obj_ptr datatable_ptr, double *x, int n_samples, int x_dim) {
		lastFuncCallError = 0;
		DataTable *dataTable = get_datatable(datatable_ptr);
		if (dataTable == nullptr) {
			lastFuncCallError = 1;
			return;
		}

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
		lastFuncCallError = 0;
		DataTable *dataTable = get_datatable(datatable_ptr);
		if (dataTable == nullptr) {
			lastFuncCallError = 1;
			return;
		}

		objects.erase(datatable_ptr);
		delete dataTable;
	}

	/* BSpline interface */
	/* Constructor */
	obj_ptr bspline_init(obj_ptr datatable_ptr, int degree) {
		lastFuncCallError = 0;
		DataTable *table = get_datatable(datatable_ptr);
		if (table == nullptr) {
			lastFuncCallError = 1;
			return nullptr;
		}

		BSplineType bsplineType;
		switch (degree) {
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
			throw new Exception("Invalid degree passed to BSpline constructor!");
		}
		}

		obj_ptr bspline = (obj_ptr) new BSpline(*table, bsplineType);
		if (bspline == nullptr) {
			lastFuncCallError = 1;
			return nullptr;
		}
		objects.insert(bspline);
		return bspline;
	}

	double bspline_eval(obj_ptr bspline_ptr, double *x, int x_dim) {
		lastFuncCallError = 0;
		BSpline *bspline = get_bspline(bspline_ptr);
		if (bspline == nullptr) {
			lastFuncCallError = 1;
			return 0.0;
		}

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
		lastFuncCallError = 0;
		BSpline *bspline = get_bspline(bspline_ptr);
		if (bspline == nullptr) {
			lastFuncCallError = 1;
			return nullptr;
		}

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
		lastFuncCallError = 0;
		BSpline *bspline = get_bspline(bspline_ptr);
		if (bspline == nullptr) {
			lastFuncCallError = 1;
			return nullptr;
		}

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
		lastFuncCallError = 0;
		BSpline *bspline = get_bspline(bspline_ptr);
		if (bspline == nullptr) {
			lastFuncCallError = 1;
			return;
		}
		objects.erase(bspline_ptr);
		delete bspline;
	}

	/* PSpline interface */
	/* Constructor */
	obj_ptr pspline_init(obj_ptr pspline_ptr, double lambda) {
		lastFuncCallError = 0;
		DataTable *table = get_datatable(pspline_ptr);
		if (table == nullptr) {
			lastFuncCallError = 1;
			return nullptr;
		}

		obj_ptr pspline = (obj_ptr) new PSpline(*table, lambda);
		if (pspline == nullptr) {
			lastFuncCallError = 1;
			return nullptr;
		}
		objects.insert(pspline);
		return pspline;
	}

	double pspline_eval(obj_ptr pspline_ptr, double *x, int x_dim) {
		lastFuncCallError = 0;
		PSpline *pspline = get_pspline(pspline_ptr);
		if (pspline == nullptr) {
			lastFuncCallError = 1;
			return 0.0;
		}

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
		lastFuncCallError = 0;
		PSpline *pspline = get_pspline(pspline_ptr);
		if (pspline == nullptr) {
			lastFuncCallError = 1;
			return nullptr;
		}

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
		lastFuncCallError = 0;
		PSpline *pspline = get_pspline(pspline_ptr);
		if (pspline == nullptr) {
			lastFuncCallError = 1;
			return nullptr;
		}

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
		lastFuncCallError = 0;
		PSpline *pspline = get_pspline(pspline_ptr);
		if (pspline == nullptr) {
			lastFuncCallError = 1;
			return;
		}
		objects.erase(pspline_ptr);
		delete pspline;
	}

	/* RadialBasisFunction interface */
	/* Constructor */
	obj_ptr rbf_init(obj_ptr rbf_ptr, int type_index, int normalized) {
		lastFuncCallError = 0;
		DataTable *table = get_datatable(rbf_ptr);
		if (table == nullptr) {
			lastFuncCallError = 1;
			return nullptr;
		}

		RadialBasisFunctionType type;
		switch (type_index) {
		case 1:
			type = RadialBasisFunctionType::THIN_PLATE_SPLINE;
			break;
		case 2:
			type = RadialBasisFunctionType::MULTIQUADRIC;
			break;
		case 3:
			type = RadialBasisFunctionType::INVERSE_QUADRIC;
			break;
		case 4:
			type = RadialBasisFunctionType::INVERSE_MULTIQUADRIC;
			break;
		case 5:
			type = RadialBasisFunctionType::GAUSSIAN;
			break;
		default:
			type = RadialBasisFunctionType::THIN_PLATE_SPLINE;
			break;
		}

		bool norm = normalized == 0 ? false : true;

		obj_ptr rbf = (obj_ptr) new RadialBasisFunction(*table, type, norm);
		if (rbf == nullptr) {
			lastFuncCallError = 1;
			return nullptr;
		}
		objects.insert(rbf);
		return rbf;
	}

	double rbf_eval(obj_ptr rbf_ptr, double *x, int x_dim) {
		lastFuncCallError = 0;
		RadialBasisFunction *rbf = get_rbf(rbf_ptr);
		if (rbf == nullptr) {
			lastFuncCallError = 1;
			return 0.0;
		}

		DenseVector vec(x_dim);
		for (int i = 0; i < x_dim; i++) {
			vec(i) = x[i];
		}

		return rbf->eval(vec);
	}

	/*
	Return double ptr so MatLab identifies the libpointer type as
	'doublePtr'. Then we can do reshape(ptr, x_dim, y_dim) in MatLab
	to get a MatLab matrix.
	*/
	double *rbf_eval_jacobian(obj_ptr rbf_ptr, double *x, int x_dim) {
		lastFuncCallError = 0;
		RadialBasisFunction *rbf = get_rbf(rbf_ptr);
		if (rbf == nullptr) {
			lastFuncCallError = 1;
			return nullptr;
		}

		DenseVector vec(x_dim);
		for (int i = 0; i < x_dim; i++) {
			vec(i) = x[i];
		}

		DenseMatrix jacobian = rbf->evalJacobian(vec);

		int numCoefficients = jacobian.cols() * jacobian.rows();
		double *res = (double *)malloc(sizeof(double) * numCoefficients);
		memcpy(res, jacobian.data(), sizeof(double) * numCoefficients);
		return res;
	}

	double *rbf_eval_hessian(obj_ptr rbf_ptr, double *x, int x_dim) {
		lastFuncCallError = 0;
		RadialBasisFunction *rbf = get_rbf(rbf_ptr);
		if (rbf == nullptr) {
			lastFuncCallError = 1;
			return nullptr;
		}

		DenseVector vec(x_dim);
		for (int i = 0; i < x_dim; i++) {
			vec(i) = x[i];
		}

		DenseMatrix hessian = rbf->evalHessian(vec);

		// The DenseMatrix lives on the stack, so we need to copy
		// the data into a newly allocated memory area.
		int numCoefficients = hessian.cols() * hessian.rows();
		double *res = (double *)malloc(sizeof(double) * numCoefficients);
		memcpy(res, hessian.data(), sizeof(double) * numCoefficients);
		return res;
	}

	void rbf_delete(obj_ptr rbf_ptr) {
		lastFuncCallError = 0;
		RadialBasisFunction *rbf = get_rbf(rbf_ptr);
		if (rbf == nullptr) {
			lastFuncCallError = 1;
			return;
		}
		objects.erase(rbf_ptr);
		delete rbf;
	}
}
