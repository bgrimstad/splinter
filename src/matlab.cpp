#include "matlab.h"
#include <datatable.h>
#include <bspline.h>
#include <pspline.h>
#include <radialbasisfunction.h>
#include <polynomialregression.h>
#include <generaldefinitions.h>
#include <set>

using namespace SPLINTER;

// 1 if the last function call caused an error, 0 else
int lastFuncCallError = 0;

const char *error_string = "No error.";

static void set_error_string(const char *new_error_string) {
	error_string = new_error_string;
	lastFuncCallError = 1;
}

// Keep a list of objects so we avoid performing operations on objects that don't exist
std::set<obj_ptr> objects = std::set<obj_ptr>();

/* Cast the obj_ptr to a DataTable * */
static DataTable *get_datatable(obj_ptr datatable_ptr) {
	lastFuncCallError = 0;
	if (objects.count(datatable_ptr) > 0) {
		return (DataTable *)datatable_ptr;
	}

	set_error_string("Invalid reference to DataTable: Maybe it has been deleted?");

	return nullptr;
}

/* Cast the obj_ptr to an Approximant * */
static Approximant *get_approximant(obj_ptr approximant_ptr) {
	lastFuncCallError = 0;
	if (objects.count(approximant_ptr) > 0) {
		return (Approximant *)approximant_ptr;
	}

	set_error_string("Invalid reference to Approximant: Maybe it has been deleted?");

	return nullptr;
}

static DenseVector get_densevector(double *x, int x_dim) {
	DenseVector xvec(x_dim);
	for (int i = 0; i < x_dim; i++) {
		xvec(i) = x[i];
	}

	return xvec;
}

extern "C"
{
	/* 1 if the last call to the library resulted in an error,
	*  0 otherwise. This is used to avoid throwing exceptions across library boundaries,
	*  and we expect the caller to manually check the value of this flag.
	*/
	int get_error() {
		return lastFuncCallError;
	}

	const char *get_error_string() {
		return error_string;
	}

	/* DataTable constructor */
	obj_ptr datatable_init() {
		obj_ptr dataTable = (obj_ptr) new DataTable();

		objects.insert(dataTable);

		return dataTable;
	}

	obj_ptr datatable_load_init(const char *filename) {
		obj_ptr dataTable = (obj_ptr) new DataTable(filename);

		objects.insert(dataTable);

		return dataTable;
	}

	void datatable_add_samples(obj_ptr datatable_ptr, double *x, int n_samples, int x_dim, int size) {
		DataTable *dataTable = get_datatable(datatable_ptr);
		if (dataTable != nullptr) {
			DenseVector vec(x_dim);
			for (int i = 0; i < n_samples; ++i) {
				for (int j = 0; j < x_dim; ++j) {
					vec(j) = x[i + j * size];
				}
				dataTable->addSample(vec, x[i + x_dim * size]);
			}
		}
	}

	unsigned int datatable_get_num_variables(obj_ptr datatable_ptr) {
		DataTable *dataTable = get_datatable(datatable_ptr);
		if (dataTable != nullptr) {
			return dataTable->getNumVariables();
		}

		return 0;
	}

	unsigned int datatable_get_num_samples(obj_ptr datatable_ptr) {
		DataTable *dataTable = get_datatable(datatable_ptr);
		if (dataTable != nullptr) {
			return dataTable->getNumSamples();
		}

		return 0;
	}

	void datatable_save(obj_ptr datatable_ptr, const char *filename) {
		DataTable *dataTable = get_datatable(datatable_ptr);
		if (dataTable != nullptr) {
			dataTable->save(filename);
		}
	}

	// Deletes the previous handle and loads a new
	obj_ptr datatable_load(obj_ptr datatable_ptr, const char *filename) {
		// Delete and reset error, as it will get set if the DataTable didn't exist.
		datatable_delete(datatable_ptr);
		lastFuncCallError = 0;

		obj_ptr dataTable = (obj_ptr) new DataTable(filename);

		objects.insert(dataTable);
		return dataTable;
	}

	void datatable_delete(obj_ptr datatable_ptr) {
		DataTable *dataTable = get_datatable(datatable_ptr);
		if (dataTable != nullptr) {
			objects.erase(datatable_ptr);
			delete dataTable;
		}
	}

	/* BSpline constructor */
	obj_ptr bspline_init(obj_ptr datatable_ptr, int degree) {
		obj_ptr bspline = nullptr;

		auto table = get_datatable(datatable_ptr);
		if (table != nullptr) {
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
					set_error_string("Invalid BSplineType!");
					return nullptr;
				}
			}

			bspline = (obj_ptr) new BSpline(*table, bsplineType);
			objects.insert(bspline);
		}

		return bspline;
	}

	obj_ptr bspline_load_init(const char *filename) {
		obj_ptr bspline = (obj_ptr) new BSpline(filename);

		objects.insert(bspline);

		return bspline;
	}

	/* PSpline constructor */
	obj_ptr pspline_init(obj_ptr datatable_ptr, double lambda) {
		obj_ptr pspline = nullptr;

		auto table = get_datatable(datatable_ptr);
		if (table != nullptr) {
			pspline = (obj_ptr) new PSpline(*table, lambda);
			objects.insert(pspline);
		}

		return pspline;
	}

	obj_ptr pspline_load_init(const char *filename) {
		obj_ptr pspline = (obj_ptr) new PSpline(filename);

		objects.insert(pspline);

		return pspline;
	}

	/* RadialBasisFunction constructor */
	obj_ptr rbf_init(obj_ptr datatable_ptr, int type_index, int normalized) {
		obj_ptr rbf = nullptr;

		auto table = get_datatable(datatable_ptr);
		if (table != nullptr) {
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

			bool norm = normalized != 0;

			rbf = (obj_ptr) new RadialBasisFunction(*table, type, norm);
			objects.insert(rbf);
		}

		return rbf;
	}

	obj_ptr rbf_load_init(const char *filename) {
		obj_ptr rbf = (obj_ptr) new RadialBasisFunction(filename);

		objects.insert(rbf);

		return rbf;
	}

	/* PolynomialRegression constructor */
	obj_ptr polynomial_regression_init(obj_ptr datatable_ptr, int *degrees, int degrees_dim) {
		obj_ptr polyfit = nullptr;

		auto table = get_datatable(datatable_ptr);
		if (table != nullptr) {
			auto degreeVec = std::vector<unsigned int>(degrees_dim);
			for(int i = 0; i < degrees_dim; ++i) {
				degreeVec.at(i) = (unsigned int) degrees[i];
			}

			polyfit = (obj_ptr) new PolynomialRegression(*table, degreeVec);
			objects.insert(polyfit);
		}

		return polyfit;
	}

	obj_ptr polynomial_regression_load_init(const char *filename) {
		obj_ptr polyfit = (obj_ptr) new PolynomialRegression(filename);

		objects.insert(polyfit);

		return polyfit;
	}


	double eval(obj_ptr approximant, double *x, int x_dim) {
		double retVal = 0.0;

		auto approx = get_approximant(approximant);
		if(approx != nullptr) {
			auto xvec = get_densevector(x, x_dim);
			retVal = approx->eval(xvec);
		}

		return retVal;
	}

	double *eval_jacobian(obj_ptr approximant, double *x, int x_dim) {
		double *retVal = nullptr;

		auto approx = get_approximant(approximant);
		if(approx != nullptr) {
			auto xvec = get_densevector(x, x_dim);
			DenseMatrix jacobian = approx->evalJacobian(xvec);

			/* Copy jacobian from stack to heap */
			int numCoefficients = jacobian.cols() * jacobian.rows();
			retVal = (double *) malloc(sizeof(double) * numCoefficients);
			memcpy(retVal, jacobian.data(), sizeof(double) * numCoefficients);
		}

		return retVal;
	}

	double *eval_hessian(obj_ptr approximant, double *x, int x_dim) {
		double *retVal = nullptr;

		auto approx = get_approximant(approximant);
		if(approx != nullptr) {
			auto xvec = get_densevector(x, x_dim);
			DenseMatrix hessian = approx->evalHessian(xvec);

			/* Copy jacobian from stack to heap */
			int numCoefficients = hessian.cols() * hessian.rows();
			retVal = (double *) malloc(sizeof(double) * numCoefficients);
			memcpy(retVal, hessian.data(), sizeof(double) * numCoefficients);
		}

		return retVal;
	}

	int get_num_variables(obj_ptr approximant) {
		int retVal = 0;

		auto approx = get_approximant(approximant);
		if(approx != nullptr) {
			retVal = approx->getNumVariables();
		}

		return retVal;
	}

	void save(obj_ptr approximant, const char *filename) {
		auto approx = get_approximant(approximant);
		if(approx != nullptr) {
			approx->save(filename);
		}
	}

	void load(obj_ptr approximant, const char *filename) {
		auto approx = get_approximant(approximant);
		if(approx != nullptr) {
			approx->load(filename);
		}
	}

	void delete_approximant(obj_ptr approximant) {
		auto approx = get_approximant(approximant);

		if(approx != nullptr) {
			objects.erase(approximant);
			delete approx;
		}
	}
}
