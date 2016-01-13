/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef SPLINTER_MATLAB_H
#define SPLINTER_MATLAB_H


#ifndef SPLINTER_API
# ifdef _MSC_VER
#  define SPLINTER_API __declspec(dllexport)
# else
#  define SPLINTER_API
# endif
#endif

// Pointer to C++ objects, passed into the C interface then cast to the correct type.
typedef void *obj_ptr;


#ifdef __cplusplus
    extern "C"
    {
#endif
/**
 * Check if the last library call resulted in an error.
 * Will reset upon call, so two consecutive calls to this function may not return the same value.
 *
 * @return 1 if error, 0 else.
 */
SPLINTER_API int get_error();

/**
 * Get a string describing the error.
 *
 * @return Error string.
 */
SPLINTER_API const char *get_error_string();

/**
 * Initialize a new DataTable.
 *
 * @return Pointer to the created DataTable.
 */
SPLINTER_API obj_ptr datatable_init();

/**
 * Load a datatable from file
 *
 * @param filename Name of the file to load. Must be a datatable that has previously been stored.
 * @return Pointer to the loaded DataTable.
 */
SPLINTER_API obj_ptr datatable_load_init(const char *filename);

/**
 * Add samples that are stored in row major order to the datatable.
 *
 * @param datatable_ptr Pointer to the datatable.
 * @param x Pointer to the start of the samples.
 * @param n_samples Number of samples to add.
 * @param x_dim The dimension of each point (that is, the sample size - 1).
 */
SPLINTER_API void datatable_add_samples_row_major(obj_ptr datatable_ptr, double *x, int n_samples, int x_dim);

/**
 * Add samples that are stored in column major order to the datatable.
 *
 * @param datatable_ptr Pointer to the datatable.
 * @param x Pointer to the start of the samples.
 * @param n_samples Number of samples to add.
 * @param x_dim The dimension of each point (that is, the sample size - 1).
 * @param size The size of each sample (that is, x_dim + 1)
 */
SPLINTER_API void datatable_add_samples_col_major(obj_ptr datatable_ptr, double *x, int n_samples, int x_dim, int size);

/**
 * Get the number of variables (dimension of the samples) in the datatable.
 *
 * @param datatable_ptr Pointer to the datatable.
 * @return The number of variables (dimension) of the samples in the datatable.
 */
SPLINTER_API int datatable_get_num_variables(obj_ptr datatable_ptr);

/**
 * Get the number of samples stored in the datatable.
 *
 * @param datatable_ptr Pointer to the datatable.
 * @return The number of samples in the datatable.
 */
SPLINTER_API int datatable_get_num_samples(obj_ptr datatable_ptr);

/**
 * Save the datatable to file.
 *
 * @param datatable_ptr Pointer to the datatable.
 * @param filename The file to store the datatable to (will be overwritten!).
 */
SPLINTER_API void datatable_save(obj_ptr datatable_ptr, const char *filename);

/**
 * Free the memory of a datatable.
 *
 * @param datatable_ptr Pointer to the datatable.
 */
SPLINTER_API void datatable_delete(obj_ptr datatable_ptr);

/**
 * Load a BSpline from file.
 *
 * @param filename The file to load the BSpline from.
 * @return Pointer to the loaded BSpline.
 */
SPLINTER_API obj_ptr bspline_load_init(const char *filename);

/**
 * Create a new RadialBasisFunction.
 *
 * @param datatable_ptr The datatable to create the RBF from.
 * @param type_index The RBF type to use (see the implementation of this function for possible types).
 * @param normalized 0 or 1, whether to normalize the data or not.
 * @return Pointer to the created RBF.
 */
SPLINTER_API obj_ptr rbf_init(obj_ptr datatable_ptr, int type_index, int normalized);

/**
 * Load a RadialBasisFunction from file.
 *
 * @param filename The file to load the RBF from.
 * @return Pointer to the loaded RBF.
 */
SPLINTER_API obj_ptr rbf_load_init(const char *filename);

/**
 * Create a new Polynomial.
 *
 * @param datatable_ptr The datatable to create the Polynomial from.
 * @param degrees Array of degrees (one for each dimension in datatable_ptr).
 * @param degrees_dim Dimension of degrees (should match the dimension of datatable_ptr).
 * @return Pointer to the created Polynomial.
 */
SPLINTER_API obj_ptr polynomial_regression_init(obj_ptr datatable_ptr, int *degrees, int degrees_dim);

/**
 * Load a Polynomial from file.
 *
 * @param filename The file to load the Polynomial from.
 * @return Pointer to the loaded Polynomial.
 */
SPLINTER_API obj_ptr polynomial_regression_load_init(const char *filename);


/**
 * Create a new BSplineBuilder.
 *
 * @param datatable_ptr The datatable to create the BSplineBuilder from.
 * @return Pointer to the created BSplineBuilder.
 */
SPLINTER_API obj_ptr bspline_builder_init(obj_ptr datatable_ptr);

/**
 * Set the degree of the BSplineBuilder.
 *
 * @param bspline_builder_ptr The BSplineBuilder to set the degree of.
 * @param degrees Array of degrees (must be of the same dimension as the BSplineBuilder).
 * @param n Dimension of degrees
 */
SPLINTER_API void bspline_builder_set_degree(obj_ptr bspline_builder_ptr, int *degrees, int n);

/**
 * Set the number of basis functions of the BSplineBuilder.
 *
 * @param bspline_builder_ptr The BSplineBuilder to set the number of basis function for.
 * @param num_basis_functions Array of numbers denoting the number of basis functions in corresponding dimensions.
 * @param n Size of num_basis_functions (must match the dimension of bspline_builder_ptr).
 */
SPLINTER_API void bspline_builder_set_num_basis_functions(obj_ptr bspline_builder_ptr, int *num_basis_functions, int n);

/**
 * Set the knot spacing of the BSplineBuilder.
 *
 * @param bspline_builder_ptr The BSplineBuilder to set the knot spacing of.
 * @param knot_spacing The knot spacing (actually an enum, see the implementation of this function for details).
 */
SPLINTER_API void bspline_builder_set_knot_spacing(obj_ptr bspline_builder_ptr, int knot_spacing);

/**
 * Set the smoothing of the BSplineBuilder.
 *
 * @param bspline_builder_ptr The BSplineBuilder to set the knot spacing of.
 * @param smoothing The smoothing to use (actually an enum, see the implementation of this function for details).
 */
SPLINTER_API void bspline_builder_set_smoothing(obj_ptr bspline_builder_ptr, int smoothing);

/**
 * Set the lambda of the BSplineBuilder.
 *
 * @param bspline_builder_ptr The BSplineBuilder to set the lambda of.
 * @param lambda The new lambda to use (must be non-negative).
 */
SPLINTER_API void bspline_builder_set_lambda(obj_ptr bspline_builder_ptr, double lambda);

/**
 * Build the BSpline with the parameters of the BSplineBuilder.
 *
 * @param bspline_builder_ptr The BSplineBuilder to "build the BSpline with".
 * @return Pointer to the created BSpline.
 */
SPLINTER_API obj_ptr bspline_builder_build(obj_ptr bspline_builder_ptr);

/**
 * Free the memory of the internal BSpline::Builder
 *
 * @param bspline_builder_ptr Pointer to the BSpline::Builder.
 */
SPLINTER_API void bspline_builder_delete(obj_ptr bspline_builder_ptr);

/**
 * Evaluate a function in one or more points. Can "batch evaluate" several points by storing the points consecutively in
 * x in row major order. If function is a two-dimensional function you can evaluate the function in x0 = [0, 1] and
 * x1 = [2, 3] in one function call by having x point to the start of an array of four doubles: 0 1 2 3, and x_len = 4.
 * This function will then return an array of 2 doubles, the first being the result of evaluating the function in [0, 1],
 * and the second being the result of evaluating the function in [2, 3].
 *
 * @param function Pointer to the function to evaluate.
 * @param x Array of doubles. Is of x_len length.
 * @param x_len Length of x.
 * @return Array of results corresponding to the points in x.
 */
SPLINTER_API double *eval_row_major(obj_ptr function, double *x, int x_len);

/**
 * Evaluate the jacobian of a function in one or more points.
 * @see eval_row_major() for further explanation of the behaviour.
 *
 * @param function Pointer to the function to evaluate.
 * @param x Array of doubles. Is of x_len length.
 * @param x_len Length of x.
 * @return Flattened array of array of results corresponding to the points in x.
 */
SPLINTER_API double *eval_jacobian_row_major(obj_ptr function, double *x, int x_len);

/**
 * Evaluate the hessian of a function in one or more points.
 * @see eval_row_major() for further explanation of the behaviour.
 *
 * @param function Pointer to the function to evaluate.
 * @param x Array of doubles. Is of x_len length.
 * @param x_len Length of x.
 * @return Flattened array of array of array of results corresponding to the points in x.
 */
SPLINTER_API double *eval_hessian_row_major(obj_ptr function, double *x, int x_len);

/**
 * Evaluate the a function in one or more points that are stored in column major order.
 * @see eval_row_major() for further explanation of the behaviour.
 *
 * @param function Pointer to the function to evaluate.
 * @param x Array of doubles. Is of x_len length.
 * @param x_len Length of x.
 * @return Array of results.
 */
SPLINTER_API double *eval_col_major(obj_ptr function, double *x, int x_len);

/**
 * Evaluate the jacobian of a function in one or more points that are stored in column major order.
 * @see eval_row_major() for further explanation of the behaviour.
 *
 * @param function Pointer to the function to evaluate.
 * @param x Array of doubles. Is of x_len length.
 * @param x_len Length of x.
 * @return Flattened array of array of results corresponding to the points in x. Stored in ROW major order.
 */
SPLINTER_API double *eval_jacobian_col_major(obj_ptr function, double *x, int x_len);

/**
 * Evaluate the hessian of a function in one or more points that are stored in column major order.
 * @see eval_row_major() for further explanation of the behaviour.
 *
 * @param function Pointer to the function to evaluate.
 * @param x Array of doubles. Is of x_len length.
 * @param x_len Length of x.
 * @return Flattened array of array of array of results corresponding to the points in x. Stored in ROW major order.
 */
SPLINTER_API double *eval_hessian_col_major(obj_ptr function, double *x, int x_len);

/**
 * Get the number of variables (dimension) of a Function.
 *
 * @param function Pointer to the Function.
 */
SPLINTER_API int function_get_num_variables(obj_ptr function);

/**
 * Save a Function to file.
 *
 * @param function Pointer to the Function.
 * @param filename File to save the Function to (will be overwritten!).
 */
SPLINTER_API void function_save(obj_ptr function, const char *filename);

/**
 * Free the memory used by a Function.
 *
 * @param function Pointer to the Function.
 */
SPLINTER_API void function_delete(obj_ptr function);

#ifdef __cplusplus
    }
#endif

#endif // SPLINTER_MATLAB_H