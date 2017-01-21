/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef SPLINTER_CINTERFACE_H
#define SPLINTER_CINTERFACE_H


#ifndef SPLINTER_API
# ifdef _MSC_VER
#  define SPLINTER_API __declspec(dllexport)
# else
#  define SPLINTER_API
# endif
#endif

// Pointer to C++ objects, passed into the C interface then cast to the correct type.
typedef void *splinter_obj_ptr;


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
SPLINTER_API int splinter_get_error();

/**
 * Get a string describing the error.
 *
 * @return Error string.
 */
SPLINTER_API const char *splinter_get_error_string();





/**
 * Initialize a new DataTable.
 *
 * @return Pointer to the created DataTable.
 */
SPLINTER_API splinter_obj_ptr splinter_datatable_init();

/**
 * Load a datatable from file
 *
 * @param filename Name of the file to load. Must be a datatable that has previously been stored.
 * @return Pointer to the loaded DataTable.
 */
SPLINTER_API splinter_obj_ptr splinter_datatable_load_init(const char *filename);

/**
 * Add samples that are stored in row major order to the datatable.
 *
 * If x0 = [x0_0, x0_1, x0_2], x1 = [x1_0, x1_1, x1_2] are two inputs with
 * y0 = [y0_0, y0_1], y1 = [y1_0, y1_1], as the corresponding outputs, then
 *
 * xs = [x0_0, x0_1, x0_2, x1_0, x1_1, x1_2]
 * x_dim = 3
 * ys = [y0_0, y0_1, y1_0, y1_1]
 * y_dim = 2
 * n_samples = 2
 *
 * @param datatable_ptr Pointer to the datatable.
 * @param xs Pointer to the start of the inputs.
 * @param x_dim The dimension of each input.
 * @param ys Pointer to the start of the outputs.
 * @param y_dim The dimension of each output.
 * @param n_samples Number of samples to add (= xs/x_dim = ys/y_dim).
 */
SPLINTER_API void splinter_datatable_add_samples_row_major(splinter_obj_ptr datatable_ptr,
                                                           double *xs, int x_dim,
                                                           double *ys, int y_dim,
                                                           int n_samples);

/**
 * Add samples that are stored in column major order to the datatable.
 *
 * @param datatable_ptr Pointer to the datatable.
 * @param x Pointer to the start of the samples.
 * @param n_samples Number of samples to add.
 * @param x_dim The dimension of each point (that is, the sample size - 1).
 */
SPLINTER_API void splinter_datatable_add_samples_col_major(splinter_obj_ptr datatable_ptr, double *x, int n_samples, int x_dim);

/**
 * Get the dimension of the domain of the samples in the datatable.
 *
 * @param datatable_ptr Pointer to the datatable.
 * @return The dimension of the domain of the samples in the datatable.
 */
SPLINTER_API int splinter_datatable_get_dim_x(splinter_obj_ptr datatable_ptr);

/**
 * Get the dimension of the codomain of the samples in the datatable.
 *
 * @param datatable_ptr Pointer to the datatable.
 * @return The dimension of the codomain of the samples in the datatable.
 */
SPLINTER_API int splinter_datatable_get_dim_y(splinter_obj_ptr datatable_ptr);

/**
 * Get the number of samples stored in the datatable.
 *
 * @param datatable_ptr Pointer to the datatable.
 * @return The number of samples in the datatable.
 */
SPLINTER_API int splinter_datatable_get_num_samples(splinter_obj_ptr datatable_ptr);

/**
 * Save the datatable to file.
 *
 * @param datatable_ptr Pointer to the datatable.
 * @param filename The file to store the datatable to (will be overwritten!).
 */
SPLINTER_API void splinter_datatable_save(splinter_obj_ptr datatable_ptr, const char *filename);

/**
 * Free the memory of a datatable.
 *
 * @param datatable_ptr Pointer to the datatable.
 */
SPLINTER_API void splinter_datatable_delete(splinter_obj_ptr datatable_ptr);





/**
 * Create a new BSpline::Builder.
 *
 * @param dim_x Dimension of domain (number of inputs)
 * @param dim_y Dimension of codomain (number of outputs)
 * @return Pointer to the created BSplineBuilder.
 */
SPLINTER_API splinter_obj_ptr splinter_bspline_builder_init(int dim_x, int dim_y);

/**
 * Set the degree of the BSplineBuilder.
 *
 * @param bspline_builder_ptr The BSplineBuilder to set the degree of.
 * @param degrees Array of degrees (must be of the same dimension as the BSplineBuilder).
 * @param n Dimension of degrees
 */
SPLINTER_API void splinter_bspline_builder_set_degree(splinter_obj_ptr bspline_builder_ptr, unsigned int *degrees, int n);

/**
 * Set the number of basis functions of the BSplineBuilder.
 *
 * @param bspline_builder_ptr The BSplineBuilder to set the number of basis function for.
 * @param num_basis_functions Array of numbers denoting the number of basis functions in corresponding dimensions.
 * @param n Size of num_basis_functions (must match the dimension of bspline_builder_ptr).
 */
SPLINTER_API void splinter_bspline_builder_set_num_basis_functions(splinter_obj_ptr bspline_builder_ptr, int *num_basis_functions, int n);

/**
 * Set the knot spacing of the BSplineBuilder.
 *
 * @param bspline_builder_ptr The BSplineBuilder to set the knot spacing of.
 * @param knot_spacing The knot spacing (actually an enum, see the implementation of this function for details).
 */
SPLINTER_API void splinter_bspline_builder_set_knot_spacing(splinter_obj_ptr bspline_builder_ptr, int knot_spacing);

/**
 * Build the BSpline with the parameters of the Builder.
 *
 * @param bspline_builder_ptr The Builder to "build the BSpline with".
 * @param datatable_ptr The datatable to create the BSpline::Builder from.
 * @param smoothing Smoothing type (actually an enum, see the implementation of this function for details)
 * @param alpha Regularization/smoothing parameter (must be non-negative).
 * @return Pointer to the created BSpline.
 */
SPLINTER_API splinter_obj_ptr splinter_bspline_builder_fit(splinter_obj_ptr bspline_builder_ptr,
                                                           splinter_obj_ptr datatable_ptr,
                                                           int smoothing,
                                                           double alpha);

/**
 * Free the memory of the internal Builder
 *
 * @param bspline_builder_ptr Pointer to the Builder.
 */
SPLINTER_API void splinter_bspline_builder_delete(splinter_obj_ptr bspline_builder_ptr);




/**
 * Construct a BSpline from parameters: coefficients, knot vectors and degrees.
 *
 * @param dim_x The number of variables
 * @param dim_y The number of outputs
 * @param control_points The B-spline coefficients
 * @param num_control_points Number of coefficients
 * @param knot_vectors The B-spline knot vectors
 * @param num_knots_per_vector Number of knots per knot vector
 * @param degrees The B-spline degrees
 * @return Pointer to the created BSpline.
 */
SPLINTER_API splinter_obj_ptr splinter_bspline_param_init(int dim_x, int dim_y, double *control_points,
                                                          int num_control_points, double *knot_vectors,
                                                          int *num_knots_per_vector, unsigned int *degrees);

/**
 * Load a BSpline from file.
 *
 * @param filename The file to load the BSpline from.
 * @return Pointer to the loaded BSpline.
 */
SPLINTER_API splinter_obj_ptr splinter_bspline_load_init(const char *filename);

/**
 * Get the sizes of the knot vectors that are returned by splinter_bspline_get_knot_vectors
 *
 * @param bspline_ptr Pointer to the BSpline
 * @return Array of splinter_bspline_get_num_variables length
 */
SPLINTER_API int *splinter_bspline_get_knot_vector_sizes(splinter_obj_ptr bspline_ptr);

/**
 * Get the knot vectors of the BSpline
 * Two-dimensional matrix where the rows are the knot vectors
 * There are splinter_bspline_get_num_variables knot vectors
 * The sizes of the knot vectors are given by splinter_bspline_get_knot_vector_sizes
 *
 * @param bspline_ptr Pointer to the BSpline
 * @return Row major array of size stated above
 */
SPLINTER_API double *splinter_bspline_get_knot_vectors(splinter_obj_ptr bspline_ptr);

/**
 * Get the number of coefficients of the BSpline
 *
 * @param bspline_ptr Pointer to the BSpline
 */
SPLINTER_API int splinter_bspline_get_num_control_points(splinter_obj_ptr bspline_ptr);

/**
 * Get the control points of the BSpline
 * Returns a two-dimensional matrix with
 * splinter_bspline_get_num_control_points rows
 * and splinter_bspline_get_num_outputs columns
 *
 * @param bspline_ptr Pointer to the BSpline
 * @return Row major flattened array of the control points
 */
SPLINTER_API double *splinter_bspline_get_control_points(splinter_obj_ptr bspline_ptr);

/**
 * Get the knot averages of the BSpline
 * Returns a two-dimensional matrix with
 * splinter_bspline_get_num_control_points rows
 * and splinter_bspline_get_num_variables columns
 *
 * @param bspline_ptr Pointer to the BSpline
 * @return Row major flattened array of the knot averages
 */
SPLINTER_API double *splinter_bspline_get_knot_averages(splinter_obj_ptr bspline_ptr);

/**
 * Get the basis degrees of the BSpline
 * Returns an array of size splinter_bspline_get_num_variables
 *
 * @param bspline_ptr Pointer to the BSpline
 * @return Array of basis degrees
 */
SPLINTER_API int *splinter_bspline_get_basis_degrees(splinter_obj_ptr bspline_ptr);

/**
 * Evaluate a BSpline in one or more points. Can "batch evaluate" several points by storing the points consecutively in
 * x in row major order. If function is a two-dimensional function you can evaluate the function in x0 = [0, 1] and
 * x1 = [2, 3] in one function call by having x point to the start of an array of four doubles: 0 1 2 3, and x_len = 4.
 * This function will then return an array of 2 doubles, the first being the result of evaluating the function in [0, 1],
 * and the second being the result of evaluating the function in [2, 3].
 *
 * If the function is multi-dimensional in the codomain, the results are flattened so that the result of evaluating
 * the first point comes first, then the second, etc.
 *
 * @param bspline_ptr Pointer to the BSpline to evaluate.
 * @param xs Array of doubles. Is of x_len length.
 * @param x_len Length of x.
 * @return Array of results corresponding to the points in x.
 */
SPLINTER_API double *splinter_bspline_eval_row_major(splinter_obj_ptr bspline_ptr, double *xs, int x_len);

/**
 * Evaluate the jacobian of a BSpline in one or more points.
 * @see splinter_bspline_eval_row_major() for further explanation of the behaviour.
 *
 * @param bspline_ptr Pointer to the BSpline to evaluate.
 * @param x Array of doubles. Is of x_len length.
 * @param x_len Length of x.
 * @return Flattened array of array of results corresponding to the points in x.
 */
SPLINTER_API double *splinter_bspline_eval_jacobian_row_major(splinter_obj_ptr bspline_ptr, double *x, int x_len);

/**
 * Evaluate the a BSpline in one or more points that are stored in column major order.
 * @see splinter_bspline_eval_row_major() for further explanation of the behaviour.
 *
 * @param bspline_ptr Pointer to the BSpline to evaluate.
 * @param x Array of doubles. Is of x_len length.
 * @param x_len Length of x.
 * @return Array of results.
 */
SPLINTER_API double *splinter_bspline_eval_col_major(splinter_obj_ptr bspline_ptr, double *x, int x_len);

/**
 * Evaluate the jacobian of a BSpline in one or more points that are stored in column major order.
 * @see splinter_bspline_eval_row_major() for further explanation of the behaviour.
 *
 * @param bspline_ptr Pointer to the BSpline to evaluate.
 * @param x Array of doubles. Is of x_len length.
 * @param x_len Length of x.
 * @return Flattened array of array of results corresponding to the points in x. Stored in ROW major order.
 */
SPLINTER_API double *splinter_bspline_eval_jacobian_col_major(splinter_obj_ptr bspline_ptr, double *x, int x_len);

/**
 * Get the number of inputs of a BSpline (dimension of domain).
 *
 * @param bspline_ptr Pointer to the BSpline.
 */
SPLINTER_API int splinter_bspline_get_dim_x(splinter_obj_ptr bspline_ptr);

/**
 * Get the number of outputs of a BSpline (dimensions of codomain).
 *
 * @param bspline_ptr Pointer to the BSpline.
 */
SPLINTER_API int splinter_bspline_get_dim_y(splinter_obj_ptr bspline_ptr);

/**
 * Save a BSpline to file.
 *
 * @param bspline_ptr Pointer to the BSpline.
 * @param filename File to save the BSpline to (will be overwritten!).
 */
SPLINTER_API void splinter_bspline_save(splinter_obj_ptr bspline_ptr, const char *filename);

/**
 * Free the memory used by a BSpline.
 *
 * @param bspline_ptr Pointer to the BSpline.
 */
SPLINTER_API void splinter_bspline_delete(splinter_obj_ptr bspline_ptr);

/**
 * Insert knots at tau of multiplicity 'multiplicity' to knot vector in variable 'dim'. The B-spline is geometrically
 * unaltered by the knot insertion.
 *
 * @param bspline_ptr Pointer to the BSpline to evaluate.
 * @param tau Knot to insert
 * @param dim Knot vector to insert knots
 * @param multiplicity Desired multiplicity of knot
 */
SPLINTER_API void splinter_bspline_insert_knots(splinter_obj_ptr bspline_ptr, double tau, unsigned int dim, unsigned int multiplicity);

/**
 * Insert knots until all knots have multiplicity equal to the B-spline degree.
 *
 * @param bspline_ptr Pointer to the BSpline to decompose to Bezier form
 */
SPLINTER_API void splinter_bspline_decompose_to_bezier_form(splinter_obj_ptr bspline_ptr);

/**
 * Make a copy of a BSpline.
 *
 * @param bspline_ptr
 * @return The new copy
 */
SPLINTER_API splinter_obj_ptr splinter_bspline_copy(splinter_obj_ptr bspline_ptr);

#ifdef __cplusplus
    }
#endif

#endif // SPLINTER_CINTERFACE_H
