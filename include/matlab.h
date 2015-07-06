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

/*#define obj_ptr void * */
typedef void *obj_ptr;

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
		/* 1 if the previous function call caused an error, 0 otherwise. */
		API int get_error();

		API const char *get_error_string();

		API obj_ptr datatable_init();

		API obj_ptr datatable_load_init(const char *filename);

		API void datatable_add_samples(obj_ptr datatable_ptr, double *x, int n_samples, int x_dim, int size);

		API unsigned int datatable_get_num_variables(obj_ptr datatable_ptr);

		API unsigned int datatable_get_num_samples(obj_ptr datatable_ptr);

		API void datatable_save(obj_ptr datatable_ptr, const char *filename);

		API obj_ptr datatable_load(obj_ptr datatable_ptr, const char *filename);

		API void datatable_delete(obj_ptr datatable_ptr);


		API obj_ptr bspline_init(obj_ptr datatable_ptr, int type);

		API obj_ptr bspline_load_init(const char *filename);

		API obj_ptr pspline_init(obj_ptr datatable_ptr, double lambda);

		API obj_ptr pspline_load_init(const char *filename);

		API obj_ptr rbf_init(obj_ptr datatable_ptr, int type_index, int normalized);

		API obj_ptr rbf_load_init(const char *filename);

		API obj_ptr polynomial_regression_init(obj_ptr datatable_ptr, int *degrees, int degrees_dim);

		API obj_ptr polynomial_regression_load_init(const char *filename);


		API double eval(obj_ptr approximant, double *x, int x_dim);

		API double *eval_jacobian(obj_ptr approximant, double *x, int x_dim);

		API double *eval_hessian(obj_ptr approximant, double *x, int x_dim);

		API int get_num_variables(obj_ptr approximant);

		API void save(obj_ptr approximant, const char *filename);

		API void load(obj_ptr approximant, const char *filename);

		API void delete_approximant(obj_ptr approximant);

#ifdef __cplusplus
	}
#endif

#endif // SPLINTER_MATLAB_H