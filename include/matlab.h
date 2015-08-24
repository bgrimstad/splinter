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

// Pointer to C++ objects, passed into the C interface then casted to the correct type.
typedef void *obj_ptr;


#ifdef __cplusplus
	extern "C"
	{
#endif
		/* 1 if the previous function call caused an error, 0 otherwise. */
SPLINTER_API int get_error();

SPLINTER_API const char *get_error_string();

SPLINTER_API obj_ptr datatable_init();

SPLINTER_API obj_ptr datatable_load_init(const char *filename);

SPLINTER_API void datatable_add_samples(obj_ptr datatable_ptr, double *x, int n_samples, int x_dim, int size);

SPLINTER_API unsigned int datatable_get_num_variables(obj_ptr datatable_ptr);

SPLINTER_API unsigned int datatable_get_num_samples(obj_ptr datatable_ptr);

SPLINTER_API void datatable_save(obj_ptr datatable_ptr, const char *filename);

SPLINTER_API obj_ptr datatable_load(obj_ptr datatable_ptr, const char *filename);

SPLINTER_API void datatable_delete(obj_ptr datatable_ptr);


SPLINTER_API obj_ptr bspline_init(obj_ptr datatable_ptr, int type);

SPLINTER_API obj_ptr bspline_load_init(const char *filename);

SPLINTER_API obj_ptr pspline_init(obj_ptr datatable_ptr, double lambda);

SPLINTER_API obj_ptr pspline_load_init(const char *filename);

SPLINTER_API obj_ptr rbf_init(obj_ptr datatable_ptr, int type_index, int normalized);

SPLINTER_API obj_ptr rbf_load_init(const char *filename);

SPLINTER_API obj_ptr polynomial_regression_init(obj_ptr datatable_ptr, int *degrees, int degrees_dim);

SPLINTER_API obj_ptr polynomial_regression_load_init(const char *filename);


SPLINTER_API double eval(obj_ptr approximant, double *x, int x_dim);

SPLINTER_API double *eval_jacobian(obj_ptr approximant, double *x, int x_dim);

SPLINTER_API double *eval_hessian(obj_ptr approximant, double *x, int x_dim);

SPLINTER_API int get_num_variables(obj_ptr approximant);

SPLINTER_API void save(obj_ptr approximant, const char *filename);

SPLINTER_API void load(obj_ptr approximant, const char *filename);

SPLINTER_API void delete_approximant(obj_ptr approximant);

#ifdef __cplusplus
	}
#endif

#endif // SPLINTER_MATLAB_H