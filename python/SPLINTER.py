# This file is part of the SPLINTER library.
# Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.


from ctypes import * # Loading libraries
from utilities import *

__handle = None


def getHandle():
	global __handle
	
	if __handle is None:
		raise Exception("Splinter is not loaded!")
	
	return __handle


# Set expected argument types and return types of all functions
def init():
	global __handle
	global c_double_p
	
	c_int_p = POINTER(c_int)
	c_double_p = POINTER(c_double)
	
	__handle.get_error_string.restype = c_char_p
	__handle.get_error_string.argtypes = []
	
	__handle.datatable_load.restype = c_void_p
	__handle.datatable_load.argtypes = [c_int, c_char_p]
	
	__handle.datatable_add_samples_row_major.restype = c_void_p
	__handle.datatable_add_samples_row_major.argtypes = [c_void_p, c_double_p, c_int, c_int]
	
	__handle.eval.restype = c_double
	__handle.eval.argtypes = [c_void_p, c_double_p, c_int]
	
	__handle.eval_jacobian.restype = c_double_p
	__handle.eval_jacobian.argtypes = [c_void_p, c_double_p, c_int]
	
	__handle.eval_hessian.restype = c_double_p
	__handle.eval_hessian.argtypes = [c_void_p, c_double_p, c_int]
	
	__handle.polynomial_regression_init.restype = c_void_p
	__handle.polynomial_regression_init.argtypes = [c_void_p, c_int_p, c_int]
	
	__handle.pspline_load_init.restype = c_void_p
	__handle.pspline_load_init.argtypes = [c_char_p]
	
	__handle.pspline_init.restype = c_void_p
	__handle.pspline_init.argtypes = [c_void_p, c_double]
	
	__handle.get_num_variables.restype = c_int
	__handle.get_num_variables.argtypes = [c_void_p]

# Try to load the library libFile
def load(libFile):
	global __handle
	if __handle is not None:
		print("Splinter is already loaded!")
		print("If you wish to unload it you must first call unload.")
	else:
		try:
			__handle = cdll.LoadLibrary(libFile)
			init()
			print("Loaded Splinter!")

		except Exception as e:
			print("Error:")
			print("Either you are trying to load a library with another architecture (32 bit/64 bit) than the Python you are using, ",end="")
			print("or the file you are trying to load (" + libFile + ") could not be found.")
			print("For reference your Python is " + str(8*sizeof(c_void_p)) + "bit.")
			print("Here is the error message:")
			print(e)
			__handle = None


# TODO:
# Is this function really necessary?
def unload():
	global __handle
	__handle = None


def call(function, *args):
	res = function(*args)
	
	if getHandle().get_error():
		# TODO: Sometimes the string is correct, sometimes not. Investigate.
		errorMsg = getPyString(getHandle().get_error_string())
		raise Exception(errorMsg)
	
	return res


if __name__ == "__main__":
	load("/home/anders/SPLINTER/build/release/libsplinter-matlab-1-4.so")