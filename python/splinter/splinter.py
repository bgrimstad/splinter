# This file is part of the SPLINTER library.
# Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.


from ctypes import * # Loading libraries
from .utilities import *


# Try to load the library libFile
# If libFile is not specified, we try to locate and load SPLINTER
def load(libFile = None):
	global __handle
	if __handle is not None:
		out("SPLINTER is already loaded!")
		out("If you wish to unload it you must first call unload.")
	else:
		if libFile is None:
			libFile = __locateSplinter()
		if libFile is None:
			raise Exception("Unable to automatically locate SPLINTER.\nYou can load it manually by doing splinter.load(\"/path/to/SPLINTER.so\")")
			
		try:
			__handle = cdll.LoadLibrary(libFile)
			__init()
			out("Loaded SPLINTER!")

		except Exception as e:
			out("Error:")
			out("Either you are trying to load a library with another architecture (32 bit/64 bit) than the Python you are using, ", True)
			out("or the file you are trying to load (" + libFile + ") could not be found.")
			out("For reference your Python is " + str(8*sizeof(c_void_p)) + "bit.")
			out("Here is the error message:")
			out(e)
			__handle = None

def isLoaded():
	global __handle
	return __handle is not None

# TODO:
# Is this function really necessary?
def unload():
	global __handle
	__handle = None


__handle = None


# Below are functions that should only be used internally
def _getHandle():
	global __handle
	
	if __handle is None:
		raise Exception("splinter is not loaded!")
	
	return __handle


# Set expected argument types and return types of all functions
def __init():
	# Define types for int * and double *
	c_int_p = POINTER(c_int)
	c_double_p = POINTER(c_double)
	
	_getHandle().get_error_string.restype = c_char_p
	_getHandle().get_error_string.argtypes = []
	
	_getHandle().datatable_load_init.restype = c_void_p
	_getHandle().datatable_load_init.argtypes = [c_char_p]
	
	_getHandle().datatable_add_samples_row_major.restype = c_void_p
	_getHandle().datatable_add_samples_row_major.argtypes = [c_void_p, c_double_p, c_int, c_int]
	
	_getHandle().eval.restype = c_double_p
	_getHandle().eval.argtypes = [c_void_p, c_double_p, c_int]
	
	_getHandle().eval_jacobian.restype = c_double_p
	_getHandle().eval_jacobian.argtypes = [c_void_p, c_double_p, c_int]
	
	_getHandle().eval_hessian.restype = c_double_p
	_getHandle().eval_hessian.argtypes = [c_void_p, c_double_p, c_int]
	
	_getHandle().polynomial_regression_init.restype = c_void_p
	_getHandle().polynomial_regression_init.argtypes = [c_void_p, c_int_p, c_int]
	
	_getHandle().pspline_load_init.restype = c_void_p
	_getHandle().pspline_load_init.argtypes = [c_char_p]
	
	_getHandle().pspline_init.restype = c_void_p
	_getHandle().pspline_init.argtypes = [c_void_p, c_double]
	
	_getHandle().get_num_variables.restype = c_int
	_getHandle().get_num_variables.argtypes = [c_void_p]
	
	_getHandle().get_num_variables.restype = c_int
	_getHandle().get_num_variables.argtypes = [c_void_p]
	
	_getHandle().datatable_get_num_variables.restype = c_int
	_getHandle().datatable_get_num_variables.argtypes = [c_void_p]


# Try to locate SPLINTER relative to this script
# Assumes the Python interface of splinter has the following directory structure:
# splinter/
# - version.txt
# - approximant.py
# - *.py
# - lib/
#   - libsplinter-x-y.so / splinter-x-y.dll
def __locateSplinter():
	import os
	import platform # Detect OS
	
	isLinux = platform.system() == 'Linux'
	isWindows = platform.system() == 'Windows'
	isMac = platform.system() == 'Darwin'
	
	fullPath = os.path.realpath(__file__) # Path to this file
	splinterPythonMainDir = os.path.dirname(fullPath)
	
	# Locate version.txt. If we cannot find it then we won't be able to find splinter either
	versionFile = os.path.join(splinterPythonMainDir, "version.txt")
	if not os.path.exists(versionFile):
		return None
	
	f = open(versionFile)
	splinterVersion = f.read().strip()
	
	splinterBaseName = "splinter-matlab-" + splinterVersion
	if isWindows:
		splinterName = splinterBaseName + ".dll"
	elif isLinux or isMac:
		splinterName = "lib" + splinterBaseName + ".so"
	else:
		raise("Unknown platform: " + platform.system())
	
	libSplinter = os.path.join(splinterPythonMainDir, "lib", splinterName)
	
	if os.path.exists(libSplinter):
		return libSplinter
	
	return None
	
def _call(function, *args):
	res = function(*args)
	
	if _getHandle().get_error():
		# TODO: Sometimes the string is correct, sometimes not. Investigate.
		errorMsg = getPyString(_getHandle().get_error_string())
		raise Exception(errorMsg)
	
	return res