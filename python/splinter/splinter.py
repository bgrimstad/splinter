# This file is part of the SPLINTER library.
# Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.


from ctypes import *  # Loading libraries
from .utilities import *


# Handle to the cdll instance
__handle = None


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
            out("Loaded SPLINTER from " + str(libFile) + "!")

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


# Below are functions that should only be used internally
def _getHandle():
    global __handle

    if __handle is None:
        raise Exception("splinter is not loaded!")

    return __handle


# Set expected argument types and return types of all functions
# For void functions the restype is None
def __init():
    # Define types for int* and double*
    c_int_p = POINTER(c_int)
    c_double_p = POINTER(c_double)
    handle_type = c_void_p

    _getHandle().splinter_get_error.restype = c_int
    _getHandle().splinter_get_error.argtypes = []

    _getHandle().splinter_get_error_string.restype = c_char_p
    _getHandle().splinter_get_error_string.argtypes = []


    _getHandle().splinter_datatable_init.restype = handle_type
    _getHandle().splinter_datatable_init.argtypes = []

    _getHandle().splinter_datatable_load_init.restype = handle_type
    _getHandle().splinter_datatable_load_init.argtypes = [c_char_p]

    _getHandle().splinter_datatable_add_samples_row_major.restype = None
    _getHandle().splinter_datatable_add_samples_row_major.argtypes = [handle_type, c_double_p, c_int, c_int]

    _getHandle().splinter_datatable_get_num_variables.restype = c_int
    _getHandle().splinter_datatable_get_num_variables.argtypes = [handle_type]

    _getHandle().splinter_datatable_get_num_samples.restype = c_int
    _getHandle().splinter_datatable_get_num_samples.argtypes = [handle_type]

    _getHandle().splinter_datatable_save.restype = None
    _getHandle().splinter_datatable_save.argtypes = [handle_type, c_char_p]

    _getHandle().splinter_datatable_delete.restype = None
    _getHandle().splinter_datatable_delete.argtypes = [handle_type]


    _getHandle().splinter_bspline_builder_init.restype = handle_type
    _getHandle().splinter_bspline_builder_init.argtypes = [handle_type]

    _getHandle().splinter_bspline_builder_set_degree.restype = None
    _getHandle().splinter_bspline_builder_set_degree.argtypes = [handle_type, c_int_p, c_int]

    _getHandle().splinter_bspline_builder_set_num_basis_functions.restype = None
    _getHandle().splinter_bspline_builder_set_num_basis_functions.argtypes = [handle_type, c_int_p, c_int]

    _getHandle().splinter_bspline_builder_set_knot_spacing.restype = None
    _getHandle().splinter_bspline_builder_set_knot_spacing.argtypes = [handle_type, c_int]

    _getHandle().splinter_bspline_builder_set_smoothing.restype = None
    _getHandle().splinter_bspline_builder_set_smoothing.argtypes = [handle_type, c_int]

    _getHandle().splinter_bspline_builder_set_alpha.restype = None
    _getHandle().splinter_bspline_builder_set_alpha.argtypes = [handle_type, c_double]

    _getHandle().splinter_bspline_builder_build.restype = handle_type
    _getHandle().splinter_bspline_builder_build.argtypes = [handle_type]

    _getHandle().splinter_bspline_builder_delete.restype = None
    _getHandle().splinter_bspline_builder_delete.argtypes = [handle_type]

    _getHandle().splinter_bspline_load_init.restype = handle_type
    _getHandle().splinter_bspline_load_init.argtypes = [c_char_p]

    _getHandle().splinter_bspline_get_knot_vector_sizes.restype = c_int_p
    _getHandle().splinter_bspline_get_knot_vector_sizes.argtypes = [handle_type]

    _getHandle().splinter_bspline_get_knot_vectors.restype = c_double_p
    _getHandle().splinter_bspline_get_knot_vectors.argtypes = [handle_type]

    _getHandle().splinter_bspline_get_num_coefficients.restype = c_int
    _getHandle().splinter_bspline_get_num_coefficients.argtypes = [handle_type]

    _getHandle().splinter_bspline_get_coefficients.restype = c_double_p
    _getHandle().splinter_bspline_get_coefficients.argtypes = [handle_type]

    _getHandle().splinter_bspline_get_control_points.restype = c_double_p
    _getHandle().splinter_bspline_get_control_points.argtypes = [handle_type]

    _getHandle().splinter_bspline_eval_row_major.restype = c_double_p
    _getHandle().splinter_bspline_eval_row_major.argtypes = [handle_type, c_double_p, c_int]

    _getHandle().splinter_bspline_eval_jacobian_row_major.restype = c_double_p
    _getHandle().splinter_bspline_eval_jacobian_row_major.argtypes = [handle_type, c_double_p, c_int]

    _getHandle().splinter_bspline_eval_hessian_row_major.restype = c_double_p
    _getHandle().splinter_bspline_eval_hessian_row_major.argtypes = [handle_type, c_double_p, c_int]

    _getHandle().splinter_bspline_get_num_variables.restype = c_int
    _getHandle().splinter_bspline_get_num_variables.argtypes = [handle_type]

    _getHandle().splinter_bspline_save.restype = None
    _getHandle().splinter_bspline_save.argtypes = [handle_type, c_char_p]

    _getHandle().splinter_bspline_delete.restype = None
    _getHandle().splinter_bspline_delete.argtypes = [handle_type]


# Try to locate SPLINTER relative to this script
# Assumes the Python interface of splinter has the following directory structure:
# splinter/
# - version.txt
# - function.py
# - *.py
# - lib/
#   - linux or windows or osx
#     - x86 or x86_64
#       - libsplinter-x-y.so / splinter-x-y.dll
def __locateSplinter():
    import os
    import platform  # Detect OS

    isLinux = platform.system() == 'Linux'
    isWindows = platform.system() == 'Windows'
    isMac = platform.system() == 'Darwin'

    fullPath = os.path.realpath(__file__)  # Path to this file
    # Go two folders up (yes, two is correct, because the first dirname gives us the directory this file resides in).
    splinterPythonMainDir = os.path.dirname(os.path.dirname(os.path.dirname(fullPath)))

    # Locate version.txt. If we cannot find it then we won't be able to find splinter either
    versionFile = os.path.join(splinterPythonMainDir, "version")
    if not os.path.exists(versionFile):
        return None

    f = open(versionFile)
    splinterVersion = f.read().strip()

    splinterBaseName = "splinter-" + splinterVersion
    if isWindows:
        splinterName = splinterBaseName + ".dll"
    elif isLinux or isMac:
        splinterName = "lib" + splinterBaseName + ".so"
    else:
        raise("Unknown platform: " + platform.system())

    operatingSystem = "unknown"
    if isLinux:
        operatingSystem = "linux"
    elif isWindows:
        operatingSystem = "windows"
    elif isMac:
        operatingSystem = "osx"

    libSplinter = os.path.join(splinterPythonMainDir, "lib", operatingSystem, getArchitecture(), splinterName)

    if os.path.exists(libSplinter):
        return libSplinter

    return None


def _call(function, *args):
    res = function(*args)

    if _getHandle().splinter_get_error():
        # TODO: Sometimes the string is correct, sometimes not. Investigate.
        errorMsg = getPyString(_getHandle().splinter_get_error_string())
        raise Exception(errorMsg)

    return res
