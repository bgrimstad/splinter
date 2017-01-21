# This file is part of the SPLINTER library.
# Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

import os  # Path manipulation for
import platform  # Detect OS
from .utilities import *
import ctypes


class SplinterBackend:
    def __init__(self):
        # Handle to the cdll instance
        self._handle = None

        self.__version__ = "Not loaded"

    def is_loaded(self) -> bool:
        """
        Check if the back-end of SPLINTER has been loaded
        :return:
        """
        return self._handle is not None

    def load(self, lib_file: str=None) -> bool:
        """
        Attempt to load the shared library back-end from the file in lib_file. If lib_file is None,
        then we will try to automatically locate it.

        Raises FileNotFoundError if we were unable to find the library automatically, or if lib_file does not exist.
        Raises
        :param lib_file: Path to attempt to load SPLINTER from.
        :return: True if SPLINTER was loaded, False otherwise
        """
        if self.is_loaded():
            return True

        if lib_file is None:
            lib_file = self._locate_splinter()
            if lib_file is None:
                raise FileNotFoundError("Unable to automatically locate SPLINTER.\n"\
                + "You can load it manually by doing splinter_py.load(\"/path/to/SPLINTER.so\")")
        else:
            if not os.path.exists(lib_file):
                raise FileNotFoundError("Unable to load SPLINTER from {}: File does not exist!".format(lib_file))

        try:
            self._handle = ctypes.cdll.LoadLibrary(lib_file)
            self._set_function_signatures()
            out("Loaded SPLINTER from " + str(lib_file) + "!")

        except Exception as e:
            out("Error:")
            out("Either you are trying to load a library with another architecture (32bit/64bit) than the Python you are using, ", True)
            out("or the file you are trying to load ({}) could not be found.".format(lib_file))
            out("For reference your Python is " + str(8*ctypes.sizeof(ctypes.c_void_p)) + "bit.")
            out("Here is the error message:")
            out(str(e))
            self._handle = None

    @property
    def handle(self):
        """
        Get the handle to the cdll instance.
        Raises Exception if the SPLINTER back-end has not been loaded.
        :return:
        """
        if self._handle is None:
            raise Exception("The SPLINTER back-end has not been loaded!\n"\
                    + "You can do it with splinter_py.load(\"/path/to/libsplinter-x-y.so\")")
        return self._handle

    def call(self, function, *args):
        """
        Make a call to the C++ back-end
        :param function: What function to call
        :param args: Arguments to the function
        :return: Return value of the function
        """
        res = function(*args)

        if self.handle.splinter_get_error():
            # TODO: Sometimes the string is correct, sometimes not. Investigate.
            errorMsg = get_py_string(self.handle.splinter_get_error_string())
            raise Exception(errorMsg)

        return res

    # Set expected argument types and return types of all functions
    def _set_function_signatures(self):
        """
        Set the return types and argument types of all functions as declared in the C interface.
        :return:
        """
        # Define types for int* and double*
        c_int = ctypes.c_int
        c_double = ctypes.c_double
        c_void_p = ctypes.c_void_p
        c_char = ctypes.c_char
        c_char_p = ctypes.c_char_p

        c_int_p = ctypes.POINTER(c_int)
        c_double_p = ctypes.POINTER(c_double)
        handle_type = c_void_p
        c_void = None

        this = self
        def set_signature(function_name: str, return_type, *parameters):
            this._handle[function_name].restype = return_type
            this._handle[function_name].argtypes = list(*parameters)

        self._handle.splinter_get_error.restype = c_int
        self._handle.splinter_get_error.argtypes = []

        self._handle.splinter_get_error_string.restype = c_char_p
        self._handle.splinter_get_error_string.argtypes = []


        self._handle.splinter_datatable_init.restype = handle_type
        self._handle.splinter_datatable_init.argtypes = []

        self._handle.splinter_datatable_load_init.restype = handle_type
        self._handle.splinter_datatable_load_init.argtypes = [c_char_p]

        self._handle.splinter_datatable_add_samples_row_major.restype = c_void
        self._handle.splinter_datatable_add_samples_row_major.argtypes = [handle_type, c_double_p, c_int, c_double_p, c_int, c_int]

        self._handle.splinter_datatable_get_dim_x.restype = c_int
        self._handle.splinter_datatable_get_dim_x.argtypes = [handle_type]

        self._handle.splinter_datatable_get_dim_y.restype = c_int
        self._handle.splinter_datatable_get_dim_y.argtypes = [handle_type]

        self._handle.splinter_datatable_get_num_samples.restype = c_int
        self._handle.splinter_datatable_get_num_samples.argtypes = [handle_type]

        self._handle.splinter_datatable_save.restype = c_void
        self._handle.splinter_datatable_save.argtypes = [handle_type, c_char_p]

        self._handle.splinter_datatable_delete.restype = c_void
        self._handle.splinter_datatable_delete.argtypes = [handle_type]


        self._handle.splinter_bspline_builder_init.restype = handle_type
        self._handle.splinter_bspline_builder_init.argtypes = [handle_type]

        self._handle.splinter_bspline_builder_set_degree.restype = c_void
        self._handle.splinter_bspline_builder_set_degree.argtypes = [handle_type, c_int_p, c_int]

        self._handle.splinter_bspline_builder_set_num_basis_functions.restype = c_void
        self._handle.splinter_bspline_builder_set_num_basis_functions.argtypes = [handle_type, c_int_p, c_int]

        self._handle.splinter_bspline_builder_set_knot_spacing.restype = c_void
        self._handle.splinter_bspline_builder_set_knot_spacing.argtypes = [handle_type, c_int]

        self._handle.splinter_bspline_builder_set_smoothing.restype = c_void
        self._handle.splinter_bspline_builder_set_smoothing.argtypes = [handle_type, c_int]

        self._handle.splinter_bspline_builder_set_alpha.restype = c_void
        self._handle.splinter_bspline_builder_set_alpha.argtypes = [handle_type, c_double]

        self._handle.splinter_bspline_builder_fit.restype = handle_type
        self._handle.splinter_bspline_builder_fit.argtypes = [handle_type, handle_type]

        self._handle.splinter_bspline_builder_delete.restype = c_void
        self._handle.splinter_bspline_builder_delete.argtypes = [handle_type]

        self._handle.splinter_bspline_param_init.restype = handle_type
        self._handle.splinter_bspline_param_init.argtypes = [c_int, c_int, c_double_p, c_int, c_double_p, c_int_p, c_int_p]

        self._handle.splinter_bspline_load_init.restype = handle_type
        self._handle.splinter_bspline_load_init.argtypes = [c_char_p]

        self._handle.splinter_bspline_get_knot_vector_sizes.restype = c_int_p
        self._handle.splinter_bspline_get_knot_vector_sizes.argtypes = [handle_type]

        self._handle.splinter_bspline_get_knot_vectors.restype = c_double_p
        self._handle.splinter_bspline_get_knot_vectors.argtypes = [handle_type]

        self._handle.splinter_bspline_get_num_control_points.restype = c_int
        self._handle.splinter_bspline_get_num_control_points.argtypes = [handle_type]

        self._handle.splinter_bspline_get_control_points.restype = c_double_p
        self._handle.splinter_bspline_get_control_points.argtypes = [handle_type]

        self._handle.splinter_bspline_get_knot_averages.restype = c_double_p
        self._handle.splinter_bspline_get_knot_averages.argtypes = [handle_type]

        self._handle.splinter_bspline_get_basis_degrees.restype = c_int_p
        self._handle.splinter_bspline_get_basis_degrees.argtypes = [handle_type]

        self._handle.splinter_bspline_eval_row_major.restype = c_double_p
        self._handle.splinter_bspline_eval_row_major.argtypes = [handle_type, c_double_p, c_int]

        self._handle.splinter_bspline_eval_jacobian_row_major.restype = c_double_p
        self._handle.splinter_bspline_eval_jacobian_row_major.argtypes = [handle_type, c_double_p, c_int]

        self._handle.splinter_bspline_get_dim_x.restype = c_int
        self._handle.splinter_bspline_get_dim_x.argtypes = [handle_type]

        self._handle.splinter_bspline_get_dim_y.restype = c_int
        self._handle.splinter_bspline_get_dim_y.argtypes = [handle_type]

        self._handle.splinter_bspline_save.restype = c_void
        self._handle.splinter_bspline_save.argtypes = [handle_type, c_char_p]

        self._handle.splinter_bspline_delete.restype = c_void
        self._handle.splinter_bspline_delete.argtypes = [handle_type]

        self._handle.splinter_bspline_insert_knots.restype = c_void
        self._handle.splinter_bspline_insert_knots.argtypes = [handle_type, c_double, c_int, c_int]

        self._handle.splinter_bspline_decompose_to_bezier_form.restype = c_void
        self._handle.splinter_bspline_decompose_to_bezier_form.argtypes = [handle_type]

        self._handle.splinter_bspline_copy.restype = handle_type
        self._handle.splinter_bspline_copy.argtypes = [handle_type]

    def _locate_splinter(self) -> str:
        is_linux = platform.system() == 'Linux'
        is_windows = platform.system() == 'Windows'
        is_mac = platform.system() == 'Darwin'

        full_path = os.path.realpath(__file__)  # Path to this file
        splinter_python_main_dir = os.path.dirname(full_path)

        # Locate version file. If we cannot find it then we won't be able to find splinter either
        version_file = os.path.join(splinter_python_main_dir, "version")
        if not os.path.exists(version_file):
            return None

        f = open(version_file)
        splinter_version = f.read().strip()
        self.__version__ = splinter_version

        splinter_basename = "splinter-" + splinter_version
        if is_windows:
            operating_system = "windows"
            splinter_name = splinter_basename + ".dll"
        elif is_linux:
            operating_system = "linux"
            splinter_name = "lib" + splinter_basename + ".so"
        elif is_mac:
            operating_system = "osx"
            splinter_name = "lib" + splinter_basename + ".dylib"
        else:
            raise "SPLINTER: Unknown platform: " + platform.system()

        lib_splinter = os.path.join(splinter_python_main_dir, "lib", operating_system, get_architecture(), splinter_name)

        if os.path.exists(lib_splinter):
            return lib_splinter

        return None

splinter_backend_obj = SplinterBackend()
