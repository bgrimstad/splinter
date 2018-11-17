# This file is part of the SPLINTER library.
# Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

import os  # Path manipulation for finding files
import ctypes
import platform  # Detect OS
import typing as ty

from .utilities import *


class SplinterBackend:
    def __init__(self):
        # Handle to the cdll instance
        self._handle = None

        full_path = os.path.realpath(__file__)  # Path to this file
        splinter_python_main_dir = os.path.dirname(full_path)

        # Locate version file. If we cannot find it then we won't be able to find splinter either
        version_file = os.path.join(splinter_python_main_dir, "version")
        if not os.path.exists(version_file):
            raise Exception(
                f"Missing version file in SPLINTER directory! "
                f"Please notify the developers of SPLINTER of this error."
            )

        with open(version_file) as f:
            splinter_version = f.read().strip()
        self.version = splinter_version

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
                raise FileNotFoundError(
                    "Unable to automatically locate SPLINTER.\n"
                    "It is possible that SPLINTER was not compiled for the operating system and architecture you are "
                    "using. If that is the case, you can compile SPLINTER yourself and load it using "
                    "'splinterpy.load(\"/path/to/SPLINTER.so\")'"
                )
        else:
            if not os.path.exists(lib_file):
                raise FileNotFoundError(
                    "Unable to load SPLINTER from {}: File does not exist!".format(lib_file)
                )

        try:
            self._handle = ctypes.cdll.LoadLibrary(lib_file)
            self._set_function_signatures()

        except Exception as e:
            out("Error:")
            out("Either you are trying to load a library with another architecture (32bit/64bit) ")
            out("than the Python you are using, ", True)
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
            raise Exception("The SPLINTER back-end has not been loaded!\n"
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
            raise Exception("Got exception when calling {}: {}".format(function.__name__, errorMsg))

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
            function = getattr(this._handle, function_name)

            setattr(function, 'restype', return_type)
            setattr(function, 'argtypes', list(parameters))

        set_signature('splinter_get_error', c_int)
        set_signature('splinter_get_error_string', c_char_p)

        set_signature('splinter_datatable_init', handle_type)
        set_signature('splinter_datatable_add_samples_row_major', c_void, handle_type, c_double_p, c_int, c_double_p, c_int, c_int)
        set_signature('splinter_datatable_get_dim_x', c_int, handle_type)
        set_signature('splinter_datatable_get_dim_y', c_int, handle_type)
        set_signature('splinter_datatable_get_num_samples', c_int, handle_type)
        set_signature('splinter_datatable_delete', c_void, handle_type)

        set_signature('splinter_bspline_from_param', handle_type, c_int, c_int, c_int_p, c_double_p, c_int_p, c_double_p, c_int)
        set_signature('splinter_bspline_from_param_zero', handle_type, c_int, c_int, c_int_p, c_double_p, c_int_p)
        set_signature('splinter_bspline_get_knot_vector_sizes', c_int_p, handle_type)
        set_signature('splinter_bspline_get_knot_vectors', c_double_p, handle_type)
        set_signature('splinter_bspline_get_num_control_points', c_int, handle_type)
        set_signature('splinter_bspline_get_control_points', c_double_p, handle_type)
        set_signature('splinter_bspline_get_knot_averages', c_double_p, handle_type)
        set_signature('splinter_bspline_get_basis_degrees', c_int_p, handle_type)
        set_signature('splinter_bspline_eval_row_major', c_double_p, handle_type, c_double_p, c_int)
        set_signature('splinter_bspline_eval_jacobian_row_major', c_double_p, handle_type, c_double_p, c_int)
        set_signature('splinter_bspline_get_dim_x', c_int, handle_type)
        set_signature('splinter_bspline_get_dim_y', c_int, handle_type)
        set_signature('splinter_bspline_to_json', c_void, handle_type, c_char_p)
        set_signature('splinter_bspline_from_json', handle_type, c_char_p)
        set_signature('splinter_bspline_delete', c_void, handle_type)
        set_signature('splinter_bspline_insert_knots', c_void, handle_type, c_double, c_int, c_int)
        set_signature('splinter_bspline_decompose_to_bezier_form', c_void, handle_type)
        set_signature('splinter_bspline_copy', handle_type, handle_type)
        set_signature('splinter_bspline_fit', handle_type, handle_type, handle_type, c_int, c_double, c_double_p, c_int)

        set_signature('splinter_bspline_interpolator', handle_type, handle_type, c_int)
        set_signature('splinter_bspline_smoother', handle_type, handle_type, c_int, c_int, c_double, c_double_p, c_int)
        set_signature('splinter_bspline_unfitted', handle_type, handle_type, c_int_p, c_int, c_int, c_int_p, c_int)

    def _locate_splinter(self) -> ty.Optional[str]:
        is_linux = platform.system() == 'Linux'
        is_windows = platform.system() == 'Windows'
        is_mac = platform.system() == 'Darwin'

        full_path = os.path.realpath(__file__)  # Path to this file
        splinter_python_main_dir = os.path.dirname(full_path)

        splinter_basename = "splinter-" + self.version
        if is_linux:
            operating_system = "linux"
            splinter_name = "lib" + splinter_basename + ".so"
        elif is_windows:
            operating_system = "windows"
            splinter_name = splinter_basename + ".dll"
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
