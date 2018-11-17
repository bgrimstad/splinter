# This file is part of the SPLINTER library.
# Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

from setuptools import setup
import os
import os.path
import sys

# Publishing to PyPI:
# mkdir build && cd build
# rm -r *
# cmake .. -DCMAKE_BUILD_TYPE=release
# make
# make install
# cd splinter-python
# python3 setup.py sdist bdist_wheel
# Upload to TestPyPI: python -m twine upload --repository-url https://test.pypi.org/legacy/ dist/*
# Upload to PyPI: python -m twine upload --repository-url https://www.pypi.org/legacy/ dist/*
# Install from TestPyPI: python3 -m pip install --index-url https://test.pypi.org/simple/ splinterpy

# Patch version of the Python interface.
# SPLINTER itself should get a patch version, but for now we only add a patch version to the Python interface
# A drawback is that running splinterpy.__version__ returns only the major and minor version numbers.
PYTHON_INTERFACE_PATCH_VERSION = "4" or "0"

try:
    int(PYTHON_INTERFACE_PATCH_VERSION)
except ValueError:
    print(f"Invalid patch version!")
    exit(1)

version_file_name = 'version'  # Name of the file where the C++ back-end version is written
interface_package_name = 'splinterpy'  # Both the name of the project and the name of the package
library_files_dir_name = 'lib'  # Path to the compiled library files

# Get the absolute path of the directory containing this file
setup_py_dir_path = os.path.dirname(os.path.realpath(__file__))
os.chdir(setup_py_dir_path)  # In case this script is being run from another directory
splinter_python_interface_path = os.path.join(setup_py_dir_path, interface_package_name)

version_file_path = os.path.join(setup_py_dir_path, interface_package_name, version_file_name)
library_files_dir_path = os.path.join(setup_py_dir_path, interface_package_name, library_files_dir_name)

# Verify that the library has been built
# If these files don't exist, it means that the `make install` step hasn't been run.
if not (os.path.exists(version_file_path) and os.path.exists(library_files_dir_path)):
    print("Error: It seems that SPLINTER has not been built.\n\n"
          + "Please note that to install the Python interface of SPLINTER from source,\n"
          + "you first need to build the library. See the guide at TODO for more information.")
    exit(1)

# Read the C++ back-end version number.
with open(version_file_path, 'r') as version_file:
    splinter_version_string = version_file.read()

splinter_major_version, splinter_minor_version = splinter_version_string.split('-')


python_version_string = '{}.{}.{}'.format(
    splinter_major_version,
    splinter_minor_version,
    PYTHON_INTERFACE_PATCH_VERSION
)


def get_available_backends(lib_path):
    """
    Gets a list of the available back-ends by looking into the lib directory.
    If there is at least one file in the directory that is supposed to contain a library file, this function
    assumes that the library is available.
    :param lib_path: Path to the lib directory
    :return: List of tuples: [(os, architecture), ...]
    """
    lib_versions_available = []
    for operating_system in os.listdir(lib_path):
        operating_system_path = os.path.join(library_files_dir_path, operating_system)
        if os.path.isdir(operating_system_path):
            for arch in os.listdir(os.path.join(lib_path, operating_system)):
                arch_path = os.path.join(operating_system_path, arch)
                if os.path.isdir(arch_path):
                    lib_versions_available.append((operating_system, arch))
    return lib_versions_available


def get_backend_binary_path(lib_path, operating_system, arch):
    dir_path = os.path.join(lib_path, operating_system, arch)
    files = os.listdir(dir_path)
    # Shouldn't be more than one file per directory
    return os.path.join(dir_path, files[0])


# Make sure all required files will be installed by adding them to the package data list
# Add version file to the package data list
package_data = [version_file_path]
# Add all available library files to the package data list
for operating_system, arch in get_available_backends(library_files_dir_path):
    package_data.append(get_backend_binary_path(library_files_dir_path, operating_system, arch))

# If the action is to upload a package, provide a list of the available back-ends so the user can verify that
# all the wanted back-ends are included
if len(sys.argv) >= 3 and sys.argv[2] == "upload":
    print("Included back-ends:")
    oss = {}
    for operating_system, arch in get_available_backends(library_files_dir_path):
        if operating_system not in oss:
            oss[operating_system] = []
        oss[operating_system].append(arch)
    for operating_system, archs in oss.items():
        print("{}:".format(operating_system))
        for arch in archs:
            print("\t{}".format(arch))
    answer = input("Proceed? (y/n): ")
    if answer.lower() != "y":
        exit(0)

setup(
    name=interface_package_name,
    version=python_version_string,
    description='SPLINTER is a library for multivariate function approximation using splines',
    url='http://github.com/bgrimstad/splinter',
    author='Bjarne Grimstad and others',
    author_email='bjarne.grimstad@gmail.com',
    license='MPL 2.0',
    keywords=['spline', 'interpolation', 'b-spline', 'approximation'],
    packages=[interface_package_name],
    # This seems to work fine, so we keep it on
    zip_safe=True,
    package_data={interface_package_name: package_data},
    install_requires=['numpy'],
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: Mozilla Public License 2.0 (MPL 2.0)',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering :: Mathematics',
        'Topic :: Scientific/Engineering :: Artificial Intelligence',
        'Topic :: Software Development :: Libraries :: Python Modules',
    ],
)
