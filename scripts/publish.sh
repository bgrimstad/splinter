#!/usr/bin/env bash
#
# This file is part of the SPLINTER library.
# Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

# Publish SPLINTER to (Test)PyPI
# Publish to TestPyPI: ./publish.sh
# Publish to PyPI: ./publish.sh release

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"
BUILD_DIR="${DIR}/../build/"

MODE_DEBUG="debug"
MODE_RELEASE="release"

MODE=${MODE_DEBUG}
if [ $# -eq 1 ] && [ $1 == ${MODE_RELEASE} ]; then
    echo "Release mode"
    MODE=${MODE_RELEASE}
else
    echo "Debug mode"
fi

rm -r ${BUILD_DIR}
mkdir ${BUILD_DIR} && cd ${BUILD_DIR}

cmake .. -DCMAKE_BUILD_TYPE=release

echo "Building SPLINTER"
make -j$(nproc)

make install

cd splinter-python/
python3 setup.py sdist bdist_wheel

if [ ${MODE} == ${MODE_RELEASE} ]; then
    python3 -m twine upload dist/*
    echo "Finished. To install SPLINTER, run 'pip install splinterpy'"
elif [ ${MODE} == ${MODE_DEBUG} ]; then
    python3 -m twine upload --repository-url https://test.pypi.org/legacy/ dist/*
    echo "Finished. To install SPLINTER, run 'pip install --index-url https://test.pypi.org/simple/ splinterpy'"
else
    echo "Invalid mode: ${MODE}"
fi
