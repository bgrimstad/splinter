#!/bin/bash
# This file is part of the SPLINTER library.
# Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

PROJECT_DIR="$(pwd)/.."

NUM_MISSING=0
NUM_OK=0
TO_CHECK="$PROJECT_DIR/CMakeLists.txt $PROJECT_DIR/src $PROJECT_DIR/include $PROJECT_DIR/matlab $PROJECT_DIR/test $PROJECT_DIR/scripts"

function check_license {
	while read p; do
		if ! grep -q "$p" $1
		then
			MISSING_LICENSE="true"
			return
		fi
	done < LICENSE_HEADER
}

function check_entry {
	if [ -f $1 ]; then
		check_file $1
	else
		check_directory $1
	fi
}

function check_file {
#	echo "Checking file: $1"

	MISSING_LICENSE="false"
	check_license $1
	if [ $MISSING_LICENSE == "true" ]
	then
		echo "$(readlink -m $1) is missing the license header"
		NUM_MISSING=$(($NUM_MISSING + 1))
	else
		NUM_OK=$(($NUM_OK + 1))
	fi
}

function check_directory {
#	echo "Checking directory: $1"

	for _ENTRY in $(ls $1)
	do
		ENTRY=$1/$_ENTRY

		check_entry $ENTRY
	done
}

check_license

for __ENTRY in $TO_CHECK
do
	check_entry $__ENTRY
done

echo "$NUM_OK/$(($NUM_OK + $NUM_MISSING)) files has the license header."
