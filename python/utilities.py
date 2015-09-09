# This file is part of the SPLINTER library.
# Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.


import sys # Python version
from ctypes import *


def isPython3():
	return sys.version_info[0] == 3

def getCString(pyString):
	if isPython3():
		return c_char_p(pyString.encode("UTF-8"))
	else:
		return c_char_p(pyString)

def getPyString(CString):
	return str(CString)

def isString(pyString):
	if isPython3():
		return isinstance(pyString, str)
	else:
		return isinstance(pyString, basestring)

def listToCArrayOfDoubles(pyList):
	return (c_double * len(pyList))(*pyList)

def listToCArrayOfInts(pyList):
	intList = [int(x) for x in pyList]
	return (c_int * len(intList))(*intList)

def CArrayToList(CArray, size):
	return [CArray[i] for i in xrange(size)]
