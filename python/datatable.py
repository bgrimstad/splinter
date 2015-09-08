# This file is part of the SPLINTER library.
# Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.


import SPLINTER
from ctypes import *


class DataTable:
	def __init__(self, fileName=None):
		self.__handle = None # Handle to instance in the library
		self.xDim = None
		self.numSamples = 0 # Number of samples not yet transferred to back end
		self.samples = []
		
		if fileName is not None:
			self.load(fileName)
		else:
			self.__handle = SPLINTER.call(SPLINTER.getHandle().datatable_init)
	
	def load(self, fileName):
		PY3 = sys.version_info[0] == 3
		if PY3:
			self.__handle = SPLINTER.call(SPLINTER.getHandle().datatable_load, 0, c_char_p(fileName.encode("UTF-8")))
		else:
			self.__handle = SPLINTER.call(SPLINTER.getHandle().datatable_load, 0, c_char_p(fileName))
	
	def addSample(self, x, y):
		if self.xDim is None:
			self.xDim = len(x)
		
		if self.xDim != len(x):
			raise Exception("Dimension of the new sample disagrees with the dimension of previous samples!\nPrevious: " + str(self.xDim) + ", new: " + str(len(x)))
			
		self.samples += list(x)
		self.samples += [y]
		self.numSamples += 1
	
	# Transfer samples to the library
	def finish(self):
		#print("Transferring " + str(self.numSamples) + " samples to backend:")
		#for i in range(self.numSamples):
		#	print(str(self.samples[i*(self.xDim+1)]) + "," + str(self.samples[i*(self.xDim+1)+1]) + " = " + str(self.samples[i*(self.xDim+1)+2]))
		
		SPLINTER.call(SPLINTER.getHandle().datatable_add_samples_row_major, self.__handle, (c_double * len(self.samples))(*self.samples), self.numSamples, self.xDim)
		
		self.samples = []
		self.numSamples = 0
	
	# Getter for the datatable for use by Approximants
	# Will make sure all samples are transferred to the back end before returning the handle to an Approximant
	def getHandle(self):
		if self.numSamples > 0:
			self.finish()
		
		return self.__handle


if __name__ == "__main__":
	SPLINTER.load("/home/anders/SPLINTER/build/release/libsplinter-matlab-1-4.so")
	#d = DataTable("test.datatable")
	d = DataTable()
	for i in range(10):
		for j in range(10):
			d.addSample([i,j], i*j)
	d.finish()