import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as fits
from jwst import datamodels
from jwst.datamodels import CubeModel

import pdb







# this function implements the outlier detection algorithm for TSOs with JWST

def outlier_det(model, ndiff=2, box_size=10, sigma=5.0, verbose=False):
	
	''' This function implements the median difference image algorithm on Webb time series observations for detection of outliers. The algorithm follows broadly the method of Nikolov et al 2014 developed for HST.
	
	Parameters
	----------
	- model: 		the input data model. this should be a JWST CubeModel, of a spectrophotometrically calibrated time series. (the _calints.file) [CubeModel]
	- ndiff:			the number of images before and after that will be used for the master median difference image. default value = 2 [integer]
	- box_size:			the size of the pixel box for the pixel spatial median method. default value = 10 [integer]
	- sigma:			the cutoff threshold for outlier flagging. default value = 5.0 [float]
	- verbose:			if you want extra text output, default = False [Boolean]
	
	'''
	
	
	# Validate the input
	# * is a valid CubeModel?
	# * is a time series?
	# * is calibrated?
	assert isinstance(model, CubeModel), "Input model is not a CubeModel"
	assert model.meta.visit.tsovisit == True, "Input model is not a TSO"
	assert model.meta.cal_step.photom == 'COMPLETE', "Input model has not been calibrated"
	
	# Validate the function parameters:
	# * check that (2*ndiff) + 1 is less than the number of integrations
	# * ndiff should be bigger than 1
	# * check that the box size is smaller than the size of the subarray passed
	# 
	assert model.meta.exposure.nints >= ((2. * ndiff) +1.), "Ndiff is too large compared to the number of integrations.Choose a smaller number"
	assert ndiff >= 1, "Ndiff should be >= 1"
	assert box_size < np.shape(model.data)[1], "Pixel box size incompatible with size of data array"
	assert box_size < np.shape(model.data)[2], "Pixel box size incompatible with size of data array"
	
	# Extract relevant parameters
	# * Number of integrations
	
	nint = model.meta.exposure.nints
	if not isinstance(ndiff, int):
		ndiff = np.int(np.round(ndiff))
	diff_ims = np.zeros((nint, np.shape(model.data)[1], np.shape(model.data)[2]))
	diffs = np.zeros((2*ndiff+1, np.shape(model.data)[1], np.shape(model.data)[2]))
	
	flags = np.ones(nint, dtype=int)
	
	# Construct master difference image:
	# for each image I in the cube, create difference of I-ndiff, I+ndiff etc
 	
	for i in range(nint):
		dd = []
		print('Diffs for image {}'.format(i))
		for ii in range(-ndiff, ndiff+1):
			if (i+ii < 0) | (ii==0): 
				if verbose:
					print('i = {0}, ii = {1}, subtraction out of range and skipping'.format(i, ii))
				continue
			if i+ii > nint-1:
				if verbose:
					print('i = {0}, ii = {1}, subtraction out of range and skipping'.format(i, ii))
				continue
			if verbose:
				print('Subtracting image {0} from image {1}'.format(i+ii, i))
			diff = model.data[i,:,:] - model.data[i+ii, :, :]
			dd.append(diff)
		if verbose:
			print('Number of valid diffs: {0}'.format(len(dd)))
		
		# if the length of the dd list is < 3 then set the flag to -1 to mark that these data are potentially of lower quality
		if len(dd) < 3:
			flags[i] = 0
		
		# conver the list into a 3-D array
		tmp_arr = np.asarray(dd)
		diff_ims[i,:,:] = np.median(tmp_arr, axis=0)
		
	
	
	fig, ax = plt.subplots(ncols=10, nrows=1, figsize=[20,8])
	
	for i in range(nint):
		ax[i].imshow(diff_ims[i,:,:], origin='lower', interpolation='None')
	
	fig.show()
	
	print(flags)
	
	return diff_ims
	
	
	