import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as fits
from jwst import datamodels
from jwst.datamodels import CubeModel

import pdb







# this function implements the outlier detection algorithm for TSOs with JWST

def make_diff_im(model, ndiff=2, box_size=10, sigma=5.0, verbose=False):
	
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
		
		# convert the list into a 3-D array and take a median along the first axis
		tmp_arr = np.asarray(dd)
		diff_ims[i,:,:] = np.median(tmp_arr, axis=0)
		
	
	
	fig, ax = plt.subplots(ncols=10, nrows=1, figsize=[20,8])
	
	for i in range(nint):
		ax[i].imshow(diff_ims[i,:,:], origin='lower', interpolation='None', cmap='gray')
	
	fig.show()
	
	print(flags)
	
	return diff_ims
	
	

def find_outliers(im, dq=None, sigma=4.0, mode='box', box_size=None, dispersion_axis=None, verbose=False, plot=False):
	
	'''
	This function will find outliers in image im
	
	Parameters
	----------
	- im:		the input image. must be a 2D array of floats.
	- dq:		an array of data quality flags (the DQ array) (optional). If no DQ array is provided the algorithm will make no use of existing flags.
	- mode:		the spatial method used for outlier detection. options are: 
					* 'box' (default): pixel value will be compared with a square box of box_size width x length
					* 'line': pixel will be compared against the full pixel row/column, with direction defined to be perpendicular to the dispersion_axis.
	- box_size:	the size of the box to which each pixel will be compared (integer; required if mode='box'). 
	- dispersion_axis:	the direction along which the spectrum is dispersed. this should correspond to the axis in the input array. (integer; required if mode='line')
	- verbose:	for added text output (boolean; default=False)
	- plot:		for visual output of outlier pixels (bolean; default=False) 
	- 
	'''
	
	# change number formats if needed
	if box_size:
		if not isinstance(box_size, int):
			box_size = np.int(np.round(box_size))
			print('Converting box size to integer value: {0}'.format(box_size))
	
	if dispersion_axis:
		if not isinstance(dispersion_axis, int):
			dispersion_axis = np.int(np.round(dispersion_axis))
			print('Converting dispersion axis to integer value: {0}'.format(dispersion_axis))
	
	# checks to start:
	# - input image is a 2D array
	# - if a DQ array is provided, that the array sizes match
	# - if mode is box, that a box size is given and box size is smaller than the smallest array dimension
	# - if mode is line, that the dispersion_axis is given
	# - that the dispersion axis is 0 or 1
	assert(np.ndim(im)) == 2, "Input image is not a 2D array"
	if dq:
		assert(np.shape(im) == np.shape(dq)), "Data and DQ array shapes do not match"
	if (mode == 'box'):
		assert box_size is not None, "You must specify a box size with mode 'box' "
		assert box_size <= np.min(np.shape(im))
	if (mode == 'line'):
		assert dispersion_axis is not None, "You must specify a dispersion axis with mode 'line' "
		assert dispersion_axis < 2, "Dispersion axis not compatible with array shape"
		assert dispersion_axis > 0, "Dispersion axis cannot be negative"
		
	
	return im
	
	
	# if plot == True, set up the output image
	
	
	
	
	# 
	