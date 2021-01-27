import numpy as np
import matplotlib.pyplot as plt
import jwst.datamodels





# this function implements the outlier detection algorithm for TSOs with JWST

def tso_outlier_detection(model, ndiff=2, box_size=10, sigma=5.0):
	
	''' This function implements the median difference image algorithm on Webb time series observations for detection of outliers. The algorithm follows broadly the method of Nikolov et al 2014 developed for HST.
	
	Parameters
	----------
	- model: 		the input data model. this should be a JWST CubeModel, of a spectrophotometrically calibrated time series. (the _calints.file) [CubeModel]
	- ndiff:			the number of images before and after that will be used for the master median difference image. default value = 2 [integer]
	- box_size:			the size of the pixel box for the pixel spatial median method. default value = 10 [integer]
	- sigma:			the cutoff threshold for outlier flagging. default value = 5.0 [float]
	
	'''
	
	
	# Validate the input
	# * is a valid CubeModel?
	# * is a time series?
	# * is calibrated?
    assert isinstance(model, datamodels.CubeModel), "Input model is not a CubeModel"
    assert model.meta.visit.tsovisit == True, "Input model is not a TSO"
    assert model.meta.cal_step.photom == 'COMPLETE', "Input model has not been calibrated"

    # Validate the function parameters:
    # * check that (2*ndiff) + 1 is less than the number of integrations
    # * ndiff should be bigger than 1
    # * check that the box size is smaller than the size of the subarray passed
    # 

    assert model.meta.exposure.nints >= ((2. * ndiff) +1.), "Ndiff is too large compared to the number of integrations.Choose a smaller number"
    assert ndiff >= 1, "Ndiff should be >= 1"
    assert box_size < np.shape(model.data)[1] & box_size < np.shape(model.data)[2], "Pixel box size incompatible with size of data array"
	
	
	
	# Extract relevant parameters
	# * Number of integrations
	
	
	
	# Construct master difference image:
	# for each image I in the cube, create difference of I-ndiff, I+ndiff etc
