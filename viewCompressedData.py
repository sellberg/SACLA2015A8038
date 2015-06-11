#!/home/doniach/dermen/epd731/bin/python

import numpy as np
import h5py as h
import matplotlib
import matplotlib.pyplot as plt
import gc
import csv, sys, os, re, shutil, subprocess, time
from optparse import OptionParser

"""
# Example:
#    ./viewCompressedData.py -f 259782.h5
# For details, type:
#    ./viewCompressedData.py -h
"""

parser = OptionParser()
parser.add_option("-f", "--inputFile", action="store", type="string", dest="inputFile", help="Input HDF5 file produced by dataCompress3.py", metavar="FILENAME", default="")
parser.add_option("-d", "--detector", action="store", type="string", dest="detector", help="Detector name that should be compressed (default: all)", metavar="DETECTORNAME", default="")
parser.add_option("-r", "--run", action="store", type="int", dest="runNumber", help="Run number to compress (default: all in input HDF5 file)", metavar="RUNNUMBER", default=0)
parser.add_option("-n", "--nbins", action="store", type="int", dest="nBins", help="Number of motor positions to bin (default: 1 = no binning)", metavar="BINS", default=1)
parser.add_option("-m", "--motor", action="store", type="string", dest="motor", help="Scan motor for optical delay (default: PM_EH2_1)", metavar="MOTORNAME", default="PM_EH2_1")
parser.add_option("-v", "--verbose", action="store_true", dest="verbose", help="Output additional information to screen", default=False)

(options, args) = parser.parse_args()

##############################################
# SET NECESSARY ABSOLUTE PATHS ON CLUSTER HERE
# --------------------------------------------
# Make sure to add tailing slashes '/'
# after the end of the path!
##############################################
csv_dir = "/work/fperakis/csv/scan_output/"
#csv_dir = ""
source_dir = "/UserData/fperakis/"
#source_dir = ""
working_dir = "/work/fperakis/output/"
#working_dir = ""
original_dir = os.getcwd() + '/'

if not os.path.exists(source_dir + options.inputFile):
    source_dir = working_dir

# PARAMETERS
defaultEphoton = 11015.336 # in eV
defaultWavelength = 1.12 # in A
defaultDetectorDistance = 0.0625 # in m
defaultPixelSize = 50E-6 # in m
defaultSystemGain = 18.1497 # in electrons/ADU
# center coordinates from Dermens AgBe ring fitting
xc = 1199.582 # in pixels
yc = 1200.499 # in pixels
# calculate y and x coordinates of detector
#x = np.array([np.arange(2399) for i in np.arange(2399)])
#y = x.transpose
#x -= xc
#y -= yc
#rad = np.sqrt(x*x + y*y)
#qRad = 2.*np.pi*2.*np.sin(0.5*np.arctan2(rad*pixelSize, detectorDistance))/wavelength

def radialProfile(image, center, mask=None, wavelength=1.12, detectorDistance=0.0625, pixelSize=50E-6):
    """
    Compute the radial intensity profile of `image`.

    Parameters
    ----------
    image : np.ndarray
        The image to perform the radial profile on.

    center : tuple of floats (x, y)
        The center to use, in pixel units

    Returns
    -------
    q : np.ndarray
        The momentum transfer (A-1)

    Iq
        The average intensities as a function of q.
    """

    # compute the radii for image[y, x]
    x = np.arange(image.shape[1])
    y = np.arange(image.shape[0])
    xx, yy = np.meshgrid(x, y)
    xx -= center[0]
    yy -= center[1]
    rad = np.sqrt(xx*xx + yy*yy)
    # convert to momentum transfer in A-1
    qRad = 2.*np.pi*2.*np.sin(0.5*np.arctan2(rad*pixelSize, detectorDistance))/wavelength
    assert qRad.shape == image.shape
    # histogram the intensities and normalize by number of pixels in each bin to obtain average intensity
    nBins = max(image.shape)/2
    if mask is None:
        bin_values, bin_edges = np.histogram(qRad, weights=image, bins=nBins)
        bin_normalizations = np.histogram(qRad, bins=bin_edges)
    else:
        bin_values, bin_edges = np.histogram(qRad[mask], weights=image[mask], bins=nBins)
        bin_normalizations = np.histogram(qRad[mask], bins=bin_edges)
    Iq = bin_values/bin_normalizations[0]
    q = np.array([(bin_edges[i] + bin_edges[i+1])/2 for i in range(len(bin_values))])
    # run into 'High memory usage error', try to delete
    del x, y
    del xx, yy
    del rad, qRad
    del bin_values, bin_edges, bin_normalizations
    return q, Iq

def photonConvert(image, Ephoton=11000, Gsys=18, eps=3.65):
    """
    Convert the 2D `image` to photons.

    Parameters
    ----------
    image : np.ndarray
        The image to perform the conversion on.

    Ephoton : np.float
        Photon energy in eV

    Gsys : np.float
        System gain in electrons/ADU [default: 18]

    Eps : np.float
        Energy to create electron-hole pair in Si in eV/electron [default: 3.65]

    Returns
    -------
    gain
        The effective gain of the detector in ADU/photon.
    """

    # compute the effective gain [ADU/photon]
    gain = Ephoton/(Gsys*eps)
    # convert image to photons
    image /= gain
    # image is an object that is passed by reference, no need to return
    return gain

########################################################
# Imaging class copied from Ingrid Ofte's pyana_misc code
########################################################
class img_class(object):
	def __init__(self, inarr, filename):
		#self.inarr = inarr*(inarr>0)
		self.inarr = inarr[:]
		self.filename = filename
		self.cmax = self.inarr.max()
		self.cmin = self.inarr.min()
	
	def on_keypress(self,event):
		if event.key == 'p':
			if not os.path.exists(write_dir + runtag):
				os.mkdir(write_dir + runtag)
			pngtag = write_dir + runtag + "/%s.png" % (self.filename)	
			print "saving image as " + pngtag 
			plt.savefig(pngtag)
		if event.key == 'r':
			plt.clim(self.cmin, self.cmax)
			plt.draw()

	def on_click(self, event):
		if event.inaxes:
			lims = self.axes.get_clim()
			colmin = lims[0]
			colmax = lims[1]
			range = colmax - colmin
			value = colmin + event.ydata * range
			if event.button is 1 :
				if value > colmin and value < colmax :
					colmin = value
			elif event.button is 2 :
				colmin, colmax = self.orglims
			elif event.button is 3 :
				if value > colmin and value < colmax:
					colmax = value
			plt.clim(colmin, colmax)
			plt.draw()

	def draw_img(self):
		fig = plt.figure()
		cid1 = fig.canvas.mpl_connect('key_press_event', self.on_keypress)
		cid2 = fig.canvas.mpl_connect('button_press_event', self.on_click)
		canvas = fig.add_subplot(111)
		canvas.set_title(self.filename)
		#plt.rc('image',origin='lower',interpolation='nearest')
		##self.axes = plt.imshow(self.inarr[:,::-1], vmin = 0, vmax = self.cmax)
		self.axes = plt.imshow(self.inarr, interpolation='nearest', vmin=0, vmax=self.cmax)
		#if not canvas.xaxis_inverted():
		#	canvas.invert_xaxis() 
		self.colbar = plt.colorbar(self.axes, pad=0.01)
		self.orglims = self.axes.get_clim()
                self.printInstructions()
		plt.show() 

	def draw_map(self):
		fig = plt.figure()
		cid1 = fig.canvas.mpl_connect('key_press_event', self.on_keypress)
		cid2 = fig.canvas.mpl_connect('button_press_event', self.on_click)
		canvas = fig.add_subplot(111)
		canvas.set_title(self.filename)
		#plt.rc('image',origin='lower',interpolation='nearest')
		##self.axes = plt.imshow(self.inarr[:,::-1], vmin = 0, vmax = self.cmax)
		#self.axes = plt.imshow(self.inarr, interpolation='nearest', vmin=self.cmin, vmax=self.cmax, extent=[-1,1,-1,1])
                print "scale in Y-direction is inverted, need to change order of extent!!"
                aspectRatio = (np.max(q_ref) - np.min(q_ref))/(np.max(detectorMotorValues) - np.min(detectorMotorValues))
		self.axes = plt.imshow(self.inarr, interpolation='nearest', vmin=self.cmin, vmax=self.cmax, extent=[np.min(q_ref),np.max(q_ref),np.max(detectorMotorValues),np.min(detectorMotorValues)], aspect=aspectRatio)
		#if not canvas.xaxis_inverted():
		#	canvas.invert_xaxis() 
		self.colbar = plt.colorbar(self.axes, pad=0.01)
		self.orglims = self.axes.get_clim()
                self.printInstructions()
		plt.show() 

        def printInstructions(self):
            print "Right-click on colorbar to set maximum scale."
            print "Left-click on colorbar to set minimum scale."
            print "Center-click on colorbar (or press 'r') to reset color scale."
            print "Interactive controls for zooming at the bottom of figure screen (zooming..etc)."
            print "Press 'p' to save PNG of image (with the current colorscales) in the appropriately named folder."
            print "Hit Ctl-\ or close all windows (Alt-F4) to terminate viewing program."


if options.inputFile != '' and os.path.exists(source_dir + options.inputFile):
    print "Reading data from '%s'..." % (options.inputFile)
    times = []
    tic = time.time()
    # open input h5 file
    f = h.File(source_dir + options.inputFile, "r")
    runList = f.keys()
    print "\tFound %d run(s) in file:" % len(runList)
    print "\t\t" + ', '.join(runList)
    if options.runNumber:
        if ("run_%d" % options.runNumber in runList):
            # reduce runList to only the run chosen
            runList = ["run_%d" % options.runNumber]
        else:
            print "Could not find run %d in '%s', aborting..." % (options.runNumber, options.inputFile)
            sys.exit(1)
    toc = time.time()
    times.append(toc - tic)
    runIndex = 0
    for r in runList:
        print "Reading averages for %s..." % r
        runKeys = np.array(f[r].keys())
        #outputRunGroup = o.require_group(inputPath)
        if options.detector != '':
             if (runKeys == options.detector).any():
                 # reduce runKeys to only the detector chosen
                 runKeys = np.array([options.detector])
             else:
                 print "\tDetector '%s' does not exist in '%s', aborting..." % (options.detector, inputPath)
                 sys.exit(1)
        # tags
        runDataTags = []
        runReferenceTags = []
        runBackgroundTags = []
        # scan motor value
        runDataMotorValues = []
        detectorIndex = 0
        for d in runKeys:
            if "detector" in d:
                print "\tReading detector data for %s..." % d
                p = r + '/' + d
                detectorKeys = np.array(f[p].keys())
                if (options.motor + '_scan') in detectorKeys:
                    scanIndex = 0
                    scanRun = True
                    detectorMotorValues = []
                    detectorCorrectedDataMean = []
                    q_data = []
                    I_data = []
                    I_data_norm = []
                    I_div = []
                    I_diff = []
                    pm = p + '/' + options.motor + '_scan'
                    for n in f[pm].keys():
                        pt = pm + '/' + n + '/motor_value'
                    	if (f.get(pt)):
                    	    detectorMotorValues.append(f[pt].value)
                    	else:
                    	    print "Could not find dataset '%s' in '%s', aborting..." % (pt, options.inputFile)
                    	    sys.exit(1)
                        pt = pm + '/' + n + '/analysis'
                    	if (f.get(pt)):
                        	detectorCorrectedDataMean.append(np.array(f[pt + "/correctedImage"]))
                    	else:
                    	    print "Could not find datagroup '%s' in '%s', aborting..." % (pt, options.inputFile)
                    	    sys.exit(1)
                        pt = pm + '/' + n + '/analysis/angularAverage'
                    	if (f.get(pt)):
                            q_data.append(np.array(f[pt + "/Q"]))
                            I_data.append(np.array(f[pt + "/raw"]))
                            I_data_norm.append(np.array(f[pt + "/normalized"]))
                            I_div.append(np.array(f[pt + "/divided"]))
                            I_diff.append(np.array(f[pt + "/normalizedDifference"]))
                    	else:
                    	    print "Could not find datagroup '%s' in '%s', aborting..." % (pt, options.inputFile)
                    	    sys.exit(1)
                        runDataMotorValues.append(detectorMotorValues)
                        scanIndex +=1
                elif 'data' in detectorKeys:
                    pt = p + '/data/analysis'
                    scanRun = False
                    if (f.get(pt)):
                        detectorCorrectedDataMean = np.array(f[pt + "/correctedImage"])
                    else:
                        print "Could not find datagroup '%s' in '%s', aborting..." % (pt, options.inputFile)
                        sys.exit(1)
                    pt = p + '/data/analysis/angularAverage'
                    if (f.get(pt)):
                        q_data = np.array(f[pt + "/Q"])
                        I_data = np.array(f[pt + "/raw"])
                        I_data_norm = np.array(f[pt + "/normalized"])
                        I_div = np.array(f[pt + "/divided"])
                        I_diff = np.array(f[pt + "/normalizedDifference"])
                    else:
                        print "Could not find datagroup '%s' in '%s', aborting..." % (pt, options.inputFile)
                        sys.exit(1)
                else:
                    print "Could not find scan motor %s or averaged data set, aborting..."  % (options.motor)
                    sys.exit(1)
                
                if 'reference' in detectorKeys:
                    pt = p + '/reference/analysis'
                    if (f.get(pt)):
                        detectorCorrectedReferenceMean = np.array(f[pt + "/correctedImage"])
                    else:
                        print "Could not find datagroup '%s' in '%s', aborting..." % (pt, options.inputFile)
                        sys.exit(1)
                    pt = p + '/reference/analysis/angularAverage'
                    if (f.get(pt)):
                        q_ref = np.array(f[pt + "/Q"])
                        I_ref = np.array(f[pt + "/raw"])
                        I_ref_norm = np.array(f[pt + "/normalized"])
                    else:
                        print "Could not find datagroup '%s' in '%s', aborting..." % (pt, options.inputFile)
                        sys.exit(1)
                # plot data
                if scanRun:
                    #cmap = matplotlib.cm.seismic
                    cmap = matplotlib.cm.gist_rainbow
                    scanLegend = []
                    runDataMotorDiffStdev = []
                    runDataMotorDivStdev = []
                    plt.figure(1)
                    for n in range(len(runDataMotorValues[runIndex])):
                        plt.plot(q_data[n], I_diff[n], color=cmap(n / float(len(runDataMotorValues[runIndex]))))
                        qRangeToSum = np.where((q_data[n] > 1.5) & (q_data[n] < 3.5))
                        runDataMotorDiffStdev.append(I_diff[n][qRangeToSum].std())
                        scanLegend.append("%.2f" % (runDataMotorValues[runIndex][n]))
                    plt.title(r + ', ' + d + ', ' + options.motor + ' difference')
                    plt.xlabel('Q [A-1]')
                    plt.ylabel('I(Q) [arb. units]')
                    plt.legend(scanLegend)

                    plt.figure(2)
                    for n in range(len(runDataMotorValues[runIndex])):
                        plt.plot(q_data[n], I_div[n], color=cmap(n / float(len(runDataMotorValues[runIndex]))))
                        qRangeToSum = np.where((q_data[n] > 1.5) & (q_data[n] < 3.5))
                        runDataMotorDivStdev.append(I_div[n][qRangeToSum].std())
                    plt.title(r + ', ' + d + ', ' + options.motor + ' division')
                    plt.xlabel('Q [A-1]')
                    plt.ylabel('I(Q) [arb. units]')
                    plt.legend(scanLegend)

                    scanLegend.append('no laser')
                    plt.figure(3)
                    for n in range(len(runDataMotorValues[runIndex])):
                        plt.plot(q_data[n], I_data[n], color=cmap(n / float(len(runDataMotorValues[runIndex]))))
                    plt.plot(q_ref, I_ref, 'k')
                    plt.title(r + ', ' + d + ', ' + options.motor + ' raw')
                    plt.xlabel('Q [A-1]')
                    plt.ylabel('I(Q) [photons/pixel]')
                    plt.legend(scanLegend)
                    
                    plt.figure(4)
                    for n in range(len(runDataMotorValues[runIndex])):
                        plt.plot(q_data[n], I_data_norm[n], color=cmap(n / float(len(runDataMotorValues[runIndex]))))
                    plt.plot(q_ref, I_ref_norm, 'k')
                    plt.title(r + ', ' + d + ', ' + options.motor + ' normalized')
                    plt.xlabel('Q [A-1]')
                    plt.ylabel('I(Q) [arb. units]')
                    plt.legend(scanLegend)

                    #plt.figure(5)
                    #print (np.abs(q_data[0] - 3.5)).argmin()
                    #I_diff = np.array(I_diff)
                    #plt.plot(runDataMotorValues[runIndex], I_diff[:, 878], 'b')
                    #plt.plot(runDataMotorValues[runIndex], I_diff[:, 747], 'g')
                    #plt.plot(runDataMotorValues[runIndex], I_diff[:, 643], 'r')
                    #plt.title(r + ' time dependence')
                    #plt.xlabel("%s [pulses]" % options.motor)
                    #plt.ylabel('I(Q)-difference [arb. units]')
                    #plt.legend(['3.5 A-1', '3.0 A-1', '2.5 A-1'])
                    
                    #plt.figure(6)
                    #I_diff_before = I_diff[:6, :].mean(axis=0)
                    #I_diff_time0 = I_diff[6, :]
                    #I_diff_after = I_diff[7:, :].mean(axis=0)
                    #plt.plot(q_data[0], I_diff_before, 'b')
                    #plt.plot(q_data[6], I_diff_time0, 'g')
                    #plt.plot(q_data[7], I_diff_after, 'r')
                    #plt.title(r + ' intensity difference')
                    #plt.xlabel("Q [A-1]")
                    #plt.ylabel('I(Q)-difference [arb. units]')
                    #plt.legend(['x-rays before laser', 'x-rays at laser', 'x-rays after laser'])
                    
                    #plt.figure(7)
                    #plt.plot(runDataMotorValues[runIndex], runDataMotorDiffStdev, 'bx')
                    #plt.plot(runDataMotorValues[runIndex], runDataMotorDivStdev, 'ro')
                    #plt.title(r + ', ' + d + ', ' + options.motor + ' change between 1.5 and 3.5 A-1')
                    #plt.xlabel("%s [pulses]" % options.motor)
                    #plt.ylabel('Std. dev. [arb. units]')
                    #plt.legend(['difference', 'division'])
                    
                    i = img_class(np.array(I_diff), r + ' difference map')
                    i.draw_map()
                else:
                    plt.figure(1)
                    plt.plot(q_data, I_diff)
                    plt.title(r + ', ' + d)
                    plt.xlabel('Q [A-1]')
                    plt.ylabel('I(Q) [arb. units]')
                    
                    plt.figure(2)
                    plt.plot(q_data, I_data, 'b')
                    plt.plot(q_data, I_ref, 'r')
                    plt.title(r + ', ' + d)
                    plt.legend(['data', 'reference'])
                    plt.xlabel('Q [A-1]')
                    plt.ylabel('I(Q) [arb. units]')
                    
                    plt.figure(3)
                    plt.plot(q_data, I_data_norm, 'b')
                    plt.plot(q_data, I_ref_norm, 'r')
                    plt.title(r + ', ' + d)
                    plt.legend(['data', 'reference'])
                    plt.xlabel('Q [A-1]')
                    plt.ylabel('I(Q) [photons/pixel]')
                    
                    plt.figure(4)
                    plt.plot(q_data, I_div)
                    plt.title(r + ', ' + d)
                    plt.xlabel('Q [A-1]')
                    plt.ylabel('I(Q) [arb. units]')
                    
                    i = img_class(detectorCorrectedDataMean, r)
                    i.draw_img()
            detectorIndex += 1
        runIndex += 1
    f.close()
else:
    if options.inputFile != '':
        print "Input file '%s' does not exist, aborting..." % options.inputFile
    else:
        print "No input file specified, aborting..."
    sys.exit(1)
