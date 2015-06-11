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
parser.add_option("-o", "--outputFile", action="store", type="string", dest="outputFile", help="Output HDF5 file that the compressed data is written to (default: same as input)", metavar="FILENAME", default="")
parser.add_option("-b", "--backgroundFile", action="store", type="string", dest="backgroundFile", help="Input HDF5 file that includes the external background (default: none)", metavar="FILENAME", default="")
parser.add_option("-d", "--detector", action="store", type="string", dest="detector", help="Detector name that should be compressed (default: all)", metavar="DETECTORNAME", default="")
parser.add_option("-r", "--run", action="store", type="int", dest="runNumber", help="Run number to compress (default: all in input HDF5 file)", metavar="RUNNUMBER", default=0)
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
# Planck's constant in eV*s
hPlanck = 4.135667516E-15
# speed of light Angstrom/s
c = 2.99792458E+18
# default parameters
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

if options.inputFile != '' and os.path.exists(source_dir + options.inputFile):
    if (options.outputFile == ''):
        outputFileName = options.inputFile
        working_dir = source_dir
        readWrite = True
    else:
        readWrite = False
    print "Reading data from '%s'..." % (options.inputFile)
    times = []
    tic = time.time()
    # open input h5 file
    if readWrite:
        f = h.File(source_dir + options.inputFile, "a")
    else:
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
        # data
        runDataMean = []
        runDataStd = []
        runDataExcludedMean = []
        runDataExcludedStd = []
        # background
        runBackgroundMean = []
        runBackgroundStd = []
        runBackgroundExcludedMean = []
        runBackgroundExcludedStd = []
        # reference
        runReferenceMean = []
        runReferenceStd = []
        runReferenceExcludedMean = []
        runReferenceExcludedStd = []
        # data temperature
        runTemperatureMean = []
        runTemperatureStd = []
        runTemperatureExcludedMean = []
        runTemperatureExcludedStd = []
        # background temperature
        runBGTemperatureMean = []
        runBGTemperatureStd = []
        runBGTemperatureExcludedMean = []
        runBGTemperatureExcludedStd = []
        # reference temperature
        runRefTemperatureMean = []
        runRefTemperatureStd = []
        runRefTemperatureExcludedMean = []
        runRefTemperatureExcludedStd = []
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
                    tic = time.time()
                    scanIndex = 0
                    scanRun = True
                    detectorDataMean = []
                    detectorDataStd = []
                    detectorMotorValues = []
                    pm = p + '/' + options.motor + '_scan'
                    for n in f[pm].keys():
                        pt = pm + '/' + n + '/data/mean'
                    	if (f.get(pt)):
                    	    detectorDataMean.append(np.array(f[pt]))
                    	else:
                    	    print "Could not find dataset '%s' in '%s', aborting..." % (pt, options.inputFile)
                    	    sys.exit(1)
                        pt = pm + '/' + n + '/data/stdev'
                    	if (f.get(pt)):
                    	    detectorDataStd.append(np.array(f[pt]))
                    	else:
                    	    print "Could not find dataset '%s' in '%s', aborting..." % (pt, options.inputFile)
                    	    sys.exit(1)
                        pt = pm + '/' + n + '/motor_value'
                    	if (f.get(pt)):
                    	    detectorMotorValues.append(f[pt].value)
                    	else:
                    	    print "Could not find dataset '%s' in '%s', aborting..." % (pt, options.inputFile)
                    	    sys.exit(1)
                        scanIndex +=1
                    runDataMean.append(detectorDataMean)
                    runDataStd.append(detectorDataStd)
                    runDataMotorValues.append(np.array(detectorMotorValues))
                    toc = time.time()
                    times.append(toc - tic)
                    print "\t\tRead %d scan positions for motor %s in %.2f s" % (scanIndex, options.motor, toc - tic)
                elif 'data' in detectorKeys:
                    pd = p + '/data/data/mean'
                    scanRun = False
                    if (f.get(pd)):
                        runDataMean.append(np.array(f[pd]))
                    else:
                        print "Could not find dataset '%s' in '%s', aborting..." % (pd, options.inputFile)
                        sys.exit(1)
                    pd = p + '/data/data/stdev'
                    if (f.get(pd)):
                        runDataStd.append(np.array(f[pd]))
                    else:
                        print "Could not find dataset '%s' in '%s', aborting..." % (pd, options.inputFile)
                        sys.exit(1)
                else:
                    print "Could not find scan motor %s or averaged data set, aborting..."  % (options.motor)
                    sys.exit(1)
                
                if 'reference' in detectorKeys:
                    pd = p + '/reference/data/mean'
                    if (f.get(pd)):
                        runReferenceMean.append(np.array(f[pd]))
                    else:
                        print "Could not find dataset '%s' in '%s', aborting..." % (pd, options.inputFile)
                        sys.exit(1)
                    pd = p + '/reference/data/stdev'
                    if (f.get(pd)):
                        runReferenceStd.append(np.array(f[pd]))
                    else:
                        print "Could not find dataset '%s' in '%s', aborting..." % (pd, options.inputFile)
                        sys.exit(1)
                else:
                    print "\tNo reference data included in %s..." % p
                    sys.exit(1)

                if 'background' in detectorKeys:
                    pd = p + '/background/data/mean'
                    if (f.get(pd)):
                        runBackgroundMean.append(np.array(f[pd]))
                    else:
                        print "Could not find dataset '%s' in '%s', aborting..." % (pd, options.inputFile)
                        sys.exit(1)
                    pd = p + '/background/data/stdev'
                    if (f.get(pd)):
                        runBackgroundStd.append(np.array(f[pd]))
                    else:
                        print "Could not find dataset '%s' in '%s', aborting..." % (pd, options.inputFile)
                        sys.exit(1)
                elif (options.backgroundFile != ''):
                    if os.path.exists(source_dir + options.backgroundFile):
                        print "\tUsing external background from '%s'" % options.backgroundFile
                        b = h.File(source_dir + options.backgroundFile, 'r')
                        pd = "run_%s/" % (re.sub('_compressed.h5', '', options.backgroundFile)) + d + '/background/data/mean'
                    	if b.get(pd):
                    	    runBackgroundMean.append(np.array(b[pd]))
                    	else:
                    	    print "Could not find dataset '%s' in '%s', aborting..." % (pd, options.backgroundFile)
                    	    sys.exit(1)
                        pd = "run_%s/" % (re.sub('_compressed.h5', '', options.backgroundFile)) + d + '/background/data/stdev'
                    	if b.get(pd):
                    	    runBackgroundStd.append(np.array(b[pd]))
                    	else:
                    	    print "Could not find dataset '%s' in '%s', aborting..." % (pd, options.backgroundFile)
                    	    sys.exit(1)    
                    else:
                        print "Background file '%s' does not exist, aborting..." % options.backgroundFile
                        sys.exit(1)
                else:
                    print "\tNo background data included in %s, add external background using '-b' flag..." % p
                    sys.exit(1)
                
                if 'detector_info' in detectorKeys:
                    pd = p + '/detector_info/absolute_gain'
                    if (f.get(pd)):
                        systemGain = f[pd].value
                        print "\t\tFound system gain of %.2f electrons/ADU for %s" % (systemGain, d)
                    else:
                        systemGain = defaultSystemGain # in electrons/ADU
                        print "\t\tUsing default system gain of %.2f electrons/ADU" % (systemGain)
                    pd = p + '/detector_info/pixel_size_in_micro_meter'
                    if (f.get(pd)):
                        # assume pixel is square
                        pixelSize = np.mean(f[pd].value)*1E-6 # in m
                        print "\t\tFound pixel size of %.1e m for %s" % (pixelSize, d)
                    else:
                        pixelSize = defaultPixelSize
                        print "\t\tUsing default pixel size of %.1e m" % (pixelSize)
                else:
                    print "\tNo detector info included in %s..." % p
                
                # data corrections
                print "\tCorrecting raw images and calculating the angular averages for %s..." % d
                tic = time.time()
                # background subtraction
                runReferenceMean[runIndex] -= runBackgroundMean[runIndex]
                # calculate mask for assembled detector image
                mask = np.where(runReferenceMean[runIndex] > 0)
                masked = np.where(runReferenceMean[runIndex] == 0)
                # photon conversion, get photon energy from run_info
                pe = r + '/run_info/sacla_config/photon_energy_in_eV'
                if f.get(pe):
                    Ephoton = f[pe].value
                    wavelength = hPlanck*c/Ephoton
                    print "\t\tFound photon energy of %.0f eV (%.2f A) for %s" % (Ephoton, wavelength, r)
                else:
                    # used default photon energy
                    Ephoton = defaultEphoton
                    wavelength = hPlanck*c/Ephoton
                    print "\t\tUsing default photon energy of %.0f eV (%.2f A)" % (Ephoton, wavelength)
                gain = photonConvert(runReferenceMean[runIndex], Ephoton=Ephoton, Gsys=systemGain, eps=3.65)
                print "\t\tDetermined photon gain to be %.1f ADU/photons at %.0f eV" % (gain, Ephoton)
                # calculate angular averages
                q_ref, I_ref = radialProfile(runReferenceMean[runIndex], (xc, yc), mask=mask, wavelength=defaultWavelength, detectorDistance=defaultDetectorDistance, pixelSize=pixelSize)
                # enable garbage collector to clean up memory and avoid 'High memory usage errors'
                gc.enable()
                if scanRun:
                    q_data = []
                    I_data = []
                    I_div = []
                    I_diff = []
                    print "\tCalculating the angular averages for %d points in %s delay scan..." % (len(runDataMean[runIndex]), options.motor)
                    for n in range(len(runDataMean[runIndex])):
                        runDataMean[runIndex][n] -= runBackgroundMean[runIndex]
                        photonConvert(runDataMean[runIndex][n], Ephoton=Ephoton, Gsys=systemGain, eps=3.65)
                        q, Iq = radialProfile(runDataMean[runIndex][n], (xc, yc), mask=mask, wavelength=defaultWavelength, detectorDistance=defaultDetectorDistance, pixelSize=pixelSize)
                        # collect unused memory to avoid 'High memory usage errors', only runs on the nodes...
                        gc.collect()
                        q_data.append(q)
                        I_data.append(Iq)
                        I_div.append(Iq/I_ref)
                        # normalize by maximum
                        I_diff.append(Iq/Iq.max() - I_ref/I_ref.max())
                        if options.verbose:
                            print "\t\tscan_%04d" % n
                else:
                    runDataMean[runIndex] -= runBackgroundMean[runIndex]
                    photonConvert(runDataMean[runIndex], Ephoton=Ephoton, Gsys=systemGain, eps=3.65)
                    q_data, I_data = radialProfile(runDataMean[runIndex], (xc, yc), mask=mask, wavelength=defaultWavelength, detectorDistance=defaultDetectorDistance, pixelSize=pixelSize)
                    I_div = I_data/I_ref
                    # normalize by maximum
                    I_diff = I_data/I_data.max() - I_ref/I_ref.max()
                toc = time.time()
                times.append(toc - tic)
                print "\tCorrected images and calculated angular averages in %.2f s" % (toc - tic)
                
                print "\tWriting new corrected data arrays and angular averages to '%s'..." % outputFileName
                tic = time.time()
                if readWrite:
                    if scanRun:
                        pm = p + '/' + options.motor + '_scan'
                        for n in range(len(runDataMotorValues[runIndex])):
                        	pd = pm + '/%04d/analysis' % n
                                if f.get(pd):
                                    del f[pd]
                                    if options.verbose:
                                        print "\t\tDeleting %s..." % pd
                        	pa = f.require_group(pd)
                        	pa.create_dataset("correctedImage", data=runDataMean[runIndex][n])
                        	paa = pa.require_group("angularAverage")
                        	paa.create_dataset("Q", data=q_data[n])
                        	paa.create_dataset("raw", data=I_data[n])
                        	paa.create_dataset("normalized", data=(I_data[n]/I_data[n].max()))
                        	paa.create_dataset("divided", data=I_div[n])
                        	paa.create_dataset("normalizedDifference", data=I_diff[n])
                    else:
                        pd = p + '/data/analysis'
                        if f.get(pd):
                            del f[pd]
                            if options.verbose:
                                print "\t\tDeleting %s..." % pd
                        pa = f.require_group(pd)
                        pa.create_dataset("correctedImage", data=runDataMean[runIndex])
                        paa = pa.require_group("angularAverage")
                        paa.create_dataset("Q", data=q_data)
                        paa.create_dataset("raw", data=I_data)
                        paa.create_dataset("normalized", data=(I_data/I_data.max()))
                        paa.create_dataset("divided", data=I_div)
                        paa.create_dataset("normalizedDifference", data=I_diff)
                pd = p + '/reference/analysis'
                if f.get(pd):
                    del f[pd]
                    if options.verbose:
                        print "\t\tDeleting %s..." % pd
                pa = f.require_group(pd)
                pa.create_dataset("correctedImage", data=runReferenceMean[runIndex])
                paa = pa.require_group("angularAverage")
                paa.create_dataset("Q", data=q_ref)
                paa.create_dataset("raw", data=I_ref)
                paa.create_dataset("normalized", data=(I_ref/I_ref.max()))
                toc = time.time()
                times.append(toc - tic)
                print "\tSaved data in %.2f s" % (toc - tic)
            detectorIndex += 1
        runIndex += 1
    f.close()
    totalTime = np.sum(times)
    if (totalTime < 60):
        print "Successfully analyzed '%s' in %.2f s" % (options.inputFile, totalTime)
    elif (totalTime < 3600):
        print "Successfully analyzed '%s' in %d min %.2f s" % (options.inputFile, totalTime/60, totalTime%60)
    else:
        print "Successfully analyzed '%s' in %d h %d min %.2f s" % (options.inputFile, totalTime/3600, (totalTime%3600)/60, totalTime%60)
else:
    if options.inputFile != '':
        print "Input file '%s' does not exist, aborting..." % options.inputFile
    else:
        print "No input file specified, aborting..."
    sys.exit(1)
