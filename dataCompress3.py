#!/usr/bin/env python

import numpy as np
import h5py as h
import glob as g
import matplotlib
import matplotlib.pyplot as plt
import csv, sys, os, re, shutil, subprocess, time
from OrderedSet import OrderedSet
from optparse import OptionParser

"""
# Example:
#    ./dataCompress3.py -f 259782.h5 -c run259782_out.csv -a
# For details, type:
#    ./dataCompress3.py -h
"""

parser = OptionParser()
parser.add_option("-f", "--inputFile", action="store", type="string", dest="inputFile", help="Input HDF5 file produced by DataConvert3", metavar="FILENAME", default="")
parser.add_option("-c", "--csvFile", action="store", type="string", dest="csvFile", help="Input CSV spreadsheet from the scan output (default: XXXXXX.csv, where XXXXXX.h5 is the input HDF5 file)", metavar="FILENAME", default="")
parser.add_option("-o", "--outputFile", action="store", type="string", dest="outputFile", help="Output HDF5/CXIDB file that the compressed data is written to (default: XXXXXX_compressed.h5/cxi, where XXXXXX.h5 is the input HDF5 file)", metavar="FILENAME", default="")
parser.add_option("-d", "--detector", action="store", type="string", dest="detector", help="Detector name that should be compressed (default: all)", metavar="DETECTORNAME", default="")
parser.add_option("-r", "--run", action="store", type="int", dest="runNumber", help="Run number to compress (default: all in input HDF5 file)", metavar="RUNNUMBER", default=0)
parser.add_option("-n", "--nshots", action="store", type="int", dest="nShots", help="Number of shots to compress (default: all in input HDF5 file)", metavar="SHOTS", default=0)
parser.add_option("-a", "--average", action="store_true", dest="average", help="Save averages to output HDF5 file", default=False)
parser.add_option("-x", "--cxidb", action="store_true", dest="cxidb", help="Save averages and individual shots to output CXIDB file", default=False)
parser.add_option("-t", "--temperature", action="store_true", dest="temperature", help="Save temperature data to output HDF5 file", default=False)
parser.add_option("-u", "--update", action="store_true", dest="update", help="Update csv directory using rsync, only works from the computing nodes!", default=False)
parser.add_option("-v", "--verbose", action="store_true", dest="verbose", help="Output additional information to screen", default=False)
#parser.add_option("-e", "--entry", action="store", type="int", dest="entryNumber", help="Entry number in CXIDB file you wish to save model/tags for (default = 1)", metavar="NUMBER", default=1)
#parser.add_option("-i", "--image", action="store", type="int", dest="imageNumber", help="Image number in CXIDB file you wish to save model/tags for (default = 1)", metavar="NUMBER", default=1)
#parser.add_option("-d", "--dat", action="store", type="string", dest="dataVersion", help="Data version inside the CXIDB file you wish to save model/tags for (default: detector_and_photon_corrected)", metavar="VERSION", default="detector_and_photon_corrected")

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

# update csv directory
os.system("rsync --ignore-existing fep:/xdaq/work/share/fperakis/scan_output/*.csv /UserData/fperakis/csv/scan_output/.")

# PARAMETERS
tagIdentifier = 'Tagnumber'
motorsToIgnore = ['PM_EH2_2']

def removeStandardParameters(parameterList):
    """
    removes the parameters from parameterList that always are not motors, listed below:
    
    'SCAN_ID', 'Run number', 'Tagnumber', 'XFEL shutter', 'Laser shutter'
    
    returns a list that only inludes the motors
    """
    standardParameters = ['SCAN_ID', 'Run number', 'Tagnumber', 'XFEL shutter', 'Laser shutter'] + motorsToIgnore
    standardParameterSet = set(standardParameters)
    standardParameterSet.update(motorsToIgnore)
    output = []
    for p in parameterList:
        standardParameterSet.add(p)
        if (len(standardParameters) != len(standardParameterSet)):
            output.append(p)
            standardParameterSet.remove(p)
    return output


if options.inputFile != '' and os.path.exists(source_dir + options.inputFile):
    i = re.split('\W+', options.inputFile)
    if (len(i) > 2 and (options.csvFile == '' or options.outputFile == '')):
        print "Input filename '%s' does not match run number to automatically produce CSV or output filename, aborting..."
        sys.exit(1)
    else:
        basename = i[0]
    if options.csvFile != '':
        csvFileName = options.csvFile
    else:
        csvFileName = basename + '.csv'
    if os.path.exists(csv_dir + csvFileName):
        times = []
        print "Reading scan parameters from '%s'..." % (csvFileName)
        tic = time.time()
        with open(csv_dir + csvFileName, 'rb') as csvFile:
            csvReader = csv.reader(csvFile)
            csvList = list(csvReader)
            csvArray = np.array(csvList[1:], dtype=np.float)
            csvHeaders = np.array(csvList[0])
            csvDict = {}
            i = 0
            for m in csvHeaders:
                csvDict[m] = csvArray[:, i]
                i += 1
            toc = time.time()
            times.append(toc - tic)
            print "\tFound %d scan parameters in %.2f s:" % (len(csvHeaders), toc - tic)
            print "\t\t" + ', '.join(csvHeaders)
            tic = time.time()
            # check difference between rows to determine scan motor
            # don't use end of array if the motor is returned to initial position
            #(csvArray[-1] - csvArray[0]) != 0
            parametersThatChanged = csvHeaders[(csvArray[len(csvArray)/2] - csvArray[0]) != 0]
            print "\tFound %d scan parameters that changed:" % len(parametersThatChanged)
            print "\t\t" + ', '.join(parametersThatChanged)
            motorsThatMoved = removeStandardParameters(parametersThatChanged)
            print "\tFound %d motor(s) that moved:" % len(motorsThatMoved)
            print "\t\t" + ', '.join(motorsThatMoved)
            motorTags = csvDict[tagIdentifier].astype(np.int64)
            if len(motorsThatMoved) == 1:
                averageAllHits = False
                motorValues = csvDict[motorsThatMoved[0]]
                # unordered set
                #motorValuesSet = set(motorValues)
                # ordered set (original order)
                motorValuesOrderedSet = OrderedSet(motorValues)
                motorValuesOrderedArray = []
                for m in motorValuesOrderedSet:
                    motorValuesOrderedArray.append(m)
                motorValuesOrderedArray = np.array(motorValuesOrderedArray)
            elif (len(motorsThatMoved > 1)):
                print "Several motors moved at once during run, aborting..."
                sys.exit(1)
            else:
                print "No motor moved during run, all hits will be averaged."
                averageAllHits = True
            xraysStatus = csvDict['XFEL shutter']
            laserStatus = csvDict['Laser shutter']
    elif options.average:
        averageAllHits = True
        times = []
        tic = time.time()
    else:
        print "CSV file '%s' does not exist, try run with -a to average all hits, aborting..." % csvFileName
        sys.exit(1)
    
    print "Reading data from '%s'..." % (options.inputFile)
    # open input h5 file
    f = h.File(source_dir + options.inputFile, "r")
    # open output file
    if options.outputFile == '':
        outputFileName = basename + '_compressed.h5'
        if options.cxidb:
            print "CXIDB not yet implemented, aborting..."
            sys.exit(1)
    else:
        outputFileName = options.outputFile
    #outputFileExists = os.path.exists(options.outputFile)
    #o = h.File(working_dir + outputFileName, "a")
    outputFileExists = False
    if outputFileExists:
        # NOT USED, FIX IF WE WANT TO USE
        print "Output file '%s' already exists, appending data to file..." % outputFileName
        if options.model and o.get(outputModelPath):
            print "\tremoving '/%s' in %s..." % (outputModelPath, options.outputFile)
            del o[outputModelPath]
        if options.tags and o.get(outputTagsPath):
            print "\tremoving '/%s' in %s..." % (outputTagsPath, options.outputFile)
            del o[outputTagsPath]
    else:
        print "Creating file structure in '%s'..." % outputFileName
    o = h.File(working_dir + outputFileName, "w")
    
    runList = np.array(f['file_info/run_number_list'])
    print "\tFound %d run(s) in file:" % len(runList)
    print "\t\t" + ', '.join(runList.astype('|S10'))
    if options.runNumber:
        if (runList == options.runNumber).any():
            # reduce runList to only the run chosen
            runList = np.array([options.runNumber])
        else:
            print "Could not find run %d in '%s', aborting..." % (options.runNumber, options.inputFile)
            sys.exit(1)
    toc = time.time()
    times.append(toc - tic)
    runIndex = 0
    #averageAllHits = True
    for r in runList:
        print "Compressing shots for run %d..." % r
        inputPath = "run_%d" % r
        runKeys = np.array(f[inputPath].keys())
        outputRunGroup = o.require_group(inputPath)
        if options.detector != '':
             if (runKeys == options.detector).any():
                 # reduce runKeys to only the detector chosen
                 runKeys = np.array([options.detector] + [d for d in runKeys if "info" in d])
             else:
                 print "\tDetector '%s' does not exist in '%s', aborting..." % (options.detector, inputPath)
                 sys.exit(1)
        if averageAllHits:
            xraysStatus = np.array(f[inputPath + '/event_info/bl_3/eh_1/xfel_pulse_selector_status'])
            laserStatus = np.array(f[inputPath + '/event_info/bl_3/lh_1/laser_pulse_selector_status'])
            beamStatus = np.array(f[inputPath + '/event_info/acc/accelerator_status'])
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
        detectorIndex = 0
        for d in runKeys:
            p = inputPath + '/' + d
            if "info" in d:
                print "\tCopying %s for %s to '%s'..." % (d, r, outputFileName)
                # copy whole detector info group to new file
                f.copy(p, outputRunGroup)
            if "detector" in d:
                detectorKeys = np.array(f[p].keys())
                outputDetectorGroup = outputRunGroup.require_group(d)
                #outputDetectorGroup = o.create_group(p)
        	# data
        	detectorDataSum = np.array([])
        	detectorDataSquaredSum = np.array([])
        	detectorDataExcludedSum = np.array([])
        	detectorDataExcludedSquaredSum = np.array([])
        	# background
        	detectorBackgroundSum = np.array([])
        	detectorBackgroundSquaredSum = np.array([])
        	detectorBackgroundExcludedSum = np.array([])
        	detectorBackgroundExcludedSquaredSum = np.array([])
        	# reference
        	detectorReferenceSum = np.array([])
        	detectorReferenceSquaredSum = np.array([])
        	detectorReferenceExcludedSum = np.array([])
        	detectorReferenceExcludedSquaredSum = np.array([])
                if options.temperature:
        		# data temperature
        		detectorTemperatureSum = np.array([])
        		detectorTemperatureSquaredSum = np.array([])
        		detectorTemperatureExcludedSum = np.array([])
        		detectorTemperatureExcludedSquaredSum = np.array([])
        		# background temperature
        		detectorBGTemperatureSum = np.array([])
        		detectorBGTemperatureSquaredSum = np.array([])
        		detectorBGTemperatureExcludedSum = np.array([])
        		detectorBGTemperatureExcludedSquaredSum = np.array([])
        		# reference temperature
        		detectorRefTemperatureSum = np.array([])
        		detectorRefTemperatureSquaredSum = np.array([])
        		detectorRefTemperatureExcludedSum = np.array([])
        		detectorRefTemperatureExcludedSquaredSum = np.array([])
                if averageAllHits:
                    # tags with one index (first element is included tags, second element is excluded tags)
                    detectorDataTags = [[], []]
                else:
                    # tags where first index is the scan motor position, second index is the included/excluded tags
                    detectorDataTags = [[[], []] for i in range(len(motorValuesOrderedArray))]
        	detectorReferenceTags = [[], []]
        	detectorBackgroundTags = [[], []]
                anyDataHits = False
                anyExcludedDataHits = False
                anyReferenceHits = False
                anyExcludedReferenceHits = False
                anyBackgroundShots = False
                anyExcludedBackgroundShots = False
                tags = []
                tagCounter = 0
                for t in detectorKeys:
                    pt = p + '/' + t
                    if "detector_info" in t:
                        print "\tCopying detector info for %s to '%s'..." % (d, outputFileName)
                        # copy whole detector info group to new file
                        f.copy(pt, outputDetectorGroup)
                        #print "\tReading detector data for %d tag(s)..." % (len(detectorKeys) - 1)
                        print "\tReading %s data..." % d
                        tic = time.time()
                    elif "tag" in t:
                        # np.int might be too small for ~9 digit integer (32-bit has ~10 digits precision), going for 64-bit
                        tagNumber = np.int64((re.sub('tag_', '', t)))
                        tags.append(tagNumber)
                        if averageAllHits:
                            tagIndex = tagCounter
                        else:
                            tagIndex = np.where(motorTags == tagNumber)
                        detectorStatus = f[pt + '/detector_status'].value
                        # save data and temperature as 64-bit for sum and squared sum to avoid rounding errors when calculating stdev
                        detectorData = np.array(f[pt + '/detector_data'], dtype=np.float64)
                        if options.temperature:
                            detectorTemperature = np.float64(f[pt + '/temperature'].value)
                        # xrays on and laser on
                        if ((detectorStatus == 1) and (xraysStatus[tagIndex] == 1) and (laserStatus[tagIndex] == 1)):
                            anyDataHits = True
                            if averageAllHits:
                                detectorDataTags[0].append(tagNumber)
                            	if detectorDataSum.any():
                            	    detectorDataSum += detectorData
                            	    detectorDataSquaredSum += detectorData*detectorData
                                    if options.temperature:
                                        detectorTemperatureSum += detectorTemperature
                                        detectorTemperatureSquaredSum += detectorTemperature*detectorTemperature
                            	else:
                                    # create regular 2D array for sum if all hits are averaged
                            	    detectorDataSum = detectorData
                            	    detectorDataSquaredSum = detectorData*detectorData
                                    if options.temperature:
                                        detectorTemperatureSum = detectorTemperature
                                        detectorTemperatureSquaredSum = detectorTemperature*detectorTemperature
                                if options.verbose:
                                    print "\t\t[%d] hit" % tagNumber
                            else:
                            	if not detectorDataSum.any():
                                    # create 3D array for sum with first index being the scan motor position
                            	    detectorDataSum = np.zeros((len(motorValuesOrderedArray),) + detectorData.shape, dtype=np.float64)
                            	    detectorDataSquaredSum = np.zeros((len(motorValuesOrderedArray),) + detectorData.shape, dtype=np.float64)
                                    if options.temperature:
                                        detectorTemperatureSum = np.zeros(len(motorValuesOrderedArray), dtype=np.float64)
                                        detectorTemperatureSquaredSum = np.zeros(len(motorValuesOrderedArray), dtype=np.float64)
                                # get rid of wrapping tuple so index works in list as well as numpy array
                                motorValuesOrderedIndex = np.where(motorValuesOrderedArray == motorValues[tagIndex])[0]
                                detectorDataTags[motorValuesOrderedIndex][0].append(tagNumber)
                                detectorDataSum[motorValuesOrderedIndex] += detectorData
                                detectorDataSquaredSum[motorValuesOrderedIndex] += detectorData*detectorData
                                if options.temperature:
                                    detectorTemperatureSum[motorValuesOrderedIndex] += detectorTemperature
                                    detectorTemperatureSquaredSum[motorValuesOrderedIndex] += detectorTemperature*detectorTemperature
                                if options.verbose:
                                    print "\t\t[%d] hit (%s = %.2e)" % (tagNumber, motorsThatMoved[0], motorValues[tagIndex])
                        # xrays on and laser on (excluded)
                        elif ((detectorStatus == 0) and (xraysStatus[tagIndex] == 1) and (laserStatus[tagIndex] == 1)):
                            anyExcludedDataHits = True
                            if averageAllHits:
                                detectorDataTags[1].append(tagNumber)
                            	if detectorDataExcludedSum.any():
                            	    detectorDataExcludedSum += detectorData
                            	    detectorDataExcludedSquaredSum += detectorData*detectorData
                                    if options.temperature:
                                        detectorTemperatureExcludedSum += detectorTemperature
                                        detectorTemperatureExcludedSquaredSum += detectorTemperature*detectorTemperature
                            	else:
                            	    detectorDataExcludedSum = detectorData
                            	    detectorDataExcludedSquaredSum = detectorData*detectorData
                                    if options.temperature:
                                        detectorTemperatureExcludedSum = detectorTemperature
                                        detectorTemperatureExcludedSquaredSum = detectorTemperature*detectorTemperature
                                if options.verbose:
                                    print "\t\t[%d] excluded hit" % tagNumber
                            else:
                            	if not detectorDataExcludedSum.any():
                                    # create 3D array for sum with first index being the scan motor position
                            	    detectorDataExcludedSum = np.zeros((len(motorValuesOrderedArray),) + detectorData.shape, dtype=np.float64)
                            	    detectorDataExcludedSquaredSum = np.zeros((len(motorValuesOrderedArray),) + detectorData.shape, dtype=np.float64)
                                    if options.temperature:
                                        detectorTemperatureExcludedSum = np.zeros(len(motorValuesOrderedArray), dtype=np.float64)
                                        detectorTemperatureExcludedSquaredSum = np.zeros(len(motorValuesOrderedArray), dtype=np.float64)
                                motorValuesOrderedIndex = np.where(motorValuesOrderedArray == motorValues[tagIndex])
                                print motorValuesOrderedIndex[0]
                                detectorDataTags[motorValuesOrderedIndex][1].append(tagNumber)
                                detectorDataExcludedSum[motorValuesOrderedIndex] += detectorData
                                detectorDataExcludedSquaredSum[motorValuesOrderedIndex] += detectorData*detectorData
                                if options.temperature:
                                    detectorTemperatureExcludedSum[motorValuesOrderedIndex] += detectorTemperature
                                    detectorTemperatureExcludedSquaredSum[motorValuesOrderedIndex] += detectorTemperature*detectorTemperature
                                if options.verbose:
                                    print "\t\t[%d] excluded hit (%s = %.2e)" % (tagNumber, motorsThatMoved[0], motorValues[tagIndex])
                        # xrays on and laser off
                        elif ((detectorStatus == 1) and (xraysStatus[tagIndex] == 1) and (laserStatus[tagIndex] == 0)):
                            anyReferenceHits = True
                            detectorReferenceTags[0].append(tagNumber)
                            if detectorReferenceSum.any():
                                detectorReferenceSum += detectorData
                                detectorReferenceSquaredSum += detectorData*detectorData
                                if options.temperature:
                                    detectorRefTemperatureSum += detectorTemperature
                                    detectorRefTemperatureSquaredSum += detectorTemperature*detectorTemperature
                            else:
                                detectorReferenceSum = detectorData
                                detectorReferenceSquaredSum = detectorData*detectorData
                                if options.temperature:
                                    detectorRefTemperatureSum = detectorTemperature
                                    detectorRefTemperatureSquaredSum = detectorTemperature*detectorTemperature
                            if options.verbose:
                                print "\t\t[%d] reference hit" % tagNumber
                        # xrays on and laser off (excluded)
                        elif ((detectorStatus == 0) and (xraysStatus[tagIndex] == 1) and (laserStatus[tagIndex] == 0)):
                            anyExcludedReferenceHits = True
                            detectorReferenceTags[1].append(tagNumber)
                            if detectorReferenceExcludedSum.any():
                                detectorReferenceExcludedSum += detectorData
                                detectorReferenceExcludedSquaredSum += detectorData*detectorData
                                if options.temperature:
                                    detectorRefTemperatureExcludedSum += detectorTemperature
                                    detectorRefTemperatureExcludedSquaredSum += detectorTemperature*detectorTemperature
                            else:
                                detectorReferenceExcludedSum = detectorData
                                detectorReferenceExcludedSquaredSum = detectorData*detectorData
                                if options.temperature:
                                    detectorRefTemperatureExcludedSum = detectorTemperature
                                    detectorRefTemperatureExcludedSquaredSum = detectorTemperature*detectorTemperature
                            if options.verbose:
                                print "\t\t[%d] excluded reference hit" % tagNumber
                        # xrays off
                        elif ((detectorStatus == 1) and (xraysStatus[tagIndex] == 0)):
                            anyBackgroundShots = True
                            detectorBackgroundTags[0].append(tagNumber)
                            if detectorBackgroundSum.any():
                                detectorBackgroundSum += detectorData
                                detectorBackgroundSquaredSum += detectorData*detectorData
                                if options.temperature:
                                    detectorBGTemperatureSum += detectorTemperature
                                    detectorBGTemperatureSquaredSum += detectorTemperature*detectorTemperature
                            else:
                                detectorBackgroundSum = detectorData
                                detectorBackgroundSquaredSum = detectorData*detectorData
                                if options.temperature:
                                    detectorBGTemperatureSum = detectorTemperature
                                    detectorBGTemperatureSquaredSum = detectorTemperature*detectorTemperature                                
                            if options.verbose:
                                print "\t\t[%d] background shot" % tagNumber
                        # xrays off (excluded)
                        elif ((detectorStatus == 0) and (xraysStatus[tagIndex] == 0)):
                            anyExcludedBackgroundShots = True
                            detectorBackgroundTags[1].append(tagNumber)
                            if detectorBackgroundExcludedSum.any():
                                detectorBackgroundExcludedSum += detectorData
                                detectorBackgroundExcludedSquaredSum += detectorData*detectorData
                                if options.temperature:
                                    detectorBGTemperatureExcludedSum += detectorTemperature
                                    detectorBGTemperatureExcludedSquaredSum += detectorTemperature*detectorTemperature
                            else:
                                detectorBackgroundExcludedSum = detectorData
                                detectorBackgroundExcludedSquaredSum = detectorData*detectorData
                                if options.temperature:
                                    detectorBGTemperatureExcludedSum = detectorTemperature
                                    detectorBGTemperatureExcludedSquaredSum = detectorTemperature*detectorTemperature
                            if options.verbose:
                                print "\t\t[%d] excluded background shot" % tagNumber
                        else:
                            if xraysStatus[tagIndex] and laserStatus[tagIndex]:
                                print "Tag %d could not be classified (detector = %d, xrays = %d, laser = %d)..." % (tagNumber, detectorStatus, xraysStatus[tagIndex], laserStatus[tagIndex])
                            else:
                                print "Tag %d could not be classified (detector = %d)..." % (tagNumber, detectorStatus)
                        tagCounter += 1
                    # this is used to break the script for testing
                    if (options.nShots > 0) and (tagCounter >= options.nShots):
                        break
                toc = time.time()
                times.append(toc - tic)
                print "\tRead detector data for %d tags in %.2f s" % (tagCounter, toc - tic)
                
                # print statistics
                hitStatsString = ""
                if anyDataHits:
                    if averageAllHits:
                        hitStatsString += "\t\t%d hits\n" % len(detectorDataTags[0])
                    else:
                        totalHits = 0
                        for m in range(len(motorValuesOrderedArray)):
                            if detectorDataTags[m][0]:
                                hitStatsString += "\t\t%d hits [%s = %.2f]\n" % (len(detectorDataTags[m][0]), motorsThatMoved[0], motorValuesOrderedArray[m])
                                totalHits += len(detectorDataTags[m][0])
                        hitStatsString += "\t\t%d total hits\n" % totalHits
                if anyExcludedDataHits:
                    if averageAllHits:
                        hitStatsString += "\t\t(%d excluded hits)\n" % len(detectorDataTags[0])
                    else:
                        totalHits = 0
                        for m in range(len(motorValuesOrderedArray)):
                            if detectorDataTags[m][1]:
                                hitStatsString += "\t\t(%d excludedhits [%s = %.2f])\n" % (len(detectorDataTags[m][1]), motorsThatMoved[0], motorValuesOrderedArray[m])
                                totalHits += len(detectorDataTags[m][1])
                        hitStatsString += "\t\t(%d total hits)\n" % totalHits
                if anyReferenceHits:
                    hitStatsString += "\t\t%d reference hits\n" % len(detectorReferenceTags[0])
                if anyExcludedReferenceHits:
                    hitStatsString += "\t\t(%d excluded reference hits)\n" % len(detectorReferenceTags[1])
                if anyBackgroundShots:
                    hitStatsString += "\t\t%d background hits\n" % len(detectorBackgroundTags[0])
                if anyExcludedBackgroundShots:
                    hitStatsString += "\t\t(%d excluded background hits)\n" % len(detectorBackgroundTags[1])
                if hitStatsString != '':
                    print hitStatsString[:-1]
                tags = np.array(tags)
                
                # calculate mean and standard deviation
                print "\tCalculating data statistics..."
                tic = time.time()
                if anyDataHits:
                    if averageAllHits:
                        # calculate mean
                        runDataMean.append(detectorDataSum/np.float(len(detectorDataTags[0])))
                        # calculate stdev from sqrt(E[X^2] - E[X]^2), only works for 1/N (non-corrected) normalization
                        runDataStd.append(np.sqrt(detectorDataSquaredSum - detectorDataSum*detectorDataSum/np.float(len(detectorDataTags[0])))/np.float(len(detectorDataTags[0])))
                        if options.temperature:
                            runTemperatureMean.append(detectorTemperatureSum/np.float(len(detectorDataTags[0])))
                            # round root-operand to 10 decimals to avoid negative numbers for really small deviations
                            runTemperatureStd.append(np.sqrt(round((detectorTemperatureSquaredSum - detectorTemperatureSum*detectorTemperatureSum/np.float(len(detectorDataTags[0])))/np.float(len(detectorDataTags[0])))))
                    else:
                        # initialize arrays to zero
                        runDataMean.append(np.zeros_like(detectorDataSum))
                        runDataStd.append(np.zeros_like(detectorDataSum))
                        if options.temperature:
                            runTemperatureMean.append(np.zeros_like(detectorTemperatureSum))
                            runTemperatureStd.append(np.zeros_like(detectorTemperatureSum))
                        for m in range(len(motorValuesOrderedArray)):
                            if detectorDataTags[m][0]:
                                # calculate mean
                                runDataMean[runIndex][m] = detectorDataSum[m]/np.float(len(detectorDataTags[m][0]))
                                # calculate stdev from sqrt(E[X^2] - E[X]^2), only works for 1/N (non-corrected) normalization
                                runDataStd[runIndex][m] = np.sqrt((detectorDataSquaredSum[m] - detectorDataSum[m]*detectorDataSum[m]/np.float(len(detectorDataTags[m][0])))/np.float(len(detectorDataTags[m][0])))
                                if options.temperature:
                                    runTemperatureMean[runIndex][m] = detectorTemperatureSum[m]/np.float(len(detectorDataTags[m][0]))
                                    # round root-operand to 10 decimals to avoid negative numbers for really small deviations
                                    runTemperatureStd[runIndex][m] = np.sqrt(round((detectorTemperatureSquaredSum[m] - detectorTemperatureSum[m]*detectorTemperatureSum[m]/np.float(len(detectorDataTags[m][0])))/np.float(len(detectorDataTags[m][0]))*1.0E10)/1.0E10)
                if anyExcludedDataHits:
                    if averageAllHits:
                        # calculate mean
                        runDataExcludedMean.append(detectorDataExcludedSum/np.float(len(detectorDataTags[1])))
                        # calculate stdev from sqrt(E[X^2] - E[X]^2), only works for 1/N (non-corrected) normalization
                        runDataExcludedStd.append(np.sqrt((detectorDataExcludedSquaredSum - detectorDataExcludedSum*detectorDataExcludedSum/np.float(len(detectorDataTags[1])))/np.float(len(detectorDataTags[1]))))
                        if options.temperature:
                            runTemperatureExcludedMean.append(detectorTemperatureExcludedSum/np.float(len(detectorDataTags[0])))
                            # round root-operand to 10 decimals to avoid negative numbers for really small deviations
                            runTemperatureExcludedStd.append(np.sqrt(round((detectorTemperatureExcludedSquaredSum - detectorTemperatureExcludedSum*detectorTemperatureExcludedSum/np.float(len(detectorDataTags[1])))/np.float(len(detectorDataTags[0]))*1.0E10)/1.0E10))
                    else:
                        # initialize arrays to zero
                        runDataExcludedMean.append(np.zeros_like(detectorDataExcludedSum))
                        runDataExcludedStd.append(np.zeros_like(detectorDataExcludedSum))
                        if options.temperature:
                            runTemperatureExcludedMean.append(np.zeros_like(detectorTemperatureExcludedSum))
                            runTemperatureExcludedStd.append(np.zeros_like(detectorTemperatureExcludedSum))
                        for m in range(len(motorValuesOrderedArray)):
                            if detectorDataTags[m][1]:
                                # calculate mean
                                runDataExcludedMean[runIndex][m] = detectorDataExcludedSum[m]/np.float(len(detectorDataTags[m][1]))
                                # calculate stdev from sqrt(E[X^2] - E[X]^2), only works for 1/N (non-corrected) normalization
                                runDataExcludedStd[runIndex][m] = np.sqrt((detectorDataExcludedSquaredSum[m] - detectorDataExcludedSum[m]*detectorDataExcludedSum[m]/np.float(len(detectorDataTags[m][1])))/np.float(len(detectorDataTags[m][1])))
                                if options.temperature:
                                    runTemperatureExcludedMean[runIndex][m] = detectorTemperatureExcludedSum[m]/np.float(len(detectorDataTags[m][1]))
                                    # round root-operand to 10 decimals to avoid negative numbers for really small deviations
                                    runTemperatureExcludedStd[runIndex][m] = np.sqrt(round((detectorTemperatureExcludedSquaredSum[m] - detectorTemperatureExcludedSum[m]*detectorTemperatureExcludedSum[m]/np.float(len(detectorDataTags[m][1])))/np.float(len(detectorDataTags[m][1]))*1.0E10)/1.0E10)
                if anyReferenceHits:
                        # calculate mean
                        runReferenceMean.append(detectorReferenceSum/np.float(len(detectorReferenceTags[0])))
                        # calculate stdev from sqrt(E[X^2] - E[X]^2), only works for 1/N (non-corrected) normalization
                        runReferenceStd.append(np.sqrt((detectorReferenceSquaredSum - detectorReferenceSum*detectorReferenceSum/np.float(len(detectorReferenceTags[0])))/np.float(len(detectorReferenceTags[0]))))
                        if options.temperature:
                            runRefTemperatureMean.append(detectorRefTemperatureSum/np.float(len(detectorReferenceTags[0])))
                            # round root-operand to 10 decimals to avoid negative numbers for really small deviations
                            runRefTemperatureStd.append(np.sqrt(round((detectorRefTemperatureSquaredSum - detectorRefTemperatureSum*detectorRefTemperatureSum/np.float(len(detectorReferenceTags[0])))/np.float(len(detectorReferenceTags[0]))*1.0E10)/1.0E10))
                if anyExcludedReferenceHits:
                        # calculate mean
                        runReferenceExcludedMean.append(detectorReferenceExcludedSum/np.float(len(detectorReferenceTags[1])))
                        # calculate stdev from sqrt(E[X^2] - E[X]^2), only works for 1/N (non-corrected) normalization
                        runReferenceExcludedStd.append(np.sqrt(round((detectorReferenceExcludedSquaredSum - detectorReferenceExcludedSum*detectorReferenceExcludedSum/np.float(len(detectorReferenceTags[1])))/np.float(len(detectorReferenceTags[1]))*1.0E10)/1.0E10))
                        if options.temperature:
                            runRefTemperatureExcludedMean.append(detectorRefTemperatureExcludedSum/np.float(len(detectorReferenceTags[0])))
                            # round root-operand to 10 decimals to avoid negative numbers for really small deviations
                            runRefTemperatureExcludedStd.append(np.sqrt((detectorRefTemperatureExcludedSquaredSum - detectorRefTemperatureExcludedSum*detectorRefTemperatureExcludedSum/np.float(len(detectorReferenceTags[1])))/np.float(len(detectorReferenceTags[0]))))
                if anyBackgroundShots:
                        # calculate mean
                        runBackgroundMean.append(detectorBackgroundSum/np.float(len(detectorBackgroundTags[0])))
                        # calculate stdev from sqrt(E[X^2] - E[X]^2), only works for 1/N (non-corrected) normalization
                        runBackgroundStd.append(np.sqrt((detectorBackgroundSquaredSum - detectorBackgroundSum*detectorBackgroundSum/np.float(len(detectorBackgroundTags[0])))/np.float(len(detectorBackgroundTags[0]))))
                        if options.temperature:
                            runBGTemperatureMean.append(detectorBGTemperatureSum/np.float(len(detectorBackgroundTags[0])))
                            # round root-operand to 10 decimals to avoid negative numbers for really small deviations
                            runBGTemperatureStd.append(np.sqrt(round((detectorBGTemperatureSquaredSum - detectorBGTemperatureSum*detectorBGTemperatureSum/np.float(len(detectorBackgroundTags[0])))/np.float(len(detectorBackgroundTags[0]))*1.0E10)/1.0E10))
                if anyExcludedBackgroundShots:
                        # calculate mean
                        runBackgroundExcludedMean.append(detectorBackgroundExcludedSum/np.float(len(detectorBackgroundTags[1])))
                        # calculate stdev from sqrt(E[X^2] - E[X]^2), only works for 1/N (non-corrected) normalization
                        runBackgroundExcludedStd.append(np.sqrt((detectorBackgroundExcludedSquaredSum - detectorBackgroundExcludedSum*detectorBackgroundExcludedSum/np.float(len(detectorBackgroundTags[1])))/np.float(len(detectorBackgroundTags[1]))))
                        if options.temperature:
                            runBGTemperatureExcludedMean.append(detectorBGTemperatureExcludedSum/np.float(len(detectorBackgroundTags[0])))
                            # round root-operand to 10 decimals to avoid negative numbers for really small deviations
                            runBGTemperatureExcludedStd.append(np.sqrt(round((detectorBGTemperatureExcludedSquaredSum - detectorBGTemperatureExcludedSum*detectorBGTemperatureExcludedSum/np.float(len(detectorBackgroundTags[1])))/np.float(len(detectorBackgroundTags[0]))*1.0E10)/1.0E10))
                toc = time.time()
                times.append(toc - tic)
                print "\tCalculated data statistics in %.2f s" % (toc - tic)

                # save data to output file
                print "\tSaving averaged data to '%s'..." % outputFileName
                tic = time.time()
                # xrays on and laser on
                if anyDataHits:
                    if averageAllHits:
                    	# create data groups
                    	od = outputDetectorGroup.create_group("data")
                    	odd = od.create_group("data")
                        if options.temperature:
                            odt = od.create_group("temperature")
                    	# save data to data groups
                    	odd.create_dataset("mean", data=runDataMean[runIndex])
                    	odd.create_dataset("stdev", data=runDataStd[runIndex])
                        if options.temperature:
                            odt.create_dataset("mean", data=runTemperatureMean[runIndex])
                            odt.create_dataset("stdev", data=runTemperatureStd[runIndex])
                    	od.create_dataset("tags", data=np.array(detectorDataTags[0], dtype=np.int64))
                    else:
                        om = outputDetectorGroup.create_group(motorsThatMoved[0] + "_scan")
                        excludedPositions = 0
                        for m in range(len(motorValuesOrderedArray)):
                            if detectorDataTags[m][0]:
                                # create data groups
                                omi = om.create_group("%04d" % (m - excludedPositions))
                                od = omi.create_group("data")
                                if options.temperature:
                                    ot = omi.create_group("temperature")
                                # save data to data groups
                                od.create_dataset("mean", data=runDataMean[runIndex][m])
                                od.create_dataset("stdev", data=runDataStd[runIndex][m])
                                if options.temperature:
                                    ot.create_dataset("mean", data=runTemperatureMean[runIndex][m])
                                    ot.create_dataset("stdev", data=runTemperatureStd[runIndex][m])
                                omi.create_dataset("tags", data=np.array(detectorDataTags[m][0], dtype=np.int64))
                                omi.create_dataset("motor_value", data=motorValuesOrderedArray[m])
                            elif not detectorDataTags[m][1]:
                                excludedPositions += 1
                # xrays on and laser on (excluded)
                if anyExcludedDataHits:
                    if averageAllHits:
                    	# create data groups
                    	oe = outputDetectorGroup.create_group("excluded_data")
                    	oed = oe.create_group("data")
                        if options.temperature:
                            oet = oe.create_group("temperature")
                    	# save data to data groups
                    	oed.create_dataset("mean", data=runDataExcludedMean[runIndex])
                    	oed.create_dataset("stdev", data=runDataExcludedStd[runIndex])
                        if options.temperature:
                            oet.create_dataset("mean", data=runTemperatureExcludedMean[runIndex])
                            oet.create_dataset("stdev", data=runTemperatureExcludedStd[runIndex])
                    	oe.create_dataset("tags", data=np.array(detectorDataTags[1], dtype=np.int64))
                    else:
                        om = outputDetectorGroup.require_group(motorsThatMoved[0] + "_scan")
                        excludedPositions = 0
                        for m in range(len(motorValuesOrderedArray)):
                            if detectorDataTags[m][1]:
                                # create data groups
                                omi = om.require_group("%04d" % (m - excludedPositions))
                                oe = omi.create_group("excluded_data")
                                oed = oe.create_group("data")
                                if options.temperature:
                                    oet = oe.create_group("temperature")
                                # save data to data groups
                                oed.create_dataset("mean", data=runDataExcludedMean[runIndex][m])
                                oed.create_dataset("stdev", data=runDataExcludedStd[runIndex][m])
                                if options.temperature:
                                    oet.create_dataset("mean", data=runTemperatureExcludedMean[runIndex][m])
                                    oet.create_dataset("stdev", data=runTemperatureExcludedStd[runIndex][m])
                                oe.create_dataset("tags", data=np.array(detectorDataTags[m][1], dtype=np.int64))
                                oe.create_dataset("motor_value", data=motorValuesOrderedArray[m])
                            elif not detectorDataTags[m][0]:
                                excludedPositions += 1
                # xrays on and laser off
                if anyReferenceHits:
                    # create data groups
                    orf = outputDetectorGroup.create_group("reference")
                    ord = orf.create_group("data")
                    if options.temperature:
                        ort = orf.create_group("temperature")
                    # save data to data groups
                    ord.create_dataset("mean", data=runReferenceMean[runIndex])
                    ord.create_dataset("stdev", data=runReferenceStd[runIndex])
                    if options.temperature:
                        ort.create_dataset("mean", data=runRefTemperatureMean[runIndex])
                        ort.create_dataset("stdev", data=runRefTemperatureStd[runIndex])
                    orf.create_dataset("tags", data=np.array(detectorReferenceTags[0], dtype=np.int64))
                # xrays on and laser off (excluded)
                if anyExcludedReferenceHits:
                    # create data groups
                    oer = outputDetectorGroup.create_group("excluded_reference")
                    oed = oer.create_group("data")
                    if options.temperature:
                        oet = oer.create_group("temperature")
                    # save data to data groups
                    oed.create_dataset("mean", data=runReferenceExcludedMean[runIndex])
                    oed.create_dataset("stdev", data=runReferenceExcludedStd[runIndex])
                    if options.temperature:
                        oet.create_dataset("mean", data=runRefTemperatureExcludedMean[runIndex])
                        oet.create_dataset("stdev", data=runRefTemperatureExcludedStd[runIndex])
                    oer.create_dataset("tags", data=np.array(detectorReferenceTags[1], dtype=np.int64))
                # xrays off
                if anyBackgroundShots:
                    # create data groups
                    ob = outputDetectorGroup.create_group("background")
                    obd = ob.create_group("data")
                    if options.temperature:
                        obt = ob.create_group("temperature")
                    # save data to data groups
                    obd.create_dataset("mean", data=runBackgroundMean[runIndex])
                    obd.create_dataset("stdev", data=runBackgroundStd[runIndex])
                    if options.temperature:
                        obt.create_dataset("mean", data=runBGTemperatureMean[runIndex])
                        obt.create_dataset("stdev", data=runBGTemperatureStd[runIndex])
                    ob.create_dataset("tags", data=np.array(detectorBackgroundTags[0], dtype=np.int64))
                # xrays off (excluded)
                if anyExcludedBackgroundShots:
                    # create data groups
                    oeb = outputDetectorGroup.create_group("excluded_background")
                    oed = oeb.create_group("data")
                    if options.temperature:
                        oet = oeb.create_group("temperature")
                    # save data to data groups
                    oed.create_dataset("mean", data=runBackgroundExcludedMean[runIndex])
                    oed.create_dataset("stdev", data=runBackgroundExcludedStd[runIndex])
                    if options.temperature:
                        oet.create_dataset("mean", data=runBGTemperatureExcludedMean[runIndex])
                        oet.create_dataset("stdev", data=runBGTemperatureExcludedStd[runIndex])
                    oeb.create_dataset("tags", data=np.array(detectorBackgroundTags[1], dtype=np.int64))
                toc = time.time()
                times.append(toc - tic)
                print "\tSaved data in %.2f s" % (toc - tic)
            detectorIndex += 1
        runIndex += 1
    f.close()
    o.close()
    totalTime = np.sum(times)
    if (totalTime < 60):
        print "Successfully compressed '%s' into '%s' in %.2f s" % (options.inputFile, outputFileName, totalTime)
    elif (totalTime < 3600):
        print "Successfully compressed '%s' into '%s' in %d min %.2f s" % (options.inputFile, outputFileName, totalTime/60, totalTime%60)
    else:
        print "Successfully compressed '%s' into '%s' in %d h %d min %.2f s" % (options.inputFile, outputFileName, totalTime/3600, (totalTime%3600)/60, totalTime%60)
else:
    if options.inputFile != '':
        print "Input file '%s' does not exist, aborting..." % options.inputFile
    else:
        print "No input file specified, aborting..."
    sys.exit(1)
