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
parser.add_option("-c", "--csvFile", action="store", type="string", dest="csvFile", help="Input CSV spreadsheet from the scan output (default: runXXXXXX_out.csv, where XXXXXX.h5 is the input HDF5 file)", metavar="FILENAME", default="")
parser.add_option("-o", "--outputFile", action="store", type="string", dest="outputFile", help="Output HDF5/CXIDB file that the compressed data is written to (default: XXXXXX_compressed.h5/cxi, where XXXXXX.h5 is the input HDF5 file)", metavar="FILENAME", default="")
parser.add_option("-d", "--detector", action="store", type="string", dest="detector", help="Detector name that should be compressed (default: all)", metavar="DETECTORNAME", default="")
parser.add_option("-r", "--run", action="store", type="int", dest="runNumber", help="Run number to compress (default: all in input HDF5 file)", metavar="RUNNUMBER", default=0)
parser.add_option("-a", "--average", action="store_true", dest="average", help="Save averages to output HDF5 file", default=False)
parser.add_option("-x", "--cxidb", action="store_true", dest="cxidb", help="Save averages and individual shots to output CXIDB file", default=False)
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
csv_dir = "/work/fperakis/output/"
#csv_dir = ""
source_dir = "/work/fperakis/output/"
#source_dir = ""
working_dir = "/work/fperakis/output/"
#working_dir = ""
original_dir = os.getcwd() + '/'

# PARAMETERS
tagIdentifier = 'Tagnumber'

def removeStandardParameters(parameterList):
    """
    removes the parameters from parameterList that always are not motors, listed below:
    
    'SCAN_ID', 'Run number', 'Tagnumber', 'XFEL shutter', 'Laser shutter'
    
    returns a list that only inludes the motors
    """
    standardParameters = ['SCAN_ID', 'Run number', 'Tagnumber', 'XFEL shutter', 'Laser shutter']
    standardParameterSet = set(standardParameters)
    output = []
    for p in parameterList:
        standardParameterSet.add(p)
        if (len(standardParameters) != len(standardParameterSet)):
            output.append(p)
            standardParameterSet.remove(p)
    return output


if options.inputFile != '' and os.path.exists(options.inputFile):
    i = re.split('\W+', options.inputFile)
    if (len(i) > 2 and (options.csvFile == '' or options.outputFile == '')):
        print "Input filename '%s' does not match run number to automatically produce CSV or output filename, aborting..."
        sys.exit(1)
    else:
        basename = i[0]
    if options.csvFile != '':
        csvFileName = options.csvFile
    else:
        csvFileName = 'run' + basename + '_out.csv'
    if os.path.exists(csvFileName):
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
            if len(motorsThatMoved) == 1:
                averageAllHits = False
                motorValues = csvDict[motorsThatMoved[0]]
                motorTags = csvDict[tagIdentifier]
            elif (len(motorsThatMoved > 1)):
                print "Several motors moved at once during run, aborting..."
                sys.exit(1)
            else:
                print "No motor moved during run, all hits will be averaged."
                averageAllHits = True
            xraysOn = csvDict['XFEL shutter'].astype(np.bool)
            laserOn = csvDict['Laser shutter'].astype(np.bool)
            #xraysOn = csvDict['XFEL shutter'].astype(np.bool)[:501]
            #laserOn = csvDict['Laser shutter'].astype(np.bool)[:501]
    else:
        print "CSV file '%s' does not exist, aborting..." % csvFileName
        sys.exit(1)
    
    print "Reading data from '%s'..." % (options.inputFile)
    # open input h5 file
    f = h.File(source_dir + options.inputFile, "r")
    # open output file
    if options.outputFile == '':
        if options.average:
            outputFileName = basename + '_compressed.h5'
        else:
            print "CXIDB not yet implemented, use average (-a) mode..."
            sys.exit(1)
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
    for r in runList:
        print "Compressing shots for run %d..." % r
        inputPath = "run_%d" % r
        runKeys = np.array(f[inputPath].keys())
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
        runTags = []
        runIndex = 0
        for d in runKeys:
            if "detector" in d:
                p = inputPath + '/' + d
                detectorKeys = np.array(f[p].keys())
                #outputDetectorGroup = outputRunGroup.require_group(d)
                outputDetectorGroup = o.create_group(p)
                detectorData = []
                detectorStatus = []
                detectorTemperature = []
                tags = []
                tagCounter = 0
                for t in detectorKeys:
                    pt = p + '/' + t
                    if "detector_info" in t:
                        print "\tCopying detector info for %s to '%s'..." % (d, outputFileName)
                        # copy whole detector info group to new file
                        f.copy(pt, outputDetectorGroup)
                        print "\tReading detector data for %d tag(s)..." % (len(detectorKeys) - 1)
                        tic = time.time()
                    elif "tag" in t:
                        if options.verbose:
                            print "\t\t%s" % t
                        # np.int might be too small for ~9 digit integer (32-bit has ~10 digits precision), going for 64-bit
                        tagNumber = np.int64((re.sub('tag_', '', t)))
                        tags.append(tagNumber)
                        tagCounter += 1
                        detectorData.append(np.array(f[pt + '/detector_data']))
                        detectorStatus.append(f[pt + '/detector_status'].value)
                        detectorTemperature.append(f[pt + '/temperature'].value)
                    #if tagCounter > 500:
                    #    break
                toc = time.time()
                times.append(toc - tic)
                print "\tRead detector data for %d tags in %.2f s" % (tagCounter, toc - tic)
                print "\tConverting data for '%s' to internal arrays..." % d
                tic = time.time()
                detectorData = np.array(detectorData, dtype=np.float)
                detectorStatus = np.array(detectorStatus, dtype=np.int)
                detectorTemperature = np.array(detectorTemperature, dtype=np.float)
                tags = np.array(tags)
                toc = time.time()
                times.append(toc - tic)
                print "\tConverted data in %.2f s" % (toc - tic)
                xraysOnAndLaserOn = (detectorStatus == 1) & xraysOn & laserOn
                xraysOnAndLaserOnExcluded = (detectorStatus == 0) & xraysOn & laserOn
                xraysOnAndLaserOff = (detectorStatus == 1) & xraysOn & (laserOn == 0)
                xraysOnAndLaserOffExcluded = (detectorStatus == 0) & xraysOn & (laserOn == 0)
                xraysOff = (detectorStatus == 1) & (xraysOn == 0)
                xraysOffExcluded = (detectorStatus == 0) & (xraysOn == 0)
                if averageAllHits:
                    # average data
                    print "\tAveraging all hits for '%s'..." % d
                    tic = time.time()
                    # xrays on and laser on
                    if (xraysOnAndLaserOn).any():
                        runDataMean.append(detectorData[np.where(xraysOnAndLaserOn)].mean(axis=0))
                        runDataStd.append(detectorData[np.where(xraysOnAndLaserOn)].std(axis=0))
                        runTemperatureMean.append(detectorTemperature[np.where(xraysOnAndLaserOn)].mean(axis=0))
                        runTemperatureStd.append(detectorTemperature[np.where(xraysOnAndLaserOn)].std(axis=0))
                    else:
                        print "\t\tNo laser hits included..."
                    if (xraysOnAndLaserOnExcluded).any():
                        runDataExcludedMean.append(detectorData[np.where(xraysOnAndLaserOnExcluded)].mean(axis=0))
                        runDataExcludedStd.append(detectorData[np.where(xraysOnAndLaserOnExcluded)].std(axis=0))
                        runTemperatureExcludedMean.append(detectorTemperature[np.where(xraysOnAndLaserOnExcluded)].mean(axis=0))
                        runTemperatureExcludedStd.append(detectorTemperature[np.where(xraysOnAndLaserOnExcluded)].std(axis=0))
                    else:
                        print "\t\tNo laser hits excluded..."
                    # xrays on and laser off
                    if (xraysOnAndLaserOff).any():
                        runReferenceMean.append(detectorData[np.where(xraysOnAndLaserOff)].mean(axis=0))
                        runReferenceStd.append(detectorData[np.where(xraysOnAndLaserOff)].std(axis=0))
                        runRefTemperatureMean.append(detectorTemperature[np.where(xraysOnAndLaserOff)].mean(axis=0))
                        runRefTemperatureStd.append(detectorTemperature[np.where(xraysOnAndLaserOff)].std(axis=0))
                    else:
                        print "\t\tNo reference hits included..."
                    if (xraysOnAndLaserOffExcluded).any():
                        runReferenceExcludedMean.append(detectorData[np.where(xraysOnAndLaserOffExcluded)].mean(axis=0))
                        runReferenceExcludedStd.append(detectorData[np.where(xraysOnAndLaserOffExcluded)].std(axis=0))
                        runRefTemperatureExcludedMean.append(detectorTemperature[np.where(xraysOnAndLaserOffExcluded)].mean(axis=0))
                        runRefTemperatureExcludedStd.append(detectorTemperature[np.where(xraysOnAndLaserOffExcluded)].std(axis=0))
                    else:
                        print "\t\tNo reference hits excluded..."
                    # xrays off
                    if (xraysOff).any():
                        runBackgroundMean.append(detectorData[np.where(xraysOff)].mean(axis=0))
                        runBackgroundStd.append(detectorData[np.where(xraysOff)].std(axis=0))
                        runBGTemperatureMean.append(detectorTemperature[np.where(xraysOff)].mean(axis=0))
                        runBGTemperatureStd.append(detectorTemperature[np.where(xraysOff)].std(axis=0))
                    else:
                        print "\t\tNo background shots included..."
                    if (xraysOffExcluded).any():
                        runBackgroundExcludedMean.append(detectorData[np.where(xraysOffExcluded)].mean(axis=0))
                        runBackgroundExcludedStd.append(detectorData[np.where(xraysOffExcluded)].std(axis=0))
                        runBGTemperatureExcludedMean.append(detectorTemperature[np.where(xraysOffExcluded)].mean(axis=0))
                        runBGTemperatureExcludedStd.append(detectorTemperature[np.where(xraysOffExcluded)].std(axis=0))
                    else:
                        print "\t\tNo background shots excluded..."
                    toc = time.time()
                    times.append(toc - tic)
                    print "\tAveraged %d hits in %.2f s" % (len(tags), toc - tic)
                    print "\tSaving averaged data to '%s'..." % outputFileName
                    tic = time.time()
                    # xrays on and laser on
                    if (xraysOnAndLaserOn).any():
                        # create data groups
                        od = outputDetectorGroup.create_group("data")
                        odd = od.create_group("data")
                        odt = od.create_group("temperature")
                        # save data to data groups
                        odd.create_dataset("mean", data=runDataMean[-1])
                        odd.create_dataset("stdev", data=runDataStd[-1])
                        odt.create_dataset("mean", data=runTemperatureMean[-1])
                        odt.create_dataset("stdev", data=runTemperatureStd[-1])
                        od.create_dataset("tags", data=tags[np.where(xraysOnAndLaserOn)])
                    if (xraysOnAndLaserOnExcluded).any():
                        # create data groups
                        oe = outputDetectorGroup.create_group("excluded_data")
                        oed = oe.create_group("data")
                        oet = oe.create_group("temperature")
                        # save data to data groups
                        oed.create_dataset("mean", data=runDataExcludedMean[-1])
                        oed.create_dataset("stdev", data=runDataExcludedStd[-1])
                        oet.create_dataset("mean", data=runTemperatureExcludedMean[-1])
                        oet.create_dataset("stdev", data=runTemperatureExcludedStd[-1])
                        oe.create_dataset("tags", data=tags[np.where(xraysOnAndLaserOnExcluded)])
                    # xrays on and laser off
                    if (xraysOnAndLaserOff).any():
                        # create data groups
                        orf = outputDetectorGroup.create_group("reference")
                        ord = orf.create_group("data")
                        ort = orf.create_group("temperature")
                        # save data to data groups
                        ord.create_dataset("mean", data=runReferenceMean[-1])
                        ord.create_dataset("stdev", data=runReferenceStd[-1])
                        ort.create_dataset("mean", data=runRefTemperatureMean[-1])
                        ort.create_dataset("stdev", data=runRefTemperatureStd[-1])
                        orf.create_dataset("tags", data=tags[np.where(xraysOnAndLaserOff)])
                    if (xraysOnAndLaserOffExcluded).any():
                        # create data groups
                        oer = outputDetectorGroup.create_group("excluded_reference")
                        oed = oer.create_group("data")
                        oet = oer.create_group("temperature")
                        # save data to data groups
                        oed.create_dataset("mean", data=runReferenceExcludedMean[-1])
                        oed.create_dataset("stdev", data=runReferenceExcludedStd[-1])
                        oet.create_dataset("mean", data=runRefTemperatureExcludedMean[-1])
                        oet.create_dataset("stdev", data=runRefTemperatureExcludedStd[-1])
                        oer.create_dataset("tags", data=tags[np.where(xraysOnAndLaserOffExcluded)])
                    # xrays off
                    if (xraysOff).any():
                        # create data groups
                        ob = outputDetectorGroup.create_group("background")
                        obd = ob.create_group("data")
                        obt = ob.create_group("temperature")
                        # save data to data groups
                        obd.create_dataset("mean", data=runBackgroundMean[-1])
                        obd.create_dataset("stdev", data=runBackgroundStd[-1])
                        obt.create_dataset("mean", data=runBGTemperatureMean[-1])
                        obt.create_dataset("stdev", data=runBGTemperatureStd[-1])
                        ob.create_dataset("tags", data=tags[np.where(xraysOff)])
                    if (xraysOffExcluded).any():
                        # create data groups
                        oeb = outputDetectorGroup.create_group("excluded_background")
                        oed = oeb.create_group("data")
                        oet = oeb.create_group("temperature")
                        # save data to data groups
                        oed.create_dataset("mean", data=runBackgroundExcludedMean[-1])
                        oed.create_dataset("stdev", data=runBackgroundExcludedStd[-1])
                        oet.create_dataset("mean", data=runBGTemperatureExcludedMean[-1])
                        oet.create_dataset("stdev", data=runBGTemperatureExcludedStd[-1])
                        oeb.create_dataset("tags", data=tags[np.where(xraysOffExcluded)])
                    toc = time.time()
                    times.append(toc - tic)
                    print "\tSaved data in %.2f s" % (toc - tic)
                else:
                    print "\tAveraging hits for '%s' as a function of %s..." % (d, motorsThatMoved[0])
                    tic = time.time()
                    # unordered set
                    #motorValueSet = set(motorValues)
                    # ordered set in original order
                    motorValuesOrderedSet = OrderedSet(motorValues)
                    motorValuesOrdered = []
                    #for m, n in zip(motorValueSet, motorValueOrderedSet):
                    #    print m, n
                    if xraysOnAndLaserOn.any():
                        taggedDetectorDataMean = []
                        taggedDetectorDataStd = []
                        taggedDetectorTemperatureMean = []
                        taggedDetectorTemperatureStd = []
                        includedTags = []
                    else:
                        print "\t\tNo laser hits included..."
                    if xraysOnAndLaserOnExcluded.any():
                        taggedDetectorDataExcludedMean = []
                        taggedDetectorDataExcludedStd = []
                        taggedDetectorTemperatureExcludedMean = []
                        taggedDetectorTemperatureExcludedStd = []
                        excludedTags = []
                    else:
                        print "\t\tNo laser hits excluded..."
                    # xrays on and laser off
                    if (xraysOnAndLaserOff).any():
                        runReferenceMean.append(detectorData[np.where(xraysOnAndLaserOff)].mean(axis=0))
                        runReferenceStd.append(detectorData[np.where(xraysOnAndLaserOff)].std(axis=0))
                        runRefTemperatureMean.append(detectorTemperature[np.where(xraysOnAndLaserOff)].mean(axis=0))
                        runRefTemperatureStd.append(detectorTemperature[np.where(xraysOnAndLaserOff)].std(axis=0))
                    else:
                        print "\t\tNo reference hits included..."
                    if (xraysOnAndLaserOffExcluded).any():
                        runReferenceExcludedMean.append(detectorData[np.where(xraysOnAndLaserOffExcluded)].mean(axis=0))
                        runReferenceExcludedStd.append(detectorData[np.where(xraysOnAndLaserOffExcluded)].std(axis=0))
                        runRefTemperatureExcludedMean.append(detectorTemperature[np.where(xraysOnAndLaserOffExcluded)].mean(axis=0))
                        runRefTemperatureExcludedStd.append(detectorTemperature[np.where(xraysOnAndLaserOffExcluded)].std(axis=0))
                    else:
                        print "\t\tNo reference hits excluded..."
                    # xrays off
                    if (xraysOff).any():
                        runBackgroundMean.append(detectorData[np.where(xraysOff)].mean(axis=0))
                        runBackgroundStd.append(detectorData[np.where(xraysOff)].std(axis=0))
                        runBGTemperatureMean.append(detectorTemperature[np.where(xraysOff)].mean(axis=0))
                        runBGTemperatureStd.append(detectorTemperature[np.where(xraysOff)].std(axis=0))
                    else:
                        print "\t\tNo background shots included..."
                    if (xraysOffExcluded).any():
                        runBackgroundExcludedMean.append(detectorData[np.where(xraysOffExcluded)].mean(axis=0))
                        runBackgroundExcludedStd.append(detectorData[np.where(xraysOffExcluded)].std(axis=0))
                        runBGTemperatureExcludedMean.append(detectorTemperature[np.where(xraysOffExcluded)].mean(axis=0))
                        runBGTemperatureExcludedStd.append(detectorTemperature[np.where(xraysOffExcluded)].std(axis=0))
                    else:
                        print "\t\tNo background shots excluded..."
                    for m in motorValuesOrderedSet:
                        if options.verbose:
                            print "\tAveraging hits where %s = %.2e..." % (motorsThatMoved[0], m)
                        motorValuesOrdered.append(m)
                        tagsWithMotorValues = motorTags[np.where(motorValues == m)]
                        indicesOfTagsWithMotorValues = np.zeros(tags.shape, dtype=np.bool)
                        # sum together all indices that return true
                        for t in tagsWithMotorValues:
                            indicesOfTagsWithMotorValues += (tags == t)
                        if xraysOnAndLaserOn.any():
                            if (indicesOfTagsWithMotorValues & xraysOnAndLaserOn).any():
                                taggedDetectorDataMean.append(detectorData[indicesOfTagsWithMotorValues & xraysOnAndLaserOn].mean(axis=0))
                                taggedDetectorDataStd.append(detectorData[indicesOfTagsWithMotorValues & xraysOnAndLaserOn].std(axis=0))
                                #taggedDetectorDataMean.append(np.mean(detectorData[np.where(indicesOfTagsWithMotorValues & xraysOnAndLaserOn)], axis=0))
                                #taggedDetectorDataStd.append(np.mean(detectorData[np.where(indicesOfTagsWithMotorValues & xraysOnAndLaserOn)], axis=0))
                                taggedDetectorTemperatureMean.append(detectorTemperature[indicesOfTagsWithMotorValues & xraysOnAndLaserOn].mean(axis=0))
                                taggedDetectorTemperatureStd.append(detectorTemperature[indicesOfTagsWithMotorValues & xraysOnAndLaserOn].std(axis=0))
                                includedTags.append([t for t in tags[indicesOfTagsWithMotorValues & xraysOnAndLaserOn]])
                            else:
                                if options.verbose:
                                    print "\t\tNo included hits for %s = %.2e, appending zeros..." % (motorsThatMoved[0], m)
                                taggedDetectorDataMean.append(np.zeros(detectorData[0].shape))
                                taggedDetectorDataStd.append(np.zeros(detectorData[0].shape))
                                taggedDetectorTemperatureMean.append(np.zeros(detectorTemperature[0].shape))
                                taggedDetectorTemperatureStd.append(np.zeros(detectorTemperature[0].shape))
                                includedTags.append(None)
                        if xraysOnAndLaserOnExcluded.any():
                            if (indicesOfTagsWithMotorValues & xraysOnAndLaserOnExcluded).any():
                                taggedDetectorDataExcludedMean.append(detectorData[indicesOfTagsWithMotorValues & xraysOnAndLaserOnExcluded].mean(axis=0))
                                taggedDetectorDataExcludedStd.append(detectorData[indicesOfTagsWithMotorValues & xraysOnAndLaserOnExcluded].std(axis=0))
                                taggedDetectorTemperatureExcludedMean.append(detectorTemperature[indicesOfTagsWithMotorValues & xraysOnAndLaserOnExcluded].mean(axis=0))
                                taggedDetectorTemperatureExcludedStd.append(detectorTemperature[indicesOfTagsWithMotorValues & xraysOnAndLaserOnExcluded].std(axis=0))
                                excludedTags.append([t for t in tags[indicesOfTagsWithMotorValues & xraysOnAndLaserOnExcluded]])
                            else:
                                if options.verbose:
                                    print "\t\tNo excluded hits for %s = %.2e, appending zeros..." % (motorsThatMoved[0], m)
                                taggedDetectorDataExcludedMean.append(np.zeros(detectorData[0].shape))
                                taggedDetectorDataExcludedStd.append(np.zeros(detectorData[0].shape))
                                taggedDetectorTemperatureExcludedMean.append(np.zeros(detectorTemperature[0].shape))
                                taggedDetectorTemperatureExcludedStd.append(np.zeros(detectorTemperature[0].shape))
                                exludedTags.append(None)
                    # save averaged data for each motor value to file
                    toc = time.time()
                    times.append(toc - tic)
                    print "\tAveraged %d hits into %d bins in %.2f s" % (len(tags), len(motorValuesOrdered), toc - tic)
                    print "\tSaving averaged data to '%s'..." % (outputFileName)
                    tic = time.time()
                    # xrays on and laser on
                    if xraysOnAndLaserOn.any():
                        om = outputDetectorGroup.create_group(motorsThatMoved[0] + "_scan")
                        excludedPositions = 0
                        for i in range(len(motorValuesOrdered)):
                            if includedTags[i]:
                                # create data groups
                                omi = om.create_group("%04d" % (i - excludedPositions))
                                od = omi.create_group("data")
                                ot = omi.create_group("temperature")
                                # save data to data groups
                                od.create_dataset("mean", data=taggedDetectorDataMean[i])
                                od.create_dataset("stdev", data=taggedDetectorDataStd[i])
                                ot.create_dataset("mean", data=taggedDetectorTemperatureMean[i])
                                ot.create_dataset("stdev", data=taggedDetectorTemperatureStd[i])
                                omi.create_dataset("tags", data=includedTags[i])
                                omi.create_dataset("motor_value", data=motorValuesOrdered[i])
                            elif (not xraysOnAndLaserOffExcluded.any()) or (not excludedTags[i]):
                                excludedPositions += 1
                    if xraysOnAndLaserOnExcluded.any():
                        om = outputDetectorGroup.require_group(motorsThatMoved[0] + "_scan")
                        excludedPositions = 0
                        for i in range(len(motorValuesOrdered)):
                            if excludedTags[i]:                            
                                # create data groups
                                omi = om.require_group('%04d' % (i - excludedPositions))
                                oe = outputDetectorGroup.create_group("excluded_data")
                                oed = oe.create_group("data")
                                oet = oe.create_group("temperature")
                                # save data to data groups
                                oed.create_dataset("mean", data=taggedDetectorDataExcludedMean[-1])
                                oed.create_dataset("stdev", data=taggedDetectorDataExcludedStd[-1])
                                oet.create_dataset("mean", data=taggedDetectorTemperatureExcludedMean[-1])
                                oet.create_dataset("stdev", data=taggedDetectorTemperatureExcludedStd[-1])
                                oe.create_dataset("tags", data=excludedTags[i])
                                omi.require_dataset("motor_value", data=motorValuesOrdered[i])
                            elif (not xraysOnAndLaserOff.any()) or (not excludedTags[i]):
                                excludedPositions += 1
                    # xrays on and laser off
                    if (xraysOnAndLaserOff).any():
                        # create data groups
                        orf = outputDetectorGroup.create_group("reference")
                        ord = orf.create_group("data")
                        ort = orf.create_group("temperature")
                        # save data to data groups
                        ord.create_dataset("mean", data=runReferenceMean[-1])
                        ord.create_dataset("stdev", data=runReferenceStd[-1])
                        ort.create_dataset("mean", data=runRefTemperatureMean[-1])
                        ort.create_dataset("stdev", data=runRefTemperatureStd[-1])
                        orf.create_dataset("tags", data=tags[np.where(xraysOnAndLaserOff)])
                    if (xraysOnAndLaserOffExcluded).any():
                        # create data groups
                        oer = outputDetectorGroup.create_group("excluded_reference")
                        oed = oer.create_group("data")
                        oet = oer.create_group("temperature")
                        # save data to data groups
                        oed.create_dataset("mean", data=runReferenceExcludedMean[-1])
                        oed.create_dataset("stdev", data=runReferenceExcludedStd[-1])
                        oet.create_dataset("mean", data=runRefTemperatureExcludedMean[-1])
                        oet.create_dataset("stdev", data=runRefTemperatureExcludedStd[-1])
                        oer.create_dataset("tags", data=tags[np.where(xraysOnAndLaserOffExcluded)])
                    # xrays off
                    if (xraysOff).any():
                        # create data groups
                        ob = outputDetectorGroup.create_group("background")
                        obd = ob.create_group("data")
                        obt = ob.create_group("temperature")
                        # save data to data groups
                        obd.create_dataset("mean", data=runBackgroundMean[-1])
                        obd.create_dataset("stdev", data=runBackgroundStd[-1])
                        obt.create_dataset("mean", data=runBGTemperatureMean[-1])
                        obt.create_dataset("stdev", data=runBGTemperatureStd[-1])
                        ob.create_dataset("tags", data=tags[np.where(xraysOff)])
                    if (xraysOffExcluded).any():
                        # create data groups
                        oeb = outputDetectorGroup.create_group("excluded_background")
                        oed = oeb.create_group("data")
                        oet = oeb.create_group("temperature")
                        # save data to data groups
                        oed.create_dataset("mean", data=runBackgroundExcludedMean[-1])
                        oed.create_dataset("stdev", data=runBackgroundExcludedStd[-1])
                        oet.create_dataset("mean", data=runBGTemperatureExcludedMean[-1])
                        oet.create_dataset("stdev", data=runBGTemperatureExcludedStd[-1])
                        oeb.create_dataset("tags", data=tags[np.where(xraysOffExcluded)])
                    toc = time.time()
                    times.append(toc - tic)
                    print "\tSaved data in %.2f s" % (toc - tic)
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
