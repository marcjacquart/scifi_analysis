# To launch this script:
# python3 efficiency.py -p /eos/experiment/sndlhc/testbeam/scifi/sndsw/ -f sndsw_raw_000001.root -g geofile_full.Ntuple-TGeant4.root

# Create a "figures" folder where the script is executed to save them all here.

# This script:
#    - Import the raw data and geomery files
#    - Use the geometry to compute cluster positions
#    - Select events with station hit number, chi2 
#      or constrained tracks within the detector.
#    - Fit and display the events.
#    - Compute the SciFi misalignment, show individual events
# Uses functions from analysisFunctions.py and plotFunctions.py

# General
from argparse import ArgumentParser
import os
import numpy as np
import matplotlib.pyplot as plt
# Root:
import ROOT
import rootUtils as ut
# SND:
import SndlhcGeo
import SndlhcTracking
# Custom functions defined in external file:
from analysisFunctions import (goodEvent, zPlaneArr, extendHits, distFit,
    crossAllPlanes, indexStationsHit, sortHitStation, testClusterProblem, 
    customFitStatus)
from plotFunctions import (display3dTrack, display2dTrack, chi2Hist, 
    planesHist, diffHist, allPlanesGauss, diffPosHist, rotationAngle)

# Script options:
displayTrack = True # Display 2d/3d track + fit plots of individual events.
fitReducedStations = False # True to loop 4-fit, 1-test, False to fit with all 5 stations.
needXYInfo = True # Both vertical and horizontal planes must be hit (for rotation)
# /!\ Be careful to select 5 stations with goodEvent() if fitReducedStations = True


# Paths+name of the data and geometry root files as script arguments:
parser = ArgumentParser()
parser.add_argument(
    "-p", 
    "--path", 
    dest = "path", 
    help = "run number",
    required = True,default="")
parser.add_argument(
    "-f", 
    "--inputFile", 
    dest = "inputFile", 
    help = "input root file",
    default = "",required=True)
parser.add_argument(
    "-g", 
    "--geoFile", 
    dest = "geoFile", 
    help = "geometry file", 
    required = True)
options = parser.parse_args()

# Open root and geometry files:
rootFile = ROOT.TFile.Open(options.path+options.inputFile)
geo = SndlhcGeo.GeoInterface(options.path+options.geoFile)

# Global variables definitions:
lsOfGlobals = ROOT.gROOT.GetListOfGlobals()
lsOfGlobals.Add(geo.modules['Scifi']) # only look at SciFi events

# Extract the TTree: rawConv is the name of the tree in the root file
eventTree = rootFile.rawConv

# Setup tracking:
trackTask = SndlhcTracking.Tracking()
trackTask.InitTask(eventTree)

# GetCurrentNavigator() to browse trough the geometry files (TGeo class).
nav = ROOT.gGeoManager.GetCurrentNavigator()

# z positions of the planes, array of 10 float
zArr = zPlaneArr(eventTree = eventTree, geo = geo)

# For rotation plot to fill with each test station iteration:
if needXYInfo:
    xPosyOff_Slope = []
    xPosyOff_SlopeErr = []
    yPosxOff_Slope = []
    yPosxOff_SlopeErr = []

# Fit on 4 stations and use the 5th one as a test:
if fitReducedStations:
    gaussFitArr = []
    fitStationsArr = [[2,3,4,5],[1,3,4,5],[1,2,4,5],[1,2,3,5],[1,2,3,4]]
    testStationArr = [1,2,3,4,5]
else: # Only one loop iteration with all the stations to fit.
    fitStationsArr = [[1,2,3,4,5]]
    testStationArr =  [0] # Doesn't exist, but still need for plot names.

for testStationNum in testStationArr:
    chi2_nDfArr = [] # To fill in the loop for histogram
    nPlanesHit = [] # Same

    horDiffArr = []
    verDiffArr = []
    horPosArr = []
    verPosArr = []

    # Loop over the individual events:
    for sTree in eventTree: # sTree == single tree for one event
        # 1: Select events with given number of stations hit:
        if goodEvent(eventTree = eventTree, nStations = 5, allowMore = False):
            # testClusterProblem(eventTree = sTree) # Test cluster spacing
            fit, fitStatus = customFitStatus(
                trackTask = trackTask, 
                FitStations = fitStationsArr[testStationNum-1])

            # Extend the fit crossing point to all the planes:
            fitHitsExt = extendHits(fittedTrack = fit, zArr = zArr)

            # 2: Select only events with trajectory crossing all the planes:
            if crossAllPlanes(fitHitsArr = fitHitsExt, geo = geo):
                
                indexHitArr = indexStationsHit(eventTree = eventTree)
                # Copy values, not pointers
                hitsMissed = [ROOT.TVector3(x[0],x[1],x[2]) for x in fitHitsExt]
                for index in indexHitArr:
                    # Remove the coordinate of the planes hit,
                    # only the coordinates of the missed hits remains.
                    # -1.12312312312: Flag to remove the events
                    hitsMissed[index] = ROOT.TVector3(-1.12312312312,0,0)
                hitsMissed = [item for item in hitsMissed if item[0] != -1.12312312312]

                # Chi2:
                if fitStatus.getNdf() == 0.0: # Fit has no degree of freedom
                    chi2_nDf = -1  # Impossible value to put them aside
                else: 
                    chi2_nDf = fitStatus.getChi2()/fitStatus.getNdf() 

                # Append for histograms: (maybe append after chi2 selection?)
                chi2_nDfArr.append(chi2_nDf)
                

                # 3: Select low chi2 tracks: good fit and no secondary events.   
                # Exception: unhandled, unknown C++ exception" When dof=0, 
                # also need to filter chi2_nDf = -1 cases.
                if chi2_nDf<30 and chi2_nDf>=0:
                    nPlanesHit.append(len(indexStationsHit(eventTree)))
                    if fitReducedStations:
                        # Compute difference between hits and fit:
                        horDiff, verDiff, horPos, verPos= distFit(
                            fitHits = fitHitsExt,
                            clusterArr = trackTask.clusters,
                            testStationNum = testStationNum)

                        # Select event if both planes are hit
                        if needXYInfo:
                            # Search within 1cm from the fit
                            if abs(horDiff) <1 and abs(verDiff) <1:
                                horDiffArr.append(horDiff)
                                horPosArr.append(horPos)
                                verDiffArr.append(verDiff)
                                verPosArr.append(verPos)
                        else: # Or at least one plane
                            if horDiff <1: # Don't append missing hits
                                horDiffArr.append(horDiff)
                                horPosArr.append(horPos)
                            if verDiff <1:
                                verDiffArr.append(verDiff)
                                verPosArr.append(verPos)

                    if displayTrack:                  
                        arrPosStart = []
                        arrPosStop = []
                        for cluster in trackTask.clusters:
                            # A: beginning, B: end of the scintillating bar.
                            A,B = ROOT.TVector3(),ROOT.TVector3()
                            # Fill A and B directly with position:
                            cluster.GetPosition(A,B)
                            arrPosStart.append(A)            
                            arrPosStop.append(B)
                        display3dTrack(
                            arrPosStart = arrPosStart, 
                            arrPosStop = arrPosStop, 
                            trackTask = trackTask,
                            fitHits = hitsMissed)
                        display2dTrack(
                            arrPosStart = arrPosStart, 
                            arrPosStop = arrPosStop,
                            trackTask = trackTask,
                            fitHits = hitsMissed)

    if fitReducedStations:
        resultFit = diffHist(
            horDiffArr = horDiffArr,
            verDiffArr = verDiffArr,
            stationNum = testStationNum)
        gaussFitArr.append(resultFit) # in mm since diffHist() output changed

        # Difference between hit and fit vs position within the same plane:
        diffPosHist(
            posArr = verPosArr,
            diffArr = verDiffArr,
            binsPos = np.linspace(-50,-5,90),
            fileName = f'OFFSET_ver_testSta{testStationNum}',
            labels = ['X cluster position [mm]', 'X offset [mm]'],
            isCrossed = False)
        diffPosHist(
            posArr = horPosArr,
            diffArr = horDiffArr,
            binsPos = np.linspace(15,60,90),
            fileName = f'OFFSET_hor_testSta{testStationNum}',
            labels = ['Y cluster position [mm]', 'X offset [mm]'],
            isCrossed = False)

        # If only hit one plane, x and y don't have same dimensions.
        # We need the stronger all planes hit condition.
        if needXYInfo: 
            # Cross terms to check rotation misalignments
            slope, slopeErr = diffPosHist(
                posArr = verPosArr,
                diffArr = horDiffArr,
                binsPos = np.linspace(-50,-5,90),
                fileName = f'ROT_horDiff_verPos_testSta{testStationNum}',
                labels = ['X cluster position [mm]', 'Y offset [mm]'],
                isCrossed = True)
            xPosyOff_Slope.append(slope)
            xPosyOff_SlopeErr.append(slopeErr)

            slope, slopeErr = diffPosHist(
                posArr = horPosArr,
                diffArr = verDiffArr,
                binsPos = np.linspace(15,60,90),
                fileName = f'ROT_verDiff_horPos_testSta{testStationNum}',
                labels = ['Y cluster position [mm]', 'X offset [mm]'],
                isCrossed = True)
            yPosxOff_Slope.append(slope)
            yPosxOff_SlopeErr.append(slopeErr)

    chi2Hist(chi2_nDfArr=chi2_nDfArr, stationNum=testStationNum)
    #planesHist(nPlanesHit=nPlanesHit) # Show number of planes hit histogram.
    print(f'Test station: {testStationNum}') # Info of script speed

if fitReducedStations:
    # Translation misalignment of each station:
    allPlanesGauss(fitArr=gaussFitArr)

if needXYInfo:
    # Rotation misalignment of each station:
    rotationAngle(
    xPosyOff_Slope = xPosyOff_Slope,
    xPosyOff_SlopeErr = xPosyOff_SlopeErr,
    yPosxOff_Slope = yPosxOff_Slope,
    yPosxOff_SlopeErr = yPosxOff_SlopeErr)