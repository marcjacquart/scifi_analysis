# To launch:
# python3 efficiency.py -p /eos/experiment/sndlhc/testbeam/scifi/sndsw/ -f sndsw_raw_000001.root -g geofile_full.Ntuple-TGeant4.root


# This script:
#    - Import the raw data and geomery files
#    - Use the geometry to compute cluster positions
#    - Select events with station hit number, chi2 
#      or constrained tracks within the detector.
#    - Fit and display the events.
# Uses functions from analysisFunctions.py and plotFunctions.py


# Imports:
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
    crossAllPlanes, indexStationsHit, sortHitStation)
from plotFunctions import display3dTrack, display2dTrack, chi2Hist, planesHist, diffHist, allPlanesGauss


# Paths+name of the data and geometry root files as script arguments:
parser = ArgumentParser()
parser.add_argument("-p", "--path", dest="path", help="run number",required=True,default="")
parser.add_argument("-f", "--inputFile", dest="inputFile", help="input root file",default="",required=True)
parser.add_argument("-g", "--geoFile", dest="geoFile", help="geometry file", required=True)
options = parser.parse_args()

# Open root and geometry files:
rootFile = ROOT.TFile.Open(options.path+options.inputFile)
geo = SndlhcGeo.GeoInterface(options.path+options.geoFile)

#Global variables definitions:
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
zArr = zPlaneArr(eventTree=eventTree, geo=geo)
#print(zArr)

chi2_nDfArr = [] # To fill in the loop for histogram
nPlanesHit = [] # Same

gaussFitArr = []
fitStationsArr = [[2,3,4,5],[1,3,4,5],[1,2,4,5],[1,2,3,5],[1,2,3,4]]
for testStationNum in [1,2,3,4,5]:
    horDiffArr = []
    verDiffArr = []
    # Loop over the individual events:
    for sTree in eventTree: # sTree == single tree for one event 
        # 1: Select events with given number of stations hit:   
        if goodEvent(eventTree = eventTree, nStations = 5, allowMore = True):
           
            # # Test cluster problem:
            # print('##################')
            # IDArr = []
            # for hit in sTree.Digi_ScifiHits:
            #     detID = int(hit.GetDetectorID())
            #     IDArr.append(detID)
            # IDArr.sort()
            # prevElement=0
            # for element in IDArr:
            #     diff = element - prevElement
            #     if not (diff ==1):
            #         if diff < 10:
            #             print(f'!!!---{element - prevElement}---!!!')
            #         else:
            #             print('---')
            #     print(element)
            #     prevElement = element


            # Fill .clusters, .trackCandidates and .event.fittedTracks:
            trackTask.ExecuteTask() # Clusters stored in trackTask.clusters

            # Get fit infos on the full hit list:
            fittedTrack = trackTask.event.fittedTracks[0] # ASK THOMAS WHY LIST
            #print(len(trackTask.event.fittedTracks)) # Why only 1?
            fitStatus = fittedTrack.getFitStatus()
            
            # Or fit only with a subset of the stations:     
            FitStations = fitStationsArr[testStationNum-1] # 0-4 array boundaries
            clusterArr = trackTask.clusters
            clusFit, clusTest = sortHitStation(clusterArr=clusterArr,stationArr=FitStations)
            # Do the manual fit with the reduced number of stations:
            # Need manual dictionnary for the fitTrack(hitlist):
            clusDict = {}
            for x in clusFit:
                A,B = ROOT.TVector3(),ROOT.TVector3()
                clusterID = x.GetFirst()
                dictEntery = {clusterID: x}
                clusDict.update(dictEntery)
            reducedFit = trackTask.fitTrack(hitlist=clusDict)
            reducesFitStatus= reducedFit.getFitStatus()

            fitHitsReduced = extendHits(fittedTrack=reducedFit, zArr=zArr)
            horDiff, verDiff = distFit(fitHits=fitHitsReduced, clusterArr=clusterArr, testStationNum=testStationNum)
            #print(f'Diff: {horDiff}, {verDiff}')
            if horDiff != 1000: # Don't append missing hits
                horDiffArr.append(horDiff)
            if verDiff != 1000:
                verDiffArr.append(verDiff)
            # Extend the fit crossing point to all the planes:
            fitHits = extendHits(fittedTrack=fittedTrack, zArr=zArr)

            # for i in range(len(fitHits)): # Check that it gives different fits:
            #     print(f'Fit diff in x:{fitHitsReduced[i][0]-fitHits[i][0]}')
            #     print(f'Fit diff in y:{fitHitsReduced[i][1]-fitHits[i][1]}')
            #     print(f'Fit diff in z:{fitHitsReduced[i][2]-fitHits[i][2]}')

            # 2: Select only trajectory crossing all the planes:
            if crossAllPlanes(fitHitsArr=fitHits, geo=geo, verbose=False):
                
                indexHitArr = indexStationsHit(eventTree = eventTree)
                hitsMissed = [ROOT.TVector3(x[0],x[1],x[2]) for x in fitHits] # Copy values, not pointers
                for index in indexHitArr:
                    # Remove the coordinate of the planes hit,
                    # only the coordinates of the missed hits remains.
                    hitsMissed[index] = ROOT.TVector3(-1.68283565985,0,0) # Flag to remove
                hitsMissed = [item for item in hitsMissed if item[0] != -1.68283565985]

                # Chi2:
                if fitStatus.getNdf() == 0.0:
                    # Perfect fit because no degree of freedom
                    # Put impossiblie chi2 -1 value to put them apart
                    chi2_nDf = -1
                else:
                    chi2_nDf = fitStatus.getChi2()/fitStatus.getNdf() 

                # Append for histtograms (nPlanes after chi2 selection?)
                chi2_nDfArr.append(chi2_nDf)
                nPlanesHit.append(len(indexStationsHit(eventTree)))

                # 3: Select only low chi2 tracks: good fit and no secondary events    
                if chi2_nDf<20: # display 3d trajectories with condition
                   
                    if False: # Plot trajectories                    
                        arrPosStart = []
                        arrPosStop = []
                        for cluster in trackTask.clusters:
                            # A: beginning, B: end of the activated scintillating bar.
                            A,B = ROOT.TVector3(),ROOT.TVector3()
                            cluster.GetPosition(A,B) # Fill A and B directly with position.
                            # hit.isVertical(): True if vertical fibers measuring x coord.
                            arrPosStart.append(A)            
                            arrPosStop.append(B)
                        display3dTrack(
                            arrPosStart = arrPosStart, 
                            arrPosStop = arrPosStop, 
                            trackTask = trackTask,
                            fitHits = hitsMissed)
                        # display2dTrack(
                        #     arrPosStart = arrPosStart, 
                        #     arrPosStop = arrPosStop,
                        #     trackTask = trackTask,
                        #     fitHits = hitsMissed)


    resultFit = diffHist(horDiffArr=horDiffArr, verDiffArr=verDiffArr, stationNum=testStationNum)
    gaussFitArr.append(resultFit)
    print(f'Test station: {testStationNum}')

allPlanesGauss(fitArr=gaussFitArr)
if False:
    pass
    # chi2Hist(chi2_nDfArr=chi2_nDfArr)
    # planesHist(nPlanesHit=nPlanesHit)
    