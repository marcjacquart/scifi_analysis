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
from analysisFunctions import (goodEvent, zPlaneArr,
    crossAllPlanes, indexStationsHit, chi2Hist, planesHist)
from plotFunctions import display3dTrack, display2dTrack


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

# z positions of the planes, array of 10 elements
zArr = zPlaneArr(eventTree=eventTree, geo=geo)
print(zArr)

chi2_nDfArr = [] # To fill for histogram
nPlanesHit = []

for sTree in eventTree: # sTree == single tree for one event
    # scifiCluster() only works when in the loop O_o?
    clusterArr = trackTask.scifiCluster() # Warning: cluster type sndCluster
    #print(f'Clustering reduced {len(sTree.Digi_ScifiHits)} hits '
    #    + f'into {len(clusterArr)} clusters.')

    trackTask.ExecuteTask() # Clusters stored in trackTask.clusters
    
    if goodEvent(eventTree = eventTree, nStations = 4, allowMore = True):
        # Put the event to the dict format to use fitTrack()
        hitList = {}
        for x in clusterArr:
            A,B = ROOT.TVector3(),ROOT.TVector3()
            #x.GetPosition(A,B)
            clusterID = x.GetFirst()
            dictEntery = {clusterID: x}
            hitList.update(dictEntery)
        fittedTrack = trackTask.fitTrack(hitlist=hitList) # type: genfit::Track
        fitStatus   = fittedTrack.getFitStatus()
        if fitStatus.getNdf() == 0.0:
            # Perfect fit because no degree of freedom
            # Put impossiblie chi2 -1 value to put them apart
            chi2_nDf = -1
        else:
            chi2_nDf = fitStatus.getChi2()/fitStatus.getNdf() 
        chi2_nDfArr.append(chi2_nDf)
        nPlanesHit.append(len(indexStationsHit(eventTree)))
        #print(f'Fit points: {fittedTrack.getPoints()}')
        #print(f'chi2_nDf: {chi2_nDf}')

        fitHits =[ROOT.TVector3()]*10 # Points where fit crosses the 10 planes.
        # for i in range(fittedTrack.getNumPointsWithMeasurement()):
        # Only the first can be used, other hits give same line
        state = fittedTrack.getFittedState(0)
        pos = state.getPos()
        mom = state.getMom()
        # linear fit: pos + lambda * mom
        for planeIndex in range(len(zArr)):
            lambdaPlane = (zArr[planeIndex] - pos[2]) / mom[2]
            fitHits[planeIndex] = pos + lambdaPlane * mom

        if crossAllPlanes(fitHitsArr=fitHits, geo=geo, verbose=False):
            indexHitArr = indexStationsHit(eventTree = eventTree)
            nPlanesMissed = 10-len(indexHitArr)

            #hitsMissed = fitHits
            hitsMissed = [ROOT.TVector3(x[0],x[1],x[2]) for x in fitHits]
            for index in indexHitArr:
                # Remove the coordinate of the planes hit,
                # only the coordinates of the missed hits remains.
                hitsMissed[index] = ROOT.TVector3(-1.68283565985,-1,50)
            hitsMissed = [item for item in hitsMissed if item[0] != -1.68283565985]
            #print(len(hitsMissed))

            arrPosStart = []
            arrPosStop = []
            for cluster in clusterArr:
                # A: beginning, B: end of the activated scintillating bar.
                A,B = ROOT.TVector3(),ROOT.TVector3()
                cluster.GetPosition(A,B) # Fill A and B directly with position.
                # hit.isVertical(): True if vertical fibers measuring x coord.
                arrPosStart.append(A)            
                arrPosStop.append(B)
                
            if chi2_nDf<20: # display 3d trajectories with condition
                if True:
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
        else:
            # if hasn't cross al planes
            pass
    # Else if don't pass event selection with number station hit:
    else:
        # print('Bad event!')
        pass


chi2Hist(chi2_nDfArr)