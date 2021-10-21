# To launch:
# python3 efficiency.py -p /eos/experiment/sndlhc/testbeam/scifi/sndsw/ -f sndsw_raw_000001.root -g geofile_full.Ntuple-TGeant4.root


# This script:
#    - Import the raw data and geomery files
#    - Use the geometry o compute cluster positions
#    - Do a linear fit using 4 planes
#    - Tries to find the 5th along the fit to measure the efficiency
# It uses the functions defined in analysisFunctions.py


# Imports:
# General
from argparse import ArgumentParser
import os
import numpy as np

# Root:
import ROOT
import rootUtils as ut

# SND:
import SndlhcGeo
import SndlhcTracking

# Custom functions defined in external file:
from analysisFunctions import display3dTrack, goodEvent


# Paths+name of the data and geometry root files passed to the script as
# mendatory arguments.
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

# GetCurrentNavigator() Returns current navigator for the calling thread. 
# A navigator is a class providing navigation API for TGeo geometries.
# Need it to browse trough the geometry files.
nav = ROOT.gGeoManager.GetCurrentNavigator()

# 3 branchs are in the data root file:
#   - EventHeader with the following infos:
#       - fEventTime about the time of the event (8e9 between pulses)
#   - Digi_ScifiHits with the following infos:
#       - fDetectorID (Digi_ScifiHits.fDetectorID is the full leaf name)
#       - signals
#       - times
#       - nphe_min / nphe_max
#       - ly_loss_params
#       - @size
#   - Digi_MuFilterHit which is empty due to only SciFi being tested here.


offset = ROOT.TVector3(47.8,-15.3,16.5) # To have coordinates between 0 and 40.
for sTree in eventTree: # sTree == single tree for one event
    # scifiCluster() only works when in the loop O_o?
    clusterArr = trackTask.scifiCluster() # Warning: cluster type sndCluster
    print(f'Clustering reduced {len(sTree.Digi_ScifiHits)} hits '
        + f'into {len(clusterArr)} clusters.')

    
    trackTask.ExecuteTask()
    # Clusters are stored in the list: trackTask.clusters
    
    # <sndCluster>.GetFirst() gives the first sipm channel number: STMRFFF
    # First digit S:       station # within the sub-detector
    # Second digit T:      type of the plane: 0-horizontal fiber plane, 
    #                      1-vertical fiber plane
    # Third digit M:       determines the mat number 0-2
    # Fourth digit S:      SiPM number  0-3
    # Last three digits F: local SiPM channel number in one mat  0-127

    # Chi2/nDOF to measure the fit (DOF: degrees of freedom of the fit)
    
    print(clusterArr)
    hitList = {}
    for x in clusterArr:
        
        A,B = ROOT.TVector3(),ROOT.TVector3()
        x.GetPosition(A,B)
        clusterID = x.GetFirst()
        #clusterid = x.fFirst
        dictEntery = {clusterID: [A,B]}
        hitList.update(dictEntery)

    fittedTrack = trackTask.fitTrack(hitlist=hitList)
    fitStatus   = fittedTrack.event.fittedTracks.getFitStatus()
    chi2 = fitStatus.getChi2()/fitStatus.getNdf() 
    print(f'chi2: {chi2}')

    arrPosStart = []
    arrPosStop = []
    for cluster in clusterArr:
        # A: beginning, B: end of the activated scintillating bar.
        A,B = ROOT.TVector3(),ROOT.TVector3()
        cluster.GetPosition(A,B) # Fill A and B directly with position.
        # hit.isVertical(): True if vertical fibers measuring x coord.
        
        # Apply offset and put in array
        arrPosStart.append(A + offset)            
        arrPosStop.append(B + offset)
    if True: # Put True to display 3d trajectories
        display3dTrack(
            arrPosStart = arrPosStart, 
            arrPosStop = arrPosStop, 
            trackTask = trackTask,
            offset = offset)
    

