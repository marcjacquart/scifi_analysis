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
#from array import array
import os
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import numpy as np

# Root:
import ROOT
import rootUtils as ut

# SND:
#import shipunit as u
import SndlhcGeo
import SndlhcTracking

# Custom functions defined in external file:
from analysisFunctions import goodEvent, valueToColor


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
lsOfGlobals.Add(geo.modules['Scifi'])
# Don't need the Digi_MuFilterHit because we only look at SciFi events.
# lsOfGlobals.Add(geo.modules['MuFilter'])


# Extract the TTree: rawConv is the name of the tree in the root file
eventTree = rootFile.rawConv

# Setup tracking:
trackTask = SndlhcTracking.Tracking()
trackTask.InitTask(eventTree)

# GetCurrentNavigator() Returns current navigator for the calling thread. 
# A navigator is a class providing navigation API for TGeo geometries.
# Need it to browse trough the geometry files.
nav = ROOT.gGeoManager.GetCurrentNavigator()


# Loop on the event at the branchs level:
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
lenForHist = []


offset = ROOT.TVector3(47.8,-15.3,16.5) # To have coordinates between 0 and 40.
for sTree in eventTree: # sTree == single tree for one event
    # digis = []
    # if sTree.FindBranch("Digi_ScifiHits"): 
    #     digis.append(sTree.Digi_ScifiHits) # Digi_ScifiHits is a branch
    #     lenForHist.append(len(sTree.Digi_ScifiHits))
    # if sTree.FindBranch("EventHeader"):
    #     T = sTree.EventHeader.GetEventTime()

    # scifiCluster() only works when in the loop O_o
    clusterArr = trackTask.scifiCluster() # Warning: cluster type sndCluster
    print(f'Clustering reduced {len(sTree.Digi_ScifiHits)} hits '
        + f'into {len(clusterArr)} clusters.')


    trackTask.ExecuteTask()
    print(trackTask.event.fittedTracks)
    
            

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
    fig= plt.figure(figsize = (10, 7))
    ax = plt.axes(projection="3d")

    for hitNumber in range(len(arrPosStart)):
        ax.plot(
            xs = [arrPosStart[hitNumber][0], arrPosStop[hitNumber][0]], 
            ys = [arrPosStart[hitNumber][1], arrPosStop[hitNumber][1]],
            zs = [arrPosStart[hitNumber][2], arrPosStop[hitNumber][2]],
            ls = '-',
            # RGB format to color different Scifi planes
            color = valueToColor(abs(arrPosStart[hitNumber][2])) )
    # Fit infos
    fitArr = []
    for aTrack in trackTask.event.fittedTracks:
        for i in range(aTrack.getNumPointsWithMeasurement()):
            state = aTrack.getFittedState(i)
            pos = state.getPos() + offset
            fitArr.append(pos)

    ax.plot(
        xs = [element[0] for element in fitArr], 
        ys = [element[1] for element in fitArr],
        zs = [element[2] for element in fitArr],
        color = 'r',
        label = 'fit')
    ax.set_xlabel('x [cm]')
    ax.set_ylabel('y [cm]')
    ax.set_zlabel('z [cm]')
    plt.legend()
    plt.show()
    plt.close()





# ---------- end of script ------------

# 3d scatter plot code:
if False:
    fig= plt.figure(figsize = (10, 7))
    ax = plt.axes(projection="3d")

    
    ax.scatter3D(
        xs = [x[0] for x in arrPosStart], 
        ys = [x[1] for x in arrPosStart],
        zs = [x[2] for x in arrPosStart],
        color = 'b',
        label = 'A')
    # [:,0] notation is only for np arrays. else use [x[0] for x in arrPosStart]
    ax.scatter3D(
        xs = [x[0] for x in arrPosStop],
        ys = [x[1] for x in arrPosStop],
        zs = [x[2] for x in arrPosStop],
        color = 'r',
        label = 'B')
    plt.legend()
    plt.show()
    plt.close()

    # Position of the single hits before clustering:
    for hit in sTree.Digi_ScifiHits:
        detID = hit.GetDetectorID()
        geo.modules['Scifi'].GetSiPMPosition(detID,A,B)