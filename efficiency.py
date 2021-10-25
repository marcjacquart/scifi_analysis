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
import matplotlib.pyplot as plt

# Root:
import ROOT
import rootUtils as ut

# SND:
import SndlhcGeo
import SndlhcTracking

# Custom functions defined in external file:
from analysisFunctions import display3dTrack, goodEvent, crossAllPlanes, indexStationsHit


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


zArr = [0,0,0,0,0,0,0,0,0,0] # Fill Z coordinate of the planes
# format: [1hor, 1ver, 2hor, ...,5hor, 5ver]
# /!\ THIS ASSUME Z FIX FOR ALL THE PLANE!
# Warning message if z values > epsilon for different SiPM of same plane
epsilon = 0.0001 
chi2_nDfArr = [] # To fill for histogram
nPlanesHit = []
for sTree in eventTree: # Need the first event with all planes hit:
    if goodEvent(eventTree = eventTree, nStations = 5, allowMore = False):
        A,B = ROOT.TVector3(),ROOT.TVector3()
        for hit in sTree.Digi_ScifiHits:
            hitID = hit.GetChannelID()
            geo.modules['Scifi'].GetSiPMPosition(hitID,A,B)
            indexArr = (hitID // 1000000-1) * 2  + (hitID // 100000) % 2
            #           add 2 for each plane      add 1 if vertical

            zVal = 0.5 * (A[2] + B[2])
            #Will overwrite if event of same plane, good to check if same z.
            if zArr[indexArr]!=0 and zArr[indexArr]-zVal>epsilon:
                print(f'WARNING: SciFi planes {indexArr} not aligned in Z direction!')
            zArr[indexArr] = zVal # Fill z array
        break
#print(f'zArr: {zArr}')


offset = ROOT.TVector3(0,0,0) #(47.8,-15.3,16.5) # To have coordinates between 0 and 40.
for sTree in eventTree: # sTree == single tree for one event
    # scifiCluster() only works when in the loop O_o?
    clusterArr = trackTask.scifiCluster() # Warning: cluster type sndCluster
    #print(f'Clustering reduced {len(sTree.Digi_ScifiHits)} hits '
    #    + f'into {len(clusterArr)} clusters.')

    trackTask.ExecuteTask()
    # Clusters are stored in the list: trackTask.clusters
    
    # <sndCluster>.GetFirst() gives the first sipm channel number: STMRFFF
    # First digit S:       station # within the sub-detector 1-5
    # Second digit T:      type of the plane: 0-horizontal fiber plane, 
    #                      1-vertical fiber plane
    # Third digit M:       determines the mat number 0-2
    # Fourth digit S:      SiPM number  0-3
    # Last three digits F: local SiPM channel number in one mat  0-127

    # Chi2/nDOF to measure the fit (DOF: degrees of freedom of the fit)
    
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
                
                # Apply offset and put in array
                arrPosStart.append(A + offset)            
                arrPosStop.append(B + offset)
            if chi2_nDf>150: # display 3d trajectories with condition
                
                display3dTrack(
                    arrPosStart = arrPosStart, 
                    arrPosStop = arrPosStop, 
                    trackTask = trackTask,
                    offset = offset,
                    fitHits = hitsMissed)
        else:
            # if hasn't cross al planes
            pass
    # Else if don't pass event selection with number station hit:
    else:
        # print('Bad event!')
        pass

binsArr = np.linspace(0,5000,5000)
fig, ax =plt.subplots(figsize=(6,8), dpi=300, tight_layout=True)
ax.hist(chi2_nDfArr, bins=binsArr)
ax.set_xlim(left=0.0,right=5000)
plt.xlabel('chi2/dof')
plt.ylabel('Number of events')
plt.show()
plt.close()

fig, ax =plt.subplots(figsize=(6,8), dpi=300, tight_layout=True)
ax.hist(nPlanesHit)
plt.xlabel('Number of planes hit')
plt.ylabel('Number of events')
plt.show()
plt.close()