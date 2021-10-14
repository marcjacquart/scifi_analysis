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

# Custom functions defined in external file:
from analysisFunctions import goodEvent


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
#lsOfGlobals.Add(geo.modules['MuFilter'])

# Don't need the Digi_MuFilterHit ?????

# Extract the TTree: rawConv is the name of the tree in the root file
eventTree = rootFile.rawConv


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
lenForHist=[]
for sTree in eventTree: # sTree == single tree for one event?
    digis = []
    if sTree.FindBranch("Digi_ScifiHits"): 
        digis.append(sTree.Digi_ScifiHits) # Digi_ScifiHits is a branch
        lenForHist.append(len(sTree.Digi_ScifiHits))
    if sTree.FindBranch("EventHeader"):
        T = sTree.EventHeader.GetEventTime()
    for D in digis:
        aTab = []
        bTab = []
        for hit in D:
            A,B = ROOT.TVector3(),ROOT.TVector3()
            detID = hit.GetDetectorID()
            geo.modules['Scifi'].GetSiPMPosition(detID,A,B)
            A = [A[0], A[1], A[2]]
            B = [B[0], B[1], B[2]]
            aTab.append(A)            
            bTab.append(B)

        fig= plt.figure(figsize = (10, 7))
        ax = plt.axes(projection="3d")
        

        ax.scatter3D(
            xs = [x[0] for x in aTab], # /!\ [:,0] notation is only ok for np arrays. else use [x[0] for x in aTab]
            ys = [x[1] for x in aTab],
            zs = [x[2] for x in aTab],
            color = 'b',
            label = 'A')
        ax.scatter3D(
            xs = [x[0] for x in bTab],
            ys = [x[1] for x in bTab],
            zs = [x[2] for x in bTab],
            color = 'r',
            label = 'B')
        plt.legend()
        plt.show()
        plt.close()
    print(f'---------------------------------------------------')

#print(f'len(digis)={len(digis)}')


# for D in digis: # D: TClonesArray
#     # print(f'len(D)={len(D)}')
#     # print(f'---------------------------------------------------')
#     for hit in D:
#         #print(f'len(digi)={len(digi)}')
#         detID = hit.GetDetectorID()
#         # print(digi)
#         geo.modules['Scifi'].GetSiPMPosition(detID,A,B)


