import ROOT
import numpy as np


def goodEvent(eventTree, nStations, allowMore):
    '''Return True if nStations (or more if allowMore) are hits in the event'''

    stations = {}
    for detectedHit in eventTree.Digi_ScifiHits:
        stations[detectedHit.GetDetectorID()//1000000] = 1
               
    # for detectedHit in eventTree.Digi_MuFilterHit:
    #     plane = 100*(detectedHit.GetDetectorID()//1000)
    #     stations[plane] = 1
    if allowMore:
        if len(stations) >= nStations: return True
    else:
        if len(stations) == nStations: return True
    return False

def indexStationsHit(eventTree):
    '''
    return array [X_1,X_2,X_3,...] of planes hit.
    X: 0-9 number of the plane hit. With same convention as zArr:
    2 * (plane number-1) + 1 for vertical plane
    i.e. plane 2 vertical would have number 3 here. 
    '''

    # Filling a dictionary allows overwritting to have one entery per plane.
    stations = {} 
    for detectedHit in eventTree.Digi_ScifiHits:
        detID = detectedHit.GetDetectorID()

        planeNumber = (detID // 1000000 - 1) * 2  + (detID // 100000) % 2
        stations[detID//100000] = planeNumber
    
    # Fill array to return from dict values:
    indexHitArr = [] 
    for planeID in stations.values():
        indexHitArr.append(planeID)
    return indexHitArr


def extendHits(fittedTrack, zArr):
    '''
    Extend the fit position to include missing planes hits.
    fittedTrack: ROOT.genfit.Track, obtained with customFitStatus()
    zArr: array containing z position (float) of each plane. 
    Return array of 10 TVector3, aech containing the fit intersection
    with one plane.
    '''
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
    return fitHits


def crossAllPlanes(fitHitsArr,geo, verbose=False):
    '''
    Return true if fit is within the detector boundaries.
    Tell where the fit exit the detector if verbose activated
    fitHitsArr: coordinates of fit-plane intersection from extendHits().
    geo: SndlhcGeo.GeoInterface(<geometry file>) for detector boundaries.
    '''
    isInside = True # Set to False if only once out of bounds
    A,B = ROOT.TVector3(),ROOT.TVector3()

    #First plane horizontal:
    geo.modules['Scifi'].GetSiPMPosition(1000000,A,B)
    if fitHitsArr[0][0]<B[0] or A[0]<fitHitsArr[0][0]:
        isInside = False
        if verbose:
            print('first x border exceeded:')
            print(f'{fitHitsArr[0][0]}<{B[0]} or {A[0]}<{fitHitsArr[0][0]}')

    #First plane vertical:
    geo.modules['Scifi'].GetSiPMPosition(1100000,A,B)
    if fitHitsArr[1][1]<A[1] or B[1]<fitHitsArr[1][1]:
        isInside = False
        if verbose:
            print('first y border exceeded')
            print(f'{fitHitsArr[1][1]}<{A[1]} or {B[1]}<{fitHitsArr[1][1]}')

    #Last plane horizontal:
    geo.modules['Scifi'].GetSiPMPosition(5000000,A,B)
    if fitHitsArr[8][0]<B[0] or A[0]<fitHitsArr[8][0]:
        isInside = False
        if verbose:
            print('last x border exceeded')
            print(f'{fitHitsArr[8][0]}<{B[0]} or {A[0]}<{fitHitsArr[8][0]}')

    #Last plane vertical:
    geo.modules['Scifi'].GetSiPMPosition(5100000,A,B)
    if fitHitsArr[9][1]<A[1] or B[1]<fitHitsArr[9][1]:
        isInside = False
        if verbose:
            print('last y border exceeded')
            print(f'{fitHitsArr[9][1]}<{A[1]} or {B[1]}<{fitHitsArr[9][1]}')

    return isInside


def zPlaneArr(eventTree,geo):
    '''
    Return array with the 10 z coordinates of the planes.
    format: [1hor, 1ver, 2hor, ...,5hor, 5ver]
    /!\ THIS ASSUME Z FIX FOR ALL THE PLANE!
    geo: SndlhcGeo.GeoInterface(<geometry file>) for planes positions.
    eventTree: event list to get one event crossing all planes 
    to get each plane position.
    '''
    zArr = [0,0,0,0,0,0,0,0,0,0] # Fill Z coordinate of the planes
    # Warning message if z values > epsilon for different SiPM of same plane
    epsilon = 0.0001 

    for sTree in eventTree: # Need the first event with all planes hit:
        if goodEvent(eventTree = eventTree, nStations = 5, allowMore = False):
            A,B = ROOT.TVector3(),ROOT.TVector3()
            for hit in sTree.Digi_ScifiHits:
                hitID = hit.GetChannelID()
                geo.modules['Scifi'].GetSiPMPosition(hitID,A,B)
                indexArr = (hitID // 1000000-1) * 2  + (hitID // 100000) % 2
                #           add 2 for each plane        add 1 if vertical

                zVal = 0.5 * (A[2] + B[2])
                #Will overwrite if event of same plane, good to check if same z.
                if zArr[indexArr]!=0 and zArr[indexArr]-zVal>epsilon:
                    print(f'WARNING: SciFi planes {indexArr} not aligned in Z direction!')
                zArr[indexArr] = zVal # Fill z array
            break
    return zArr

def sortHitStation(clusterArr,stationArr):
    '''
    return the array of clusters from clusterArr that belong
    to one of the selected stations of stationArr
    stationArr: Array with selected stations number, ex: [1,3,4,5]
    '''
    
    clusFit = []
    clusTest = []
    for cluster in clusterArr:
        # If the cluster station is inside our list to fit:
        if (cluster.GetFirst()//1000000) in stationArr:
            clusFit.append(cluster)
        else:
            clusTest.append(cluster)
    return clusFit, clusTest


def distFit(fitHits, clusterArr, testStationNum):
    '''
    Takes the cluster list of an event (clusterArr), the fit
    (constructed with the 4 stations) hits on all the the plane
    fitHits which is a 10 elements array and the 5th test station
    (testStationNum) to compute the distance between the test station's
    hits and the predicted position by the fit on the same plane.
    Missing hits are set to 1000, easy to remove afterward.

    fitHits: array of TVector3 plane-fit intersection from extendHits()
    Return horizontal fit-hit difference, vertical fit-hit difference,
    horizontal position of the selected hit, vertical position of the selected hit.
    '''
    horIndex = 2 * (testStationNum - 1)
    verIndex = horIndex +1
    horClusters = [1000] # Large numer so it will be the minimum
    verClusters = [1000] # only if the array is empty
    #Fill horizontal and vertical cluster lists of the test plane:
    for cluster in clusterArr:
        if (cluster.GetFirst()//1000000) == testStationNum:
            A,B = ROOT.TVector3(),ROOT.TVector3()
            cluster.GetPosition(A,B)
            if ((cluster.GetFirst()//100000) % 2 ) == 0: # horizontal
                horClusters.append(A[1])
            else:
                verClusters.append(A[0])

    horFit = fitHits[horIndex][1]
    verFit = fitHits[verIndex][0]
    # We want the differece (no abs) of the closest point(in absolute diff)
    horDiffIndex = np.argmin([abs(x - horFit) for x in horClusters]) # Minimal distance
    verDiffIndex = np.argmin([abs(y - verFit) for y in verClusters])
    horPos = horClusters[horDiffIndex]
    verPos = verClusters[verDiffIndex]
    horDiff = horPos - horFit
    verDiff = verPos - verFit
    return horDiff, verDiff, horPos , verPos

def testClusterProblem(eventTree):
    ''' 
    Test for cluster separated by only one unactivated SiPM channel.
    Print the result for easily see the neighboring clusters.
    eventTree: Event list from root file.
    return: None
    '''
    print('##################')
    IDArr = []
    for hit in sTree.Digi_ScifiHits:
        detID = int(hit.GetDetectorID())
        IDArr.append(detID)
    IDArr.sort()
    prevElement=0
    for element in IDArr:
        diff = element - prevElement
        if not (diff ==1):
            if diff < 10:
                print(f'!!!---{element - prevElement}---!!!')
            else:
                print('---')
        print(element)
        prevElement = element
def customFitStatus(trackTask, FitStations):
    '''
    Do manually the trackTask.ExecuteTask() so it can use arbitrary number
    of stations for the fit.
    Return the fit and fit status object containing the fitted track.
    stationArr: Array with selected stations number for the fit, ex: [1,3,4,5]
    '''
    trackTask.clusters = trackTask.scifiCluster() # Build cluster list
    clusFit, clusTest = sortHitStation(clusterArr=trackTask.clusters,stationArr=FitStations)
    # Do the manual fit with the reduced number of stations:
    # Need manual dictionnary for the fitTrack(hitlist):
    clusDict = {}
    for x in clusFit:
        A,B = ROOT.TVector3(),ROOT.TVector3()
        clusterID = x.GetFirst()
        dictEntery = {clusterID: x}
        clusDict.update(dictEntery)
    fit = trackTask.fitTrack(hitlist=clusDict) # Type ROOT.genfit.Track
    fitStatus= fit.getFitStatus()
    trackTask.event.fittedTracks = [fit] # Array to keep sndsw format
    return fit, fitStatus