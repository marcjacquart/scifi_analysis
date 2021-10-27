
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
    '''Extend the fit position to include missing planes hit'''
    
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
    zArr = [0,0,0,0,0,0,0,0,0,0] # Fill Z coordinate of the planes
    # format: [1hor, 1ver, 2hor, ...,5hor, 5ver]
    # /!\ THIS ASSUME Z FIX FOR ALL THE PLANE!
    # Warning message if z values > epsilon for different SiPM of same plane
    epsilon = 0.0001 

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
    return zArr

def sortHitStation(clusterArr,stationArr):
    '''
    return the array of hits from hitArr that belongs
    to one of the stations from stationArr
    '''
    
    clusFit = []
    clusTest = []
    for cluster in clusterArr:
        # If the cluster station is inside our list to fit:
        if (cluster.GetFirst()//1000000) in stationArr:
            # print(f'{cluster.GetFirst()//1000000} is OK for fit!')
            clusFit.append(cluster)
        else:
            clusTest.append(cluster)
            # print(f'{cluster.GetFirst()//1000000} is NOT ok for fit! Used to test')
    return clusFit, clusTest


def distFit(fitHits, clusterArr, testStationNum):
    '''
    Takes the cluster list of an event (clusterArr), the fit
    (constructed with the 4 stations) hits on all the the plane
    fitHits which is a 10 elements array and the 5th test station
    (testStationNum) to compute the distance between the test station's
    hits and the predicted position by the fit on the same plane.
    Missing hits are set to 1000, easy to remove afterward.
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
            if ((cluster.GetFirst()//100000) % 2 ) == 0:
                #print(f'Test must be 0: {A[1]-B[1]}')
                horClusters.append(A[1])
            else:
                verClusters.append(A[0])
                #print(f'Test must be 0: {A[0]-B[0]}')

    horFit = fitHits[horIndex][1]
    verFit = fitHits[verIndex][0]
    # We want the differece (no abs) of the closest point(in absolute diff)
    horDiffIndex = np.argmin([abs(x - horFit) for x in horClusters]) # Minimal distance
    verDiffIndex = np.argmin([abs(y - verFit) for y in verClusters])
    horDiff = horClusters[horDiffIndex] - horFit
    verDiff = verClusters[verDiffIndex] - verFit
    return horDiff, verDiff