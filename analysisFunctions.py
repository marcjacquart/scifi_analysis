
import ROOT



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
