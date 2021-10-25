import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import ROOT

# Custom legend:
from matplotlib.patches import Rectangle
from matplotlib.legend_handler import HandlerBase

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


def valueToColor(value, cmap_name='nipy_spectral', vmin=-18, vmax=22):
    '''Colormap from z coordinate to distinguish between SciFi planes.'''
    norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
    cmap = matplotlib.cm.get_cmap(cmap_name)
    rgb = cmap(norm(abs(value)))[:3]  # rgba[:3] -> rgb
    color = matplotlib.colors.rgb2hex(rgb)
    return color


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




class HandlerColormap(HandlerBase):
    '''
    Custom handles for cluster legend,
    Adapted from stackoverflow.com/questions/55501860/
    '''
    def __init__(self, cmap, num_stripes=10, **kw):
        HandlerBase.__init__(self, **kw)
        self.cmap = cmap
        self.num_stripes = num_stripes
    def create_artists(self, legend, orig_handle, 
                       xdescent, ydescent, width, height, fontsize, trans):
        stripes = []
        for i in range(self.num_stripes):
            s = Rectangle([xdescent + i * width / self.num_stripes, ydescent], 
                          width / self.num_stripes, 
                          height, 
                          fc=self.cmap((2 * i + 1) / (2 * self.num_stripes)), 
                          transform=trans)
            stripes.append(s)
        return stripes


def display3dTrack(arrPosStart, arrPosStop, trackTask, offset, fitHits):
    '''
    arrPosStart/stop: position of the activated fibers, A&B from the getPos() function
    trackTask is the fit object from SndlhcTracking.Tracking()
    offset to have nicer [0,40] coordinates, unused for now
    fitHits to display individually the missed hits of the fit on the planes.

    Uses matplotlib to display a 3d plot of the track, the fit and missing hits.

    '''

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
    # trackTask.event.fittedTracks are filled with the fitTrack() function
    for aTrack in trackTask.event.fittedTracks:
        for i in range(aTrack.getNumPointsWithMeasurement()):
            state = aTrack.getFittedState(i)
            pos = state.getPos() + offset
            fitArr.append(pos)

    # Extend the fit display to missed hits points:
    for pos in fitHits:
        fitArr.append(pos)

    ax.plot(
        xs = [element[0] for element in fitArr], 
        ys = [element[1] for element in fitArr],
        zs = [element[2] for element in fitArr],
        color = 'k',
        label = 'Fit')

    ax.scatter3D(
        xs = [hit[0] for hit in fitHits], 
        ys = [hit[1] for hit in fitHits],
        zs = [hit[2] for hit in fitHits],
        color = 'b',
        marker = '^',
        label = 'Missed hits')

    ax.set_xlabel('x [cm]')
    ax.set_ylabel('y [cm]')
    ax.set_zlabel('z [cm]')
    
    # Legend fields before adding custom entery.
    handles, labels = ax.get_legend_handles_labels()
    #handler_map = matplotlib.legend.Legend.get_legend_handler_map()
    handler_map = matplotlib.legend.Legend.get_default_handler_map()
    # Define custom legend for cluster hits.
    cmap = plt.cm.nipy_spectral_r
    cmap_handles = [Rectangle((0, 0), 1, 1)]
    handler_rainbow = dict(zip(cmap_handles, [HandlerColormap(cmap)]))
    label_rainbow = ['Fiber clusters']

    # Append new entery to legend and display it.
    # handles.append(cmap_handles)
    # labels.append(label_rainbow)

    # handler_map.update(handler_rainbow)
    legend1 = plt.legend(
        loc = 'upper right',
        handles = handles, 
        labels = labels, 
        handler_map = handler_map)
    legend2 = plt.legend(
        loc = 'upper left',
        handles = cmap_handles, 
        labels = label_rainbow, 
        handler_map = handler_rainbow)
    plt.gca().add_artist(legend1)
    plt.gca().add_artist(legend2)
    plt.show()
    plt.close()


def display2dTrack(arrPosStart,arrPosStop):
    verStart = []
    verStop = []
    horStart = []
    horStop = []

    print('Verify cluster proxymity:')
    for element in arrPosStart:
        print(f'Start: x:{element[0]}, y:{element[1]}, z:{element[2]}')
    print('---------------------------')
    
    epsilon = 0.00001
    for i in range(len(arrPosStart)):
        delta = arrPosStart[i] - arrPosStop[i]
        print(f'Delta: x:{delta[0]}, y:{delta[1]}, z:{delta[2]}')
        if delta[0] < epsilon: #/!\ Change the condition is geometry not aligned anymore.
            verStart.append(arrPosStart[i])
            verStop.append(arrPosStop[i])
        elif delta[1] < epsilon:
            horStart.append(arrPosStart[i])
            horStop.append(arrPosStop[i])
        else:
            print('ERROR: fiber neither vertical nor horizontal!')
            print('Check geometry alignment or change hor/ver conditions.')

    fig, (ax1, ax2) = plt.subplots(2)
    fig.set_size_inches(12, 8)
    fig.suptitle(f'Track 2d projections',
                 fontsize='x-large',
                 fontweight='bold')

    # z-x plane:
    # Horizontal fibers are lines in this plane
    for i in range(len(horStart)):
        ax1.vlines(
            x= horStart[i][2],
            ymin = min(horStart[i][0],horStop[i][0]),
            ymax = max(horStart[i][0],horStop[i][0]),
            colors = 'b')

    # Vertical lines are only points in this plane
    ax1.scatter(
        x=[point[2] for point in verStart],
        y=[point[0] for point in verStart],
        color = 'b',
        marker = '.')

    ax1.set_xlabel('z [cm]')
    ax1.set_ylabel('x [cm]')

    # y-z plane:
    # Vertical fibers are lines in this plane
    for i in range(len(verStart)):
        ax2.vlines(
            x= verStart[i][2],
            ymin = min(verStart[i][1],verStop[i][1]),
            ymax = max(verStart[i][1],verStop[i][1]),
            colors = 'b')
    # Horizontal lines are only points in this plane
    ax2.scatter(
        x=[point[2] for point in horStart],
        y=[point[1] for point in horStart],
        color = 'b',
        marker = '.')

    ax2.set_xlabel('z [cm]')
    ax2.set_ylabel('y [cm]')

    plt.show()
    plt.close()