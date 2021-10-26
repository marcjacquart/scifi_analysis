import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

# Custom legend:
from matplotlib.patches import Rectangle
from matplotlib.legend_handler import HandlerBase

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


def valueToColor(value, cmap_name='nipy_spectral', vmin=-17.1, vmax=22):
    '''Colormap from z coordinate to distinguish between SciFi planes.'''
    norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
    cmap = matplotlib.cm.get_cmap(cmap_name)
    rgb = cmap(norm(abs(value)))[:3]  # rgba[:3] -> rgb
    color = matplotlib.colors.rgb2hex(rgb)
    return color


def extendedFitArr(trackTask, fitHits):
    '''
    Put fit intersection with planes in an array.
    Add by linear extrapolation the missing hits.
    return the full array of 10 ROOT.TVector3().
    '''

    # Fit infos
    fitArr = []
    # trackTask.event.fittedTracks are filled with the fitTrack() function
    for aTrack in trackTask.event.fittedTracks:
        for i in range(aTrack.getNumPointsWithMeasurement()):
            state = aTrack.getFittedState(i)
            pos = state.getPos()
            fitArr.append(pos)

    # Extend the fit display to missed hits points:
    for pos in fitHits:
        fitArr.append(pos) 
    return fitArr


def display3dTrack(arrPosStart, arrPosStop, trackTask, fitHits):
    '''
    arrPosStart/stop: position of the activated fibers, A&B from the getPos() function
    trackTask is the fit object from SndlhcTracking.Tracking()
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

    fitArr = extendedFitArr(trackTask=trackTask, fitHits=fitHits)
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


def display2dTrack(arrPosStart, arrPosStop, trackTask, fitHits):
    verStart = []
    verStop = []
    horStart = []
    horStop = []

    #print('Verify cluster proxymity:')
    for element in arrPosStart:
        pass
        #print(f'Start: x:{element[0]}, y:{element[1]}, z:{element[2]}')
    #print('---------------------------')
    
    epsilon = 0.00001
    for i in range(len(arrPosStart)):
        delta = arrPosStart[i] - arrPosStop[i]
        #print(f'Delta: x:{delta[0]}, y:{delta[1]}, z:{delta[2]}')
        if delta[0] < epsilon: #/!\ Change the condition is geometry not aligned anymore.
            verStart.append(arrPosStart[i])
            verStop.append(arrPosStop[i])
        elif delta[1] < epsilon:
            horStart.append(arrPosStart[i])
            horStop.append(arrPosStop[i])
        else:
            print('ERROR: fiber neither vertical nor horizontal!')
            print('Check geometry alignment or change hor/ver conditions.')


    fitArr = extendedFitArr(trackTask=trackTask, fitHits=fitHits)

    fig, (ax1, ax2) = plt.subplots(2)
    fig.set_size_inches(8, 8)
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
        marker = '.',
        label = 'Clusters')
    # Add fit:
    # Can use sort() only because it is a straight line
    ax1.plot(
        [vect[2] for vect in fitArr],
        [vect[0] for vect in fitArr],
        color = 'r',
        label = 'Fit')
    ax1.set_xlabel('z [cm]')
    ax1.set_ylabel('x [cm]')
    ax1.legend()
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
    ax2.plot(
        [vect[2] for vect in fitArr],
        [vect[1] for vect in fitArr],
        color = 'r')
    ax2.set_xlabel('z [cm]')
    ax2.set_ylabel('y [cm]')
    plt.show()
    plt.close()


def chi2Hist(chi2_nDfArr):
    binsArr = np.linspace(0,5000,5000)
    fig, ax = plt.subplots(figsize=(6,8), dpi=300, tight_layout=True)
    ax.hist(chi2_nDfArr, bins=binsArr)
    ax.set_xlim(left=0.0,right=5000)
    plt.xlabel('chi2/dof')
    plt.ylabel('Number of events')
    plt.show()
    plt.close()

def planesHist(nPlanesHit):   
    fig, ax = plt.subplots(figsize=(6,8), dpi=300, tight_layout=True)
    ax.hist(nPlanesHit)
    plt.xlabel('Number of planes hit')
    plt.ylabel('Number of events')
    plt.show()
    plt.close()

def diffHist(horDiffArr, verDiffArr):
    fig, (ax1, ax2) = plt.subplots(2)
    fig.set_size_inches(8, 8)
    fig.suptitle(f'Difference between hit and fit',
                 fontsize='x-large',
                 fontweight='bold')
    binsArr = np.linspace(-0.25,0.25,500)
    ax1.hist(verDiffArr, bins=binsArr)
    ax1.set_xlabel(r'$\Delta$ x [cm]')
    ax1.set_ylabel('Number of events')

    ax2.hist(horDiffArr, bins=binsArr)
    ax2.set_xlabel(r'$\Delta$ y [cm]')
    ax2.set_ylabel('Number of events')
    
    plt.show()
    plt.close()