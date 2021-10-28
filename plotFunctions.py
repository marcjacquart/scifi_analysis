import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from scipy.optimize import curve_fit

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


def chi2Hist(chi2_nDfArr, stationNum=0):
    '''Chi2/nDOF histogram.'''
    binsArr = np.linspace(0,40,400)
    fig, ax = plt.subplots(figsize=(8,6), dpi=300, tight_layout=True)
    ax.hist(chi2_nDfArr, bins=binsArr)
    ax.set_xlim(left=0.0,right=40)
    plt.xlabel('chi2/dof')
    plt.ylabel('Number of events')
    plt.savefig(f'figures/chi2Hist_{stationNum}.png')
    #plt.show()
    plt.close()

def planesHist(nPlanesHit):
    '''Historam of number of planes hit.'''
    fig, ax = plt.subplots(figsize=(5,5), dpi=300, tight_layout=True)
    binsArr = np.linspace(2,11,10)
    ax.hist(nPlanesHit,bins=binsArr)
    plt.xlabel('Number of planes hit')
    plt.ylabel('Number of events')
    plt.savefig(f'figures/nPlanesHist.png')
    plt.show()
    plt.close()


def gauss(x, A, x0, sigma):
    ''' Gaussian function.'''
    return  A * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2))

def diffHist(horDiffArr, verDiffArr, stationNum):
    '''
    Histogram of position difference between the fit on 4 planes and
    the hit on th 5th one. Gaussian fit.
    Return [x0, err_x0, sigma, err_sigma] of gaussian fit.
    '''

    fig, (ax1, ax2) = plt.subplots(2)
    fig.set_size_inches(8, 8)
    fig.suptitle(f'Difference between hit and fit, test station {stationNum}.',
                 fontsize='x-large',
                 fontweight='bold')
    binsArr = np.linspace(-0.25,0.25,500)
    hist1n, hist1Bins, hist1Patches = ax1.hist(verDiffArr, bins=binsArr)
    ax1.set_xlabel(r'$\Delta$ x [cm]')
    ax1.set_ylabel('Number of events')

    hist2n, hist2Bins, hist2Patches = ax2.hist(horDiffArr, bins=binsArr)
    ax2.set_xlabel(r'$\Delta$ y [cm]')
    ax2.set_ylabel('Number of events')

    diff = hist1Bins[1] - hist1Bins[0]
    binCenters = [hist1Bins[i] + diff for i in range(len(hist1Bins)-1)]
    
    # Gaussian fits:
    param1, cov1 = curve_fit(f=gauss, xdata=binCenters, ydata = hist1n)
    errA1=np.sqrt(cov1[0][0])
    errX01=np.sqrt(cov1[1][1])
    errSigma1=np.sqrt(cov1[2][2])
    ax1.plot(
        binCenters, 
        [gauss(
            x = binCenters[i],
            A = param1[0],
            x0 = param1[1],
            sigma = param1[2]) 
            for i in range(len(binCenters))],
        color = 'r',
        label = (f'Gaussian fit: x0 = {param1[1]:.2} ± {errX01:.2}'
               + f'\n                     sigma = {abs(param1[2]):.2} ± {errSigma1:.2}'))
    ax1.legend()

    param2, cov2 = curve_fit(f=gauss, xdata=binCenters, ydata = hist2n)
    errA2 = np.sqrt(cov2[0][0])
    errX02 = np.sqrt(cov2[1][1])
    errSigma2 = np.sqrt(cov2[2][2])
    ax2.plot(
        binCenters, 
        [gauss(
            x = binCenters[i],
            A = param2[0],
            x0 = param2[1],
            sigma = param2[2]) 
            for i in range(len(binCenters))],
        color = 'r',
        label = (f'Gaussian fit: x0 = {param2[1]:.2} ± {errX02:.2}'
               + f'\n                     sigma = {abs(param2[2]):.2} ± {errSigma2:.2}'))
    ax2.legend()
    plt.savefig(f'figures/diffHistGauss_testStation{stationNum}.png')
    #plt.show()
    plt.close()
    # Sigma can be fitted negative, but return abs() by convention.
    resultFit = [param1[1], errX01, abs(param1[2]), errSigma1, param2[1], errX02, abs(param2[2]), errSigma2]
    return resultFit

def allPlanesGauss(fitArr):
    '''
    Global plot for the 5 test stations.
    '''
    stationsArr = [1,2,3,4,5]
    zeroArr = [0,0,0,0,0]
    # Extract data from array: See order of return of diffHist()
    x0_x = [x[0] for x in fitArr]
    err_x0_x = [10 * x[1] for x in fitArr]
    sigma_x = [x[2] for x in fitArr]
    err_sigma_x = [10 * x[3] for x in fitArr]
    x0_y = [x[4] for x in fitArr]
    err_x0_y = [10 * x[5] for x in fitArr]
    sigma_y = [x[6] for x in fitArr]
    err_sigma_y = [10 * x[7] for x in fitArr]

    fig, ((ax1, ax3),(ax2, ax4)) = plt.subplots(
        nrows = 2,
        ncols = 2,
        sharex = 'col',
        tight_layout = True)
    fig.set_size_inches(8, 8)
    # fig.suptitle(f'Difference between hit and fit',
    #              fontsize='x-large',
    #              fontweight='bold')
    ax1.errorbar(
        x = stationsArr,
        y = x0_x,
        xerr = zeroArr,
        yerr = err_x0_x,
        ls = '',
        marker = 'x',
        markeredgecolor = 'k')
    ax2.errorbar(
        x = stationsArr,
        y = sigma_x,
        xerr = zeroArr,
        yerr = err_sigma_x,
        ls = '',
        marker = 'x',
        markeredgecolor = 'k')
    ax3.errorbar(
        x = stationsArr,
        y = x0_y,
        xerr = zeroArr,
        yerr = err_x0_y,
        ls = '',
        marker = 'x',
        markeredgecolor = 'k')
    ax4.errorbar(
        x = stationsArr,
        y = sigma_y,
        xerr = zeroArr,
        yerr = err_sigma_y,
        ls = '',
        marker = 'x',
        markeredgecolor = 'k')
    ax1.set_ylabel('X offset [cm]')
    ax2.set_ylabel(r'$\sigma_x$ [cm]')
    ax3.set_ylabel('Y offset [cm]')
    ax4.set_ylabel(r'$\sigma_y$ [cm]')
    ax2.set_xlabel('Test station')
    ax4.set_xlabel('Test station')
    plt.savefig('figures/FullStationsDiff.png')
    plt.show()

    plt.close()

def diffPosHist(diffArr, posArr, isVetical, testStationNum):
    fig, ax = plt.subplots(figsize=(12,8),dpi=500)

    plt.rcParams.update({'font.size': 10})
    if isVetical:
        verHorStr = 'vertical'
        biny = np.linspace(-50,-5,90)
    else:
        verHorStr = 'horizontal'
        biny = np.linspace(15,60,90)
    # fig.suptitle(f'Residual and position,{verHorStr} test station {testStationNum}.',
    #          fontsize='large',
    #          fontweight='bold')
    xMin = -0.25
    xMax = 0.25
    yMin = -40
    yMax = 40
    nBins = [60,30]
    binx = np.linspace(-0.25,0.25,50)
    
    
    plt.hist2d(
        diffArr,
        posArr,
        bins = [binx,biny],
        #range = [[xMin,xMax],[yMin,yMax]],
        cmap = plt.get_cmap('gist_heat_r'))

    cbar = plt.colorbar()
    cbar.ax.set_ylabel('Number of events')
    if isVetical:
        plt.xlabel('X offset [cm]')
        plt.ylabel('X cluster position [cm]')
    else:
        plt.xlabel('Y offset [cm]')
        plt.ylabel('Y cluster position [cm]')
    plt.savefig(f'figures/diffPosHist_{testStationNum}_{verHorStr}.png')
    #plt.show()
    plt.close()