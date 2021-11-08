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


def valueToColor(value, cmap_name='nipy_spectral', vmin=-171, vmax=220):
    '''
    Colormap from z coordinate to distinguish between SciFi planes.
    v_min, v_max: color range value in mm.
    value: z plane value to set the color.
    '''
    
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
    print(fitHits)
    fig= plt.figure(figsize = (7.5, 6),dpi=500, tight_layout=True)
    ax = plt.axes(projection="3d")

    fitArr = extendedFitArr(trackTask=trackTask, fitHits=fitHits)
    
    # cm to mm conversion *10:
    arrPosStart = [10 * element for element in arrPosStart]
    arrPosStop = [10 * element for element in arrPosStop]
    fitArr = [10 * element for element in fitArr]
    fitHits = [10 * element for element in fitHits]

    for hitNumber in range(len(arrPosStart)):
        ax.plot(
            xs = [arrPosStart[hitNumber][0], arrPosStop[hitNumber][0]], 
            ys = [arrPosStart[hitNumber][1], arrPosStop[hitNumber][1]],
            zs = [arrPosStart[hitNumber][2], arrPosStop[hitNumber][2]],
            ls = '-',
            # RGB format to color different Scifi planes
            color = valueToColor(abs(arrPosStart[hitNumber][2])) )

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

    ax.set_xlabel('x [mm]')
    ax.set_ylabel('y [mm]')
    ax.set_zlabel('z [mm]')
    
    # Legend fields before adding custom entery.
    handles, labels = ax.get_legend_handles_labels()
    #handler_map = matplotlib.legend.Legend.get_legend_handler_map()
    handler_map = matplotlib.legend.Legend.get_default_handler_map()
    # Define custom legend for cluster hits.
    cmap = plt.cm.nipy_spectral_r
    cmap_handles = [Rectangle((0, 0), 1, 1)]
    handler_rainbow = dict(zip(cmap_handles, [HandlerColormap(cmap)]))
    label_rainbow = ['Channel clusters']

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
    plt.savefig('figures/3d.png')
    plt.show()
    plt.close()


def display2dTrack(arrPosStart, arrPosStop, trackTask, fitHits):
    '''
    x-z and y-z 2d projections of the 3d track.
    arrPosStart, arrPosStop: Array each containing one end of the event clusters.
    trackTask: SndlhcTracking.Tracking() object containing the fit infos.
    fitHits: Array of TVector3 of hits to display.
    '''
    verStart = []
    verStop = []
    horStart = []
    horStop = []
    
    epsilon = 0.00001
    for i in range(len(arrPosStart)):
        delta = arrPosStart[i] - arrPosStop[i]
        if delta[0] < epsilon: #/!\ Change the condition if geometry not aligned anymore.
            verStart.append(arrPosStart[i])
            verStop.append(arrPosStop[i])
        elif delta[1] < epsilon:
            horStart.append(arrPosStart[i])
            horStop.append(arrPosStop[i])
        else:
            print('ERROR: fiber neither vertical nor horizontal!')
            print('Check geometry alignment or change hor/ver conditions.')


    fitArr = extendedFitArr(trackTask = trackTask, fitHits = fitHits)

    fig, (ax1, ax2) = plt.subplots(
        2,
        figsize = (5,7),
        dpi = 500,
        tight_layout = True)

    ax1.grid(axis = 'y')
    ax2.grid(axis = 'y')

    # cm to mm *10 conversion:
    horStart = [10 * element for element in horStart]
    horStop = [10 * element for element in horStop]
    verStart = [10 * element for element in verStart]
    fitArr = [10 * element for element in fitArr]
    verStop = [10 * element for element in verStop]    

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
    ax1.set_xlabel('z [mm]')
    ax1.set_ylabel('x [mm]')
   
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
    ax2.set_xlabel('z [mm]')
    ax2.set_ylabel('y [mm]')

    plt.savefig('figures/2d.png')
    plt.show()
    plt.close()


def chi2Hist(chi2_nDfArr, stationNum=0):
    '''
    Chi2/nDOF histogram.
    chi2_nDfArr: Array filled with the value for each event.
    '''
    binsArr = np.linspace(0,40,400)
    fig, ax = plt.subplots(figsize=(8,6), dpi=500, tight_layout=True)
    ax.hist(chi2_nDfArr, bins=binsArr)
    ax.set_xlim(left=0.0,right=40)
    plt.xlabel('chi2/dof')
    plt.ylabel('Number of events')
    plt.savefig(f'figures/chi2Hist_{stationNum}.png')
    #plt.show()
    plt.close()

def planesHist(nPlanesHit):
    '''
    Historam of number of planes hit.
    nPlanesHit: Array filled with number of planes hit (0-10) for each event.
    '''
    fig, ax = plt.subplots(figsize=(5,5), dpi=500, tight_layout=True)
    binsArr = np.linspace(2,11,10)
    ax.hist(nPlanesHit,bins=binsArr)
    plt.xlabel('Number of planes hit')
    plt.ylabel('Number of events')
    plt.savefig(f'figures/nPlanesHist.png')
    plt.show()
    plt.close()


def gauss(x, A, x0, sigma):
    '''
    Gaussian function f(x).
    Magnitude A, offset x0 and width sigma.
    '''
    return  A * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2))


def diffHist(horDiffArr, verDiffArr, stationNum):
    '''
    Histogram of position difference between the fit on 4 planes and
    the hit on th 5th one. Gaussian fit.
    horDiffArr, verDiffArr: vertical/horizontal difference between hit and fit,
    each event gives one element of the array.
    stationNum: The test station where the difference hit-fit is measured.
    Return [x0_x, err_x0_x, sigma_x, err_sigma_x, x0_y, err_x0_y, sigma_y, err_sigma_y]
    of gaussian fit for the full-stations plot.
    '''

    # conversion cm in mm *10:
    horDiffArr = [10 * element for element in horDiffArr]
    verDiffArr = [10 * element for element in verDiffArr]

    fig, (ax1, ax2) = plt.subplots(2)
    ax1.grid(axis = 'y')
    ax2.grid(axis = 'y')
    fig.set_size_inches(8, 8)
    fig.suptitle(f'Difference between hit and fit, test station {stationNum}.',
                 fontsize='x-large',
                 fontweight='bold')
    binsArr = np.linspace(-2.5,2.5,500)
    hist1n, hist1Bins, hist1Patches = ax1.hist(verDiffArr, bins=binsArr)
    ax1.set_xlabel(r'$\Delta$ x [mm]')
    ax1.set_ylabel('Number of events')

    hist2n, hist2Bins, hist2Patches = ax2.hist(horDiffArr, bins=binsArr)
    ax2.set_xlabel(r'$\Delta$ y [mm]')
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
        label = (f'Gaussian fit: x0 = ({param1[1]:.2} ± {errX01:.2}) mm'
               + f'\n                     sigma = ({abs(param1[2]):.2} ± {errSigma1:.2}) mm'))
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
        label = (f'Gaussian fit: x0 = ({param2[1]:.2} ± {errX02:.2}) mm'
               + f'\n                     sigma = ({abs(param2[2]):.2} ± {errSigma2:.2}) mm'))
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
    fitArr: array of (return format from diffHist()) for each station
    '''

    # fitArr already in mm
    stationsArr = [1,2,3,4,5]
    zeroArr = [0,0,0,0,0]
    # Extract data from array: See order of return of diffHist()
    x0_x = [x[0] for x in fitArr]
    err_x0_x = [x[1] for x in fitArr]
    sigma_x = [x[2] for x in fitArr]
    err_sigma_x = [x[3] for x in fitArr]
    x0_y = [x[4] for x in fitArr]
    err_x0_y = [x[5] for x in fitArr]
    sigma_y = [x[6] for x in fitArr]
    err_sigma_y = [x[7] for x in fitArr]

    fig, ((ax1, ax3),(ax2, ax4)) = plt.subplots(
        nrows = 2,
        ncols = 2,
        sharex = 'col',
        tight_layout = True)
    fig.set_size_inches(8, 8)
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
    ax1.set_ylabel('X offset [mm]')
    ax2.set_ylabel(r'$\sigma_x$ [mm]')
    ax3.set_ylabel('Y offset [mm]')
    ax4.set_ylabel(r'$\sigma_y$ [mm]')
    ax2.set_xlabel('Test station')
    ax4.set_xlabel('Test station')
    ax1.set_ylim(bottom = -0.7, top = 0.7)
    ax3.set_ylim(bottom = -0.7, top = 0.7)
    ax2.set_ylim(bottom = 0.05, top = 0.35)
    ax4.set_ylim(bottom = 0.05, top = 0.35)
    ax1.grid()
    ax2.grid()
    ax3.grid()
    ax4.grid()
    plt.savefig('figures/FullStationsDiff.png')
    plt.show()
    plt.close()


def diffPosHist(posArr, diffArr, binsPos, labels, fileName, isCrossed):
    '''
    Difference fit-hit versus position on plane histogram.
    diffArr, posArr, binsPos must be given in cm, conversion in mm inside function.

    posArr: cluster position, one element of the array for each event.
    diffArr: difference fit-hit, one element of the array for each event.
    binsPos: histogram bins for the cluster position. Must set the limits.
    according to the plane geometry.
    labels: [xlabel, ylabel] for the axes.
    fileName: file name without extension, to save the plot.
    isCrossed: True if X-Y axis are mixed in the same plot to check rotations.
    Return: slope and its uncertainty for rotation global plot, else None.
    '''

    fig, ax = plt.subplots(figsize=(8,6),dpi=500,tight_layout=True)

    plt.rcParams.update({'font.size': 10})

    # cm to mm conversion *10:
    diffArr = [10 * element for element in diffArr]
    posArr = [10 * element for element in posArr]
    binsPos = [10 * element for element in binsPos]

    # Histogram limits:
    yMin = -2.5
    yMax = 2.5
    xMin = -400
    xMax = 400
    nBins = [60,30]
    if isCrossed: # Difference in y axis for better horizontal fit.
        binx = binsPos
        biny = np.linspace(-2.5,2.5,50)
        xHist = posArr
        yHist = diffArr
    else:
        binx = np.linspace(-2.5,2.5,50)
        biny = binsPos
        xHist = diffArr
        yHist = posArr
    
    plt.hist2d(
        xHist,
        yHist,
        bins = [binx,biny],
        cmap = plt.get_cmap('gist_heat_r'))

    cbar = plt.colorbar()
    cbar.ax.set_ylabel('Number of events')

    # Linear fit for rotation offset:
    if isCrossed:
        linFitModel, cov = np.polyfit(posArr, diffArr, 1,cov = True) # Compute slope and intercept
        linFitFunc = np.poly1d(linFitModel) # Make the line equation
        #print(linFitModel)
        ax.plot(
            binx,
            linFitFunc(binx),
            color = 'b',
            label = f'Slope: {np.format_float_scientific(linFitModel[0], precision = 2,)}±{np.format_float_scientific(np.sqrt(cov[0][0]), precision = 2,)}')
        plt.legend()
    plt.xlabel(labels[0])
    plt.ylabel(labels[1])
    
    plt.savefig(f'figures/{fileName}.png')
    #plt.show()
    plt.close()
    if isCrossed:
        return linFitModel[0], np.sqrt(cov[0][0]) # Return slope and its uncertainty.

def rotationAngle(xPosyOff_Slope, xPosyOff_SlopeErr, yPosxOff_Slope, yPosxOff_SlopeErr):
    '''
    Global rotation plot for each stations, 
    Use what is returned by diffPosHist().
    xPosyOff: Cross the X position and the Y offset
    yPosxOff: Cross the Y position and the X offset
    '''
    x1 = [0.95, 1.95, 2.95, 3.95, 4.95]
    x2 = [1.05, 2.05, 3.05, 4.05, 5.05]
    zeroes = [0, 0, 0, 0, 0]

    fig, ax = plt.subplots(figsize=(8,6),dpi=500,tight_layout=True)
    # Both plane slopes must be the same in absolute value.
    xPosyOff_Slope = [abs(element) for element in xPosyOff_Slope]
    yPosxOff_Slope = [abs(element) for element in yPosxOff_Slope]
    ax.errorbar(
        x = x1,
        y = xPosyOff_Slope,
        xerr = zeroes,
        yerr = xPosyOff_SlopeErr,
        ls = '',
        marker = 'x',
        markeredgecolor = 'b',
        label = 'X position, Y offset')

    ax.errorbar(
        x = x2,
        y = yPosxOff_Slope,
        xerr = zeroes,
        yerr = yPosxOff_SlopeErr,
        ls = '',
        marker = 'x',
        markeredgecolor = 'r',
        label = 'Y position, X offset')

    ax.grid(axis = 'y')
    plt.tick_params(bottom=False)
    plt.xlabel('Station number')
    plt.ylabel('Relative rotation angle [rad]')
    plt.legend()
    plt.savefig('figures/angles.png')
    plt.close()
    