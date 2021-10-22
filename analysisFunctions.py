import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d


def goodEvent(eventTree, nStations, allowMore):
    '''Return True if nStations (or more if allowMore) are hits in the event'''

    stations = {}
    for detectedHit in eventTree.Digi_ScifiHits:
        stations[detectedHit.GetDetectorID()//1000000] = 1
               
    for detectedHit in eventTree.Digi_MuFilterHit:
        plane = 100*(detectedHit.GetDetectorID()//1000)
        stations[plane] = 1
    if allowMore:
        if len(stations) >= nStations: return True
    else:
        if len(stations) == nStations: return True
    return False


def valueToColor(value, cmap_name='nipy_spectral', vmin=-18, vmax=22):
    '''Colormap from z coordinate to distinguish between SciFi planes.'''
    norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
    cmap = matplotlib.cm.get_cmap(cmap_name)
    rgb = cmap(norm(abs(value)))[:3]  # rgba[:3] -> rgb
    color = matplotlib.colors.rgb2hex(rgb)
    return color



def display3dTrack(arrPosStart, arrPosStop, trackTask, offset, fitHits):
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

    ax.plot(
        xs = [element[0] for element in fitArr], 
        ys = [element[1] for element in fitArr],
        zs = [element[2] for element in fitArr],
        color = 'r',
        label = 'fit')

    ax.scatter3D(
        xs = [hit[0] for hit in fitHits], 
        ys = [hit[1] for hit in fitHits],
        zs = [hit[2] for hit in fitHits],
        color = 'b',
        marker = '^',
        label = 'fit-plane intersection')

    ax.set_xlabel('x [cm]')
    ax.set_ylabel('y [cm]')
    ax.set_zlabel('z [cm]')
    plt.legend()
    plt.show()
    plt.close()

