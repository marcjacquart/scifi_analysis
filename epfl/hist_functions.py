import matplotlib.pyplot as plt
import numpy as np
import awkward as ak
from general_functions import create_folder_if_missing

def single_hist(
        sensor, attribute, output_name, output_folder, bins, axis_labels):
    ''' Histogram of chanel hits for a given sensor'''

    create_folder_if_missing(path=output_folder)
    attribute_tab = ak.flatten(sensor[attribute])

    fig, ax = plt.subplots(figsize=(16,9),dpi=300)
    ax.hist(
        x = attribute_tab,
        bins = bins,
        histtype = 'step')
    plt.xlabel(axis_labels[0])
    plt.ylabel(axis_labels[1])
    plt.savefig(f'{output_folder}{output_name}.pdf')
    plt.close()




def double_hist(
        sensor, attributes, output_name, output_folder, binsTab, axis_labels):
    ''' 
    Histogram of chanel hits for a given sensor
    axis_labels = [xlabel, ylabel, cbar_label]
    '''

    create_folder_if_missing(path=output_folder)
    attribute_tab_x = ak.flatten(sensor[attributes[0]])
    attribute_tab_y = ak.flatten(sensor[attributes[1]])

    fig, ax = plt.subplots(figsize=(16,9),dpi=300)
    heatmap, xedges, yedges = np.histogram2d(
        x = attribute_tab_x,
        y = attribute_tab_y,
        bins = binsTab)
    extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
    plt.clf() # Clear previous plots in memory
    plt.imshow(
        heatmap.T,
        extent=extent,
        interpolation='nearest',
        aspect='auto',
        origin='lower',
        cmap=plt.get_cmap('gist_heat_r'))
    cbar = plt.colorbar()
    cbar.ax.set_ylabel(axis_labels[2])
    plt.xlabel(axis_labels[0])
    plt.ylabel(axis_labels[1])
    plt.savefig(f'{output_folder}{output_name}.pdf')
    plt.close()
