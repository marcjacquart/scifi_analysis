import matplotlib.pyplot as plt
import uproot
import awkward as ak
import numpy as np
import os

scifi_folder = '/home/marc/Physique/EPFL/Master/scifi/'
# Tuto file:
#filename = scifi_folder + 'scifi-reconstruction-tutorial/reconstructed_data/run_000025.root'
# cern run:
filename = scifi_folder + 'data/run_000000/reconstructed/run_000000.root'
folder_analysis = scifi_folder + 'analysis/plots/'
file = uproot.open(filename)


n_sensors = 10 # 5 planes * 2 x-y
n_channels = 1536 # number of SiPM channels of each sensor

# Tab containing the 10 sensors:
sensor_tab = []
for sensor_index in range(n_sensors):
    sensor_tab.append(file[f'sensor_{sensor_index}/hits'].arrays())


def create_folder_if_missing(path):
    if not os.path.isdir(path):
        os.mkdir(path)
        print('New directory created: '+path)


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
    plt.imshow(heatmap.T, 
        extent=extent,interpolation='nearest', aspect='auto', 
        origin='lower',
        cmap=plt.get_cmap('gist_heat_r'))
    cbar = plt.colorbar()
    cbar.ax.set_ylabel(axis_labels[2])
    plt.xlabel(axis_labels[0])
    plt.ylabel(axis_labels[1])
    plt.savefig(f'{output_folder}{output_name}.pdf')
    plt.close()




    # ax.hist(
    #     x = attribute_tab,
    #     bins = bins,
    #     histtype = 'step')
    # plt.xlabel(axis_labels[0])
    # plt.ylabel(axis_labels[1])
    # plt.savefig(f'{output_folder}{output_name}.pdf')
    # plt.close()

if True:
    # Channel histograms:
    for n_sensor in range(10):
        single_hist(
            sensor = sensor_tab[n_sensor],
            attribute = 'channel',
            output_name = f'channel_hist_{n_sensor}',
            output_folder = folder_analysis+'run_000000/channel_hist/',
            bins = np.linspace(0,n_channels-1,n_channels),
            axis_labels = ['Channel number','Number of hits'])


if True:
    # Value hist
    for n_sensor in range(10):
        single_hist(
            sensor = sensor_tab[n_sensor],
            attribute = 'value',
            output_name = f'value_hist_{n_sensor}',
            output_folder = folder_analysis+'run_000000/value_hist/',
            bins = np.linspace(0,12,101),
            axis_labels = ['Value','Number of hits'])

if True:
    for n_sensor in range(10):
        double_hist(
            sensor = sensor_tab[n_sensor],
            attributes = ['channel', 'value'],
            output_name = f'channel_value_hist_{n_sensor}',
            output_folder = folder_analysis+'run_000000/channel_value_hist/',
            binsTab = [np.linspace(0,n_channels-1,n_channels),np.linspace(0,8,101)],
            axis_labels = ['Channel number','Value','Number of hits'])