import uproot
import numpy as np
from hist_functions import single_hist, double_hist
from general_functions import create_folder_if_missing

scifi_folder = '/home/marc/Physique/EPFL/Master/scifi/'
# cern run:
run_name = 'run_000001'
filename = scifi_folder + f'data/{run_name}/reconstructed/{run_name}.root'
folder_analysis = f'{scifi_folder}analysis/plots/{run_name}'
file = uproot.open(filename)


n_sensors = 10 # 5 planes * 2 x-y
n_channels = 1536 # number of SiPM channels of each sensor

# Tab containing the 10 sensors:
sensor_tab = []
for sensor_index in range(n_sensors):
    sensor_tab.append(file[f'sensor_{sensor_index}/hits'].arrays())

# Create folder to put all the results
create_folder_if_missing(path=folder_analysis)

# Select plots to draw by calling the histogram functions:
# if True/False is an easy way to enable or not the multi-line function call.
if True:
    # Channel histograms:
    for n_sensor in range(10):
        single_hist(
            sensor = sensor_tab[n_sensor],
            attribute = 'channel',
            output_name = f'channel_hist_{n_sensor}',
            output_folder = f'{folder_analysis}/channel_hist/',
            bins = np.linspace(0,n_channels-1,n_channels),
            axis_labels = ['Channel number','Number of hits'])


if True:
    # Value hist
    for n_sensor in range(10):
        single_hist(
            sensor = sensor_tab[n_sensor],
            attribute = 'value',
            output_name = f'value_hist_{n_sensor}',
            output_folder = f'{folder_analysis}/value_hist/',
            bins = np.linspace(0,12,101),
            axis_labels = ['Value','Number of hits'])

if True:
    for n_sensor in range(10):
        double_hist(
            sensor = sensor_tab[n_sensor],
            attributes = ['channel', 'value'],
            output_name = f'channel_value_hist_{n_sensor}',
            output_folder = f'{folder_analysis}/channel_value_hist/',
            binsTab = [np.linspace(0,n_channels-1,n_channels),
                np.linspace(0,8,101)],
            axis_labels = ['Channel number','Value','Number of hits'])
        