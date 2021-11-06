#!/bin/bash

# Reconstruction script to automate the sfr-convert/sfr-reconstruct command

# Take the run name as input. Format: run_XXXXXX (X:0-9)
run_name=$1

# Define paths:
initial_path=$PWD
folder_with_paramfiles='/home/marc/Physique/EPFL/Master/scifi/scifi-reconstruction-tutorial/'
data_folder="/home/marc/Physique/EPFL/Master/scifi/data/$run_name/"
raw_folder="${data_folder}raw/"
converted_folder="${data_folder}converted/"
reconstructed_folder="${data_folder}reconstructed/"

# Create subfolder:
mkdir $raw_folder
mkdir $converted_folder
mkdir $reconstructed_folder

# Move downloaded data in raw folder:
mv ${data_folder}boards.csv ${raw_folder}boards.csv
mv ${data_folder}channels.csv ${raw_folder}channels.csv
mv ${data_folder}data.root ${raw_folder}data.root
mv ${data_folder}fe.csv ${raw_folder}fe.csv
mv ${data_folder}qdc_cal.csv ${raw_folder}qdc_cal.csv
mv ${data_folder}tdc_cal.csv ${raw_folder}tdc_cal.csv


# Go into tutorial directory to use the parameter files:
cd $folder_with_paramfiles
# Convert raw data in more readale format:
echo "Converting the root file..."
sfr-convert -u all ${raw_folder} ${converted_folder}${run_name} --verbose

# Reconstruct tracks:
echo "Reconstructing the tracks..."
sfr-reconstruct -u 3d ${converted_folder}${run_name}.root ${reconstructed_folder}${run_name} --verbose

# Go back in initial folder:
cd $initial_path
