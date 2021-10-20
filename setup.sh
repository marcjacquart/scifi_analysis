#!/bin/bash

# Script to setup the environement for the snd software and updating it from github.

# Call "setup.sh first" for first installation (=initialize and download git reop),
# then only "setup.sh" for updating the repo and enterint the alice environnement. 
if [ "$1" == "first" ];then
    echo "*************************************************************"
    echo "First instalation: user/private/sndsw/ folder will be created"
    echo "The snd software git repository is then cloned from" 
    echo "https://github.com/SND-LHC/sndsw"
    echo "Then built and enter alienv to use the scripts."
    echo "*************************************************************"
    
    mkdir /afs/cern.ch/user/m/mjacqua2/private/sndsw/
    cd /afs/cern.ch/user/m/mjacqua2/private/sndsw/
    source /cvmfs/sndlhc.cern.ch/latest/setUp.sh
    git init
    git clone https://github.com/SND-LHC/sndsw
    aliBuild build sndsw -c $SNDDIST --always-prefer-system
    git pull --rebase

    echo "Alienv uses sndsw/latest-release folder."
    echo "Replace the setup.sh last command with output above if cannot enter alienv."
    #Replace the one below:
    alienv enter sndsw/latest-release

else
    cd /afs/cern.ch/user/m/mjacqua2/private/sndsw/
    source /cvmfs/sndlhc.cern.ch/latest/setUp.sh
    git pull https://github.com/SND-LHC/sndsw
fi
aliBuild build sndsw -c $SNDDIST --always-prefer-system

echo "Alienv uses sndsw/latest-release folder."
echo "Replace the setup.sh last command with output above if cannot enter alienv"
#Replace the one below:
alienv enter sndsw/latest-master-release

