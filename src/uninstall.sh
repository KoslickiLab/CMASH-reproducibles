#!/bin/bash

###### Pipeline information:
# v1.01, last update 07/23/2020. 
# Uninstall all dependencies to release space
# Only need to run once
######


### get paths for pipe
pipe_path="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )" 
repo=$(echo ${pipe_path%/src})


### remove conda path
conda env remove -p ${pipe_path}/Simulation_Env_py27
conda env remove -p ${pipe_path}/CMASH_Env_py37


### rm CAMISIM
rm -rf ${pipe_path}/CAMISIM

### rm NCBI annotation
rm NCBI_GenBank_*txt


echo "All dependencies have been removed."
date
