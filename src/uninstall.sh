#!/bin/bash

###### Pipeline information:
# Last update: V1.00, 04/19/2021
# Uninstall all dependencies to release space
# Only need to run once
######


### get paths for pipe
pipe_path="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )" 
repo=$(echo ${pipe_path%/src})



echo "All dependencies have been removed."
date
