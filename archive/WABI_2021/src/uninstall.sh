#!/bin/bash

###### Pipeline information:
# Last update: V1.00, 04/19/2021
# Uninstall all dependencies to release space
# Only need to run once
######


### get paths for pipe
pipe_path="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )" 
repo=$(echo ${pipe_path%/src})



### uninstall conda and git
cd ${pipe_path}
rm -r GenBank_db 2>/dev/null
rm -rf ./git_repo 2>/dev/null
rm -r conda_env 2>/dev/null



### delete running folders from WABI 2021 CMash manuscript
cd ${pipe_path}/../1_WABI_2021_CMash_manuscript/
ls -d CMash_out_* 2>/dev/null && rm -r CMash_out_*



echo "All dependencies have been removed."
date
