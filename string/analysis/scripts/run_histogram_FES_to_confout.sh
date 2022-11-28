#!/bin/bash

## CHECKLIST!

## 1. 
## create confout_files directory       - done H
## change indir name to this       -  done H

## 2. 
# create ../string_sims/rGLUT5_mutation/xxx 
# run string_MSM_deeptime_dev.ipynb to get cv_MSM and F_MSM
# change input names in histogram_FES_to_confout.py:
#      1 - MSM_sim_dir variable  - need this for finding F for plotting (just for your personal reference in test.png) 
#                and also for cv_proj, to find actual IC and EC values in each grid
#      2 - name_sim variable  - for naming in FES_grids_confouts_selection
#      3 - indir variable (lower down on the page) - for finding where the confout.gro files are stored


## max # grids possible is 2916
## to make sure i don't do anything wrong, I am just going chunk by chunk for now..
## 0 - 33 done (nothing to do)
for grid_number in $(seq 0 1); do
    echo $grid_number
    indir=../confout_files/FES_grids_confouts/H386F_BFRU/histogram_$grid_number
    if [ ! -d "$indir" ]
    then
        mkdir $indir
    else
       echo 'directory exists: ' $indir
    fi
    python histogram_FES_to_confout.py $grid_number
   

    if grep -q no "$indir/runjob.txt"
    then
	    echo "no runjob, skipping FES_grids_confouts_selection.sh"
    else
	    echo "no runjob fail, running confouts selection now"
            bash FES_grids_confouts_selection.sh
    fi

    rm $indir/*clu*pdb
done

