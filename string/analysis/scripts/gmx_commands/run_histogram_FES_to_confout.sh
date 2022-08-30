#!/bin/bash

## max # grids possible is 2916
## to make sure i don't do anything wrong, I am just going chunk by chunk for now..
## 0 - 33 done (nothing to do)
for grid_number in $(seq 0 2916); do
    echo $grid_number
    indir=../confout_files/FES_grids_confouts/influx_apo_gate_CV/histogram_$grid_number
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

