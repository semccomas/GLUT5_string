#!/usr/bin/env python
# coding: utf-8

# The purpose of this script is to parse the output of `gmx hbond` and output a score into a small numpy array
# 
# This will eventually be converted to a python script, so that it can be incorporated into a master bash script where I run the hbond, and then immediately parse the results. The actual analysis and plotting will take place in `sugar_plot_hbond_score`


import numpy as np
import MDAnalysis as mda
from itertools import groupby
import sys
import pandas as pd



#histogram_square = sys.argv[1]
histogram_square = 'BFRU_transition_occ'

indir = f'../textfiles_out/gmx_hbond_files/{histogram_square}'
binding_site_res = {18:'Q160', 19:'N287'}
all_residues_res_descriptions = []

for residue in binding_site_res:
    hbond_ndx = open(f"{indir}/{residue}.BFRU.hbond.ndx").read()
    hbond_pixmap = open(f"{indir}/{residue}.BFRU.hbmap.xpm").read().splitlines()
    
    
    # Parsing the index and pixmap file. 
    # **index:** As each is shaped like [header] then we just take the final header group, and split line by line (since these will match out pixmap file) As there would be one line after the header, we skip the first newline.
    # 
    # **pixmap** Each header should be exactly the same size, as it would contain the same amount of info
    
    hbond_groups = hbond_ndx.split(']')[-1].splitlines()[1:]
    
    total_mean_res_time = []
    total_std_res_time = []
    total_n_res_time = []
    
    #don't bother including interactions where only 2 H bonds occured during the whole frame
    ## This cutoff is coupled to a cutoff of the mean interaction time, in the case that we have only 
    ## one interaction but that it takes up like tons of frames
    cutoff = 2 
    
    for interaction_group in hbond_pixmap[14:]:
        
        groups = groupby(interaction_group)
        continuous_counts = [(label, sum(1 for _ in group)) for label, group in groups]
       
        individual_bond_res_time = []
        for line in continuous_counts:
            if 'o' in line[0]:
                individual_bond_res_time.append(line[1])
        
        if len(individual_bond_res_time) > cutoff or np.mean(individual_bond_res_time) > 1:
            total_mean_res_time.append(np.mean(individual_bond_res_time))
            total_std_res_time.append(np.std(individual_bond_res_time))
            total_n_res_time.append(len(individual_bond_res_time))
        
        num_frames = (interaction_group.count(' ') + interaction_group.count('o'))
    
    total_mean_res_time = np.mean(total_mean_res_time)
    total_std_res_time = np.std(total_std_res_time)
    total_n_res_time = sum(total_n_res_time)
    
    
    res_description = {"mean_hbond_frames":total_mean_res_time,
                      "std_hbond_frames":total_std_res_time,
                      "n_interactions_total":total_n_res_time,
                      "n_frames":num_frames,
                      "resid":binding_site_res[residue]}
    
    all_residues_res_descriptions.append(res_description)


final_out = pd.DataFrame(all_residues_res_descriptions)
final_out = final_out.set_index('resid')
final_out.to_pickle(f'{indir}/hbond_by_res.pkl')
