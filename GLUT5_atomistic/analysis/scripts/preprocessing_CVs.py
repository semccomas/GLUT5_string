#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 23 11:26:36 2020

@author: semccomas
"""

import mdtraj as md
import numpy as np
from sklearn.preprocessing import MinMaxScaler
scalar = MinMaxScaler()

### from old demystifying where I had all sims together
#traj = md.load_xtc('../input_f/protein_only/clipped/GLUT5_all.clipped.xtc', top='../input_f/protein_only/clipped/GLUT5_in.clipped.protein.start.pdb')

### now here is the new demystifying where I just have in and out open (9 
traj = md.load_xtc('../input_f/protein_only/clipped/GLUT5_in_GLUT5_out_ONLY.protein.clipped.xtc', top='../input_f/protein_only/clipped/GLUT5_in.clipped.protein.start.pdb')
top = traj.topology
contacts = md.compute_contacts(traj, scheme = 'ca')

## since indexing will be off with clipped trajectory, we need to relate the resids to their index number
## so we loop through these to make a dict, which we can reassign to the contacts output, so we have correct
## residue numbers
resid_dict = {}
for resid in top.residues:
    print(resid.resSeq, resid.index)    # == resid # 
    resid_dict[resid.index] = resid.resSeq

feature_to_resid = contacts[1]
feature_to_resid = np.vectorize(resid_dict.get)(feature_to_resid) #take the value for each of these keys

### again, all sims
#np.save('samples_features_arr/samples.clipped.ca.npy', contacts[0])
## for in out only
np.save('samples_features_arr_IN_OUT_ONLY/samples.clipped.ca.npy', contacts[0])

### again, all sims
#np.save('samples_features_arr/feature_to_resids.clipped.ca.npy', feature_to_resid)
## for in out only
np.save('samples_features_arr_IN_OUT_ONLY/feature_to_resids.clipped.ca.npy', feature_to_resid)


scalar.fit(contacts[0])
norm_contacts = scalar.transform(contacts[0])

### finally, all sims
#np.save('samples_features_arr/samples.clipped.ca.normalized.npy', norm_contacts)
### in out only
np.save('samples_features_arr_IN_OUT_ONLY/samples.clipped.ca.normalized.npy', norm_contacts)


'''

#### uncomment me for picking distances ONLY based on structures / models, not sims
traj = md.load_pdb('../input_f/protein_only/clipped/GLUT5_all.clipped.protein.start.pdb')
top = traj.topology
contacts = md.compute_contacts(traj, scheme = 'ca')

np.save('samples_features_arr/samples.MODELSONLY.clipped.ca.npy', contacts[0])

scalar.fit(contacts[0])
norm_contacts = scalar.transform(contacts[0])
np.save('samples_features_arr/samples.MODELSONLY.clipped.ca.normalized.npy', norm_contacts)
'''
