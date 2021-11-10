#!/usr/bin/env python
# coding: utf-8

# The point of this script is to gather all the ways that we could make numpy arrays for FES calculations. Generally, these take quite a long time and only need to be run once.
# This is loosely based on `parse_gate_COM.py` and `gate_COM_xvg.sh` which I used in previous string sims to get the gate COM out of string sims using other CVs
# 
# **The shape of these CVs must be**: `(n_iterations*n_swarms_per_iter*n_beads, n_frames_per_iter, n_cvs)`
# In the first dimension, this is also organized in a certain way... Iteration by iteration, then bead by bead, then swarm by swarm
# ex: 20 iterations, 10 beads, 3 swarms/ bead = (600,2,2) for 2 cvs
# iteration1, bead1, swarm1
# iteration1, bead1, swarm2...
# iteration19, bead10, swarm3
# iteration20, bead1, swarm1,
# iteration20, bead1, swarm2
# 
# Generally, we have been using 16 beads, with 2 fixed ends, so the calculated CVs will be for **14 beads**, and we have also used **32 swarms**

## this is adapted from an ipython script but it doesn't really make sense to have it there because these calculations take a very long time. 
## so now the idea is that we just un-comment out what we want to run at this time and send this off in the background
### Then we will parse them for the FES calculations in prep_other_coords.ipynb

import MDAnalysis as mda
import MDAnalysis.analysis.distances as d
import numpy as np
import subprocess

########################################
# # Darko's contacts
########################################
contacts_cv1_groups = {
                    'OUTWARD_OPEN':[('174_HN','31_HA'), ('196_HB1','103_HG11'), ('186_HD2','112_HG1')\
                                    , ('372_O','348_HB2'), ('318_OD1','286_HG11')],
                     "OUTWARD_OCCLUDED":[('420_HA1','68_HG13'), ('290_HB1','30_HB2'), ('315_HA1','171_HG1')],
                     "INWARD_OCCLUDED":[('308_HE21','174_HA1'), ('315_HA1','171_HG1')]
                      }

contacts_cv2_groups = {
                'OUTWARD_OPEN':[('308_HB1','174_O'), ('381_O','160_HB1'), ('381_HB2','164_HG22'),\
                                ('128_O','71_CD1'), ('315_HA1','171_HG1')],
                 "OUTWARD_OCCLUDED":[('390_HD11','329_HG22'), ('345_HB2','284_HB2'), ('236_HG22','224_HN'), \
                                     ('318_OD1','286_HG23'), ('416_OD1','29_HG21'), ('383_HA2','337_O'), ('372_O','348_HB1')],
                 "INWARD_OCCLUDED":[('198_O','19_HD2'), ('239_O','231_HB1'), ('387_HB','337_HD21'), ('386_HD1','325_HN'),\
                                    ('206_CD2','142_OH'), ('236_HG21','224_HG2'), ('122_HG21','67_OG1'), \
                                    ('390_HD22','329_HG21'), ('372_HB1','348_HB2')]
                        }


reffile = '../confout_files/tpr_files/influx_apo_gate_CV.wholesys.tpr'

start_iteration = 1
end_iteration = 745


def get_darko_arr(start_iteration, end_iteration, cv_group, u):
    n_iterations = end_iteration - start_iteration + 2
    distance_arr_cv1 = np.zeros((len(u.trajectory), len(contacts_cv1_groups[cv_group])))
    distance_arr_cv2 = np.zeros((len(u.trajectory), len(contacts_cv2_groups[cv_group])))

    ## CV1
    distance_index = 0
    for contact1,contact2 in contacts_cv1_groups[cv_group]:
        resid1,atomid1 = contact1.split('_')
        resid2,atomid2 = contact2.split('_')

        r1 = u.select_atoms(f"resid {int(resid1) - 1} and name {atomid1}") #MDA HAS 0 INDEXING, MINUS 1!!
        #print(r1, resid1, atomid1)
        r2 = u.select_atoms(f"resid {int(resid2) - 1} and name {atomid2}") #MDA HAS 0 INDEXING
        #print(r2, resid2, atomid2)
        for time_index, frame in enumerate(u.trajectory):
            distance_arr_cv1[time_index, distance_index] = d.dist(r1, r2)[2]

        distance_index = distance_index + 1



    ## CV2
    distance_index = 0
    for contact1,contact2 in contacts_cv2_groups[cv_group]:
        resid1,atomid1 = contact1.split('_')
        resid2,atomid2 = contact2.split('_')

        r1 = u.select_atoms(f"resid {int(resid1) - 1} and name {atomid1}") #MDA HAS 0 INDEXING, MINUS 1!!
        #print(r1, resid1, atomid1)
        r2 = u.select_atoms(f"resid {int(resid2) - 1} and name {atomid2}") #MDA HAS 0 INDEXING
        #print(r2, resid2, atomid2)


        for time_index, frame in enumerate(u.trajectory):
            distance_arr_cv2[time_index, distance_index] = d.dist(r1, r2)[2]

        distance_index = distance_index + 1

    
    return distance_arr_cv1, distance_arr_cv2


###### OUT OPEN
##OutOpen_big_distance_arr_cv1 = []
##OutOpen_big_distance_arr_cv2 = []
##
##for iteration in range(start_iteration, end_iteration + 1):
##    u = mda.Universe(reffile, f'../../string_sims/TMD_initial_path/influx_apo_gate_CV/md/{iteration}/{iteration}.all_beads_swarms.xtc')
##    
##    cv1, cv2 = get_darko_arr(start_iteration= start_iteration,
##                           end_iteration=end_iteration,
##                           cv_group = 'OUTWARD_OPEN',
##                           u = u)
##    OutOpen_big_distance_arr_cv1.append(cv1)
##    OutOpen_big_distance_arr_cv2.append(cv2)
##    
##    print(iteration)
##
##np.save('../textfiles_out/darko_dists_np/influx_apo_gate_CV/OutOpen_cv1.npy', OutOpen_big_distance_arr_cv1)
##np.save('../textfiles_out/darko_dists_np/influx_apo_gate_CV/OutOpen_cv2.npy', OutOpen_big_distance_arr_cv2)



### OUT OCC
OutOcc_big_distance_arr_cv1 = []
OutOcc_big_distance_arr_cv2 = []

for iteration in range(start_iteration, end_iteration + 1):
    u = mda.Universe(reffile, f'../../string_sims/TMD_initial_path/influx_apo_gate_CV/md/{iteration}/{iteration}.all_beads_swarms.xtc')
    
    cv1, cv2 = get_darko_arr(start_iteration= start_iteration,
                           end_iteration=end_iteration,
                           cv_group = 'OUTWARD_OCCLUDED',
                           u = u)
    OutOcc_big_distance_arr_cv1.append(cv1)
    OutOcc_big_distance_arr_cv2.append(cv2)
    
    print(iteration)


np.save('../textfiles_out/darko_dists_np/influx_apo_gate_CV/OutOcc_cv1.npy', OutOcc_big_distance_arr_cv1)
np.save('../textfiles_out/darko_dists_np/influx_apo_gate_CV/OutOcc_cv2.npy', OutOcc_big_distance_arr_cv2)




### IN OCC
InOcc_big_distance_arr_cv1 = []
InOcc_big_distance_arr_cv2 = []

for iteration in range(start_iteration, end_iteration + 1):
    u = mda.Universe(reffile, f'../../string_sims/TMD_initial_path/influx_apo_gate_CV/md/{iteration}/{iteration}.all_beads_swarms.xtc')
    
    cv1, cv2 = get_darko_arr(start_iteration= start_iteration,
                           end_iteration=end_iteration,
                           cv_group = 'INWARD_OCCLUDED',
                           u = u)
    InOcc_big_distance_arr_cv1.append(cv1)
    InOcc_big_distance_arr_cv2.append(cv2)
    
    print(iteration)


np.save('../textfiles_out/darko_dists_np/influx_apo_gate_CV/InOcc_cv1.npy', InOcc_big_distance_arr_cv1)
np.save('../textfiles_out/darko_dists_np/influx_apo_gate_CV/InOcc_cv2.npy', InOcc_big_distance_arr_cv2)





