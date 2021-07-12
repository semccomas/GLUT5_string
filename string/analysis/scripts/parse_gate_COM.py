#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  7 15:03:35 2021

@author: semccomas
"""


## read .xvg files parsed by gate_COM_xvg.sh
# files are structured like:    md/$iteration/$bead/$swarm/EC_gate_COM.xvg

import numpy as np

mdname = '../../string_sims/influx_BFRU_demyst_CV'

mddir = f'{mdname}/md'
min_iteration = 51
max_iteration = 186

n_iterations = (max_iteration - min_iteration) + 1
n_beads = 34
n_swarms = 32

EC_dist_arr = np.zeros((n_iterations * (n_beads - 2) * n_swarms, 2))
IC_dist_arr = np.zeros((n_iterations * (n_beads - 2) * n_swarms, 2))

position = 0
for iteration in np.arange(min_iteration, max_iteration + 1):  #want to include max_iterations in the final calc
    for bead in np.arange(1, n_beads - 1):  #34 beads, assuming fixed ends where bead 0 & 33 shouldn't be calculated
        for swarm in np.arange(0, n_swarms):
            #print(iteration, bead, swarm)
            with open(f"{mddir}/{iteration}/{bead}/s{swarm}/EC_gate_COM.xvg") as f:
               data = np.array([line.strip().split()[1] for line in f], float)
               EC_dist_arr[position] = data
               
            with open(f"{mddir}/{iteration}/{bead}/s{swarm}/IC_gate_COM.xvg") as f:
               data = np.array([line.strip().split()[1] for line in f], float)
               IC_dist_arr[position] = data
               
            position = position + 1
            print(position)

gate_coordinates = np.dstack((IC_dist_arr, EC_dist_arr))
                
outdir = f'{mdname}/postprocessing'
np.save(f'{outdir}/gate_coordinates.npy', gate_coordinates)
