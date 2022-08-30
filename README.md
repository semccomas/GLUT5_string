# GLUT5 string simulations project - July 2020 - September 2022



# Pipeline:
## 1. GLUT5 atomistic simulations 
First, simulate each state in conventional MD
Goal is to see sampling of the 5 states in the CV space, and measure that they are stable
1. [plot_gate_dists_jupyter.ipynb](https://github.com/semccomas/GLUT5_string/blob/master/GLUT5_atomistic/analysis/scripts/plot_gate_dists_jupyter.ipynb)
   - Plot gate distances as a kernel density estimation. Add a small point for each starting configuration to see drift
   - Found in Figure 2B
2. [RMSD_to_model.ipynb](https://github.com/semccomas/GLUT5_string/blob/master/GLUT5_atomistic/analysis/scripts/RMSD_to_model.ipynb)
   - Calculate RMSD of each simulation to the respective model that they started at
   - Additional RMSD calculations for relevant features, like ICH5, or protein with no loops
   - Found in Figure S2 A, B, C


Auxilary scripts:
1. [analyze_end_states.ipynb](https://github.com/semccomas/GLUT5_string/blob/master/GLUT5_atomistic/analysis/scripts/analyze_end_states.ipynb)
   - In the end I ran 3 more replicas of the outward open and inward open states to ensure that I had sampled the end states correctly. Here, I just see if any changes occur in other control replicas.








## 2. Targeted MD

### Running steered/ targeted MD
I used PLUMED for TMD. I wanted to interpolate state-by-state but occasionally this required some extra pushing, and took time to monitor. I developed a simple way to increase the force with every script submitted. Scripts submitted were short and to Beskow, usually 30-40m.
1. [run_TMD.sh](https://github.com/semccomas/GLUT5_string/blob/master/steered/state_by_state_running/targeted_MD/plumed_master/run_TMD.sh)
   - Gmx run file for plumed on Beskow
2. [plumed1.dat](https://github.com/semccomas/GLUT5_string/blob/master/steered/state_by_state_running/targeted_MD/plumed_master/plumed1.dat)
   - Original plumed1.dat file to use in first run file
3. [plumed_continuous.sh](https://github.com/semccomas/GLUT5_string/blob/master/steered/state_by_state_running/targeted_MD/plumed_master/plumed_continuous.sh)
   - Script to automatically increase force used on next steps of TMD. Will read from plumed_template
4. [plumed_template.dat](https://github.com/semccomas/GLUT5_string/blob/master/steered/state_by_state_running/targeted_MD/plumed_master/plumed_template.dat)
   - Mock plumed file to update forces for plumed2.dat, plumed3.dat, etc. Must have same target input as plumed1.dat


### Analyzing steered / targeted MD
1. [skipping_occluded_for_paper](https://github.com/semccomas/GLUT5_string/blob/master/steered/analysis/scripts/skipping_occluded_for_paper.ipynb)
   - Code for comparing TMD gate projections with or without the occluded state
   - Found in Figure 2C and S3.
2. [TMD_to_strings.ipynb](https://github.com/semccomas/GLUT5_string/blob/master/steered/analysis/scripts/TMD_to_strings.ipynb)
   - Projecting all TMD frames in gate space, then selecting representative frames for string simulations
   - Corresponds to Figure 2D and S4A.


### Auxilary scripts:
1. [steered_analysis](https://github.com/semccomas/GLUT5_string/blob/master/steered/analysis/scripts/steered_analysis.ipynb)
   - Older script where I was intially comparing features of the steering/ TMD that I was doing. Not really maintained
2. [feature_selection](https://github.com/semccomas/GLUT5_string/blob/master/steered/analysis/scripts/feature_selection.ipynb)
   - Older script from when I was using demystefying to find CVs. Here I analyzed the CVs dictated important and clustered based on residue groups
