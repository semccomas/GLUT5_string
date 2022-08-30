# GLUT5 string simulations project - July 2020 - September 2022



# Pipeline:
## 1. GLUT5 atomistic simulations 
First, simulate each state in conventional MD
Goal is to see sampling of the 5 states in the CV space, and measure that they are stable
1. [plot_gate_dists_jupyter.ipynb](https://github.com/semccomas/GLUT5_string/blob/master/GLUT5_atomistic/analysis/scripts/plot_gate_dists_jupyter.ipynb)
   - Plot gate distances as a kernel density estimation. Add a small point for each starting configuration to see drift
2. [RMSD_to_model.ipynb](https://github.com/semccomas/GLUT5_string/blob/master/GLUT5_atomistic/analysis/scripts/RMSD_to_model.ipynb)
   - Calculate RMSD of each simulation to the respective model that they started at
   - Additional RMSD calculations for relevant features, like ICH5, or protein with no loops

