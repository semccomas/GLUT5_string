# GLUT5 string simulations project - July 2020 - September 2022



# Pipeline:
## 1. GLUT5 atomistic simulations 

### Running atomistic
First, simulate each state in conventional MD
Goal is to see sampling of the 5 states in the CV space, and measure that they are stable
1. [plot_gate_dists_jupyter.ipynb](https://github.com/semccomas/GLUT5_string/blob/master/GLUT5_atomistic/analysis/scripts/plot_gate_dists_jupyter.ipynb)
   - Plot gate distances as a kernel density estimation. Add a small point for each starting configuration to see drift
   - Found in Figure 2B
2. [RMSD_to_model.ipynb](https://github.com/semccomas/GLUT5_string/blob/master/GLUT5_atomistic/analysis/scripts/RMSD_to_model.ipynb)
   - Calculate RMSD of each simulation to the respective model that they started at
   - Additional RMSD calculations for relevant features, like ICH5, or protein with no loops
   - Found in Figure S2 A, B, C


### Auxilary scripts:
* [analyze_end_states.ipynb](https://github.com/semccomas/GLUT5_string/blob/master/GLUT5_atomistic/analysis/scripts/analyze_end_states.ipynb)
   - In the end I ran 3 more replicas of the outward open and inward open states to ensure that I had sampled the end states correctly. Here, I just see if any changes occur in other control replicas.








## 2. Targeted MD

### Running steered/ targeted MD
I used PLUMED for TMD. I wanted to interpolate state-by-state but occasionally this required some extra pushing, and took time to monitor. I developed a simple way to increase the force with every script submitted. Scripts submitted were short and sent to Beskow, usually 30-40m jobs.
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
* [steered_analysis](https://github.com/semccomas/GLUT5_string/blob/master/steered/analysis/scripts/steered_analysis.ipynb)
   - Older script where I was intially comparing features of the steering/ TMD that I was doing. Not really maintained
* [feature_selection](https://github.com/semccomas/GLUT5_string/blob/master/steered/analysis/scripts/feature_selection.ipynb)
   - Older script from when I was using demystefying to find CVs. Here I analyzed the CVs dictated important and clustered based on residue groups




## 3. String simulations with swarms of trajectories
### Setting up strings (after beads have been chosen). Scripts are in the [example folder](https://github.com/semccomas/GLUT5_string/tree/master/string/string_sims/TMD_initial_path/influx_BFRU_gate_CV)
1. Use [steered_confout_dump](https://github.com/semccomas/GLUT5_string/blob/master/string/string_sims/TMD_initial_path/influx_BFRU_gate_CV/steered_confout_dump.sh) to get frames from TMD
2. Run [input_maker](https://github.com/semccomas/GLUT5_string/blob/master/string/string_sims/TMD_initial_path/influx_BFRU_gate_CV/input_maker.ipynb)
   - This code has been modified from the [Delemotte lab GitHub] (https://github.com/delemottelab/string-method-swarms-trajectories/tree/master/examples/start-up)
   - You should already have the `md/0/*` `confout.gro` files from `steered_confout_dump` before running this
   - Follow startup instructions on Delemotte lab github for more complete info
3. Check initial string with [analyze_initial_string](https://github.com/semccomas/GLUT5_string/blob/master/string/string_sims/TMD_initial_path/influx_BFRU_gate_CV/analyze_initial_string.ipynb)
4. Submit string simulations to Beskow with [slurm_string_beskow](https://github.com/semccomas/GLUT5_string/blob/master/string/string_sims/TMD_initial_path/influx_BFRU_gate_CV/slurm_string_beskow.ipynb)
   - Using the method on Beskow requires the right environment. Again, info can be found on Delemotte lab github as linked above. Environment to set is [environment_update_string_method](https://github.com/semccomas/GLUT5_string/blob/master/string/analysis/scripts/conda_env/environment.yml)

### Measuring string progress
1. This is primarily done with [string_analysis](https://github.com/semccomas/GLUT5_string/blob/master/string/analysis/scripts/string_analysis.ipynb)
   - This will check string progress along CV's, and can give a rough FES using detailed balance method
2. Can also measure features along bead or strings (to see drift from starting points)
   - Per bead analysis = look at progress of one bead through time
      - Use: [per_bead_confout](https://github.com/semccomas/GLUT5_string/blob/master/string/analysis/scripts/gmx_commands/per_bead_confout_SLURM.sh) to get frames. Analyze with [per_bead_analysis](https://github.com/semccomas/GLUT5_string/blob/master/string/analysis/scripts/per_bead_analysis.ipynb)
   - Per iteration analysis = look at progress of several entire strings along a few iterations (like max 10)
      - Use: [per_iteration_confout](https://github.com/semccomas/GLUT5_string/blob/master/string/analysis/scripts/gmx_commands/per_iteration_confout_SLURM.sh) to get frames. Analyze with [per_iteration_analysis](https://github.com/semccomas/GLUT5_string/blob/master/string/analysis/scripts/per_iteration_analysis.ipynb)

### Free energy surfaces
As you start to get better sampling, you can begin to look at the free energy of your system
1. [string_MSM_analysis](https://github.com/semccomas/GLUT5_string/blob/master/string/analysis/scripts/string_MSM_analysis_TICs_deep_time-dev.ipynb)
   - Also developed by Delemotte lab, original script can be found on their GitHub. This can also calculate a 1D FES. This is the FES method calculation used in the paper (both 2D and 1D). 
   - This uses [string_tica_msm](https://github.com/semccomas/GLUT5_string/blob/master/string/analysis/scripts/string_tica_msm.py) as a function to call


## 4. Analysis of strings and features for paper
### Histogram grid extraction
1. [histogram_FES_to_confout](https://github.com/semccomas/GLUT5_string/blob/master/string/analysis/scripts/histogram_FES_to_confout.py)
   - I run this on a cluster to speed up. It will **take time** to run this. For ~500 iterations, 16 beads, this took about 4 days. This will parse whatever piece of the histogram you tell it to.
   - In this, it will call [FES_grids_confouts_selection](https://github.com/semccomas/GLUT5_string/blob/master/string/analysis/scripts/FES_grids_confouts_selection.sh) to get the appropriate confout files from relevant histogram boxes, and collate them into a trajectory
   - To automate this, I run [run_histogram_FES_to_confout](https://github.com/semccomas/GLUT5_string/blob/master/string/analysis/scripts/gmx_commands/run_histogram_FES_to_confout.sh). This will just run in order, 1-?? histogram boxes. A lot of this is a bit fragile still and requires a lot of tinkering to work
2. If you need to align your trajectory on itself (for clustering of sugar poses t.ex.): first run [run_trjconv_for_cluster_on_histogram](https://github.com/semccomas/GLUT5_string/blob/master/string/analysis/scripts/gmx_commands/run_trjconv_for_cluster_on_histograms.sh)

### Analysis of sugar position
1. First, sugar poses must be clustered. Use [run_cluster](https://github.com/semccomas/GLUT5_string/blob/master/string/analysis/scripts/gmx_commands/run_cluster.sh)
2. Then, sugar cluster positions are parsed with [sugar_coordination_analysis](https://github.com/semccomas/GLUT5_string/blob/master/string/analysis/scripts/sugar_coordination_analysis.ipynb)
   - This generates the figure seen in Fig 3B
3. If you have frames that you want to take from certain histogram bins, including finding bins nearest to a certain state (as defined by an input structure, but you could change that), run [sugar_snapshots_prep](https://github.com/semccomas/GLUT5_string/blob/master/string/analysis/scripts/sugar_snapshots_prep.ipynb). It essentially just takes the frames found in the top cluster and puts them into a pdb.
   - This generates the inserts seen in Fig 3B next to the map, and Fig 3C. 
4. If you have a certain histogram in mind and want to measure distances to protein residues, use [distances_for_sugar_coupling](https://github.com/semccomas/GLUT5_string/blob/master/string/analysis/scripts/distances_for_sugar_coupling.ipynb).
   - Very straightforward, just treating the histo grid like a normal trajectory.
   - This generates Fig 3D


### Movement of TM7b, TM10b, and salt bridge residues
1. [TM7b_TM10](https://github.com/semccomas/GLUT5_string/blob/master/string/analysis/scripts/TM7b_TM10.ipynb) will calculate (if specified) the TM7b angle, TM10b RMSD, and salt-bridge distances for each histogram. This will take time, so the calculations are saved to disc and then it's best to set the calculation prompt to `False`
   - Plots for Figure 4, Figure 5 B,C come from here



### Free energy surface plotting
* [plot_FES_for_paper](https://github.com/semccomas/GLUT5_string/blob/master/string/analysis/scripts/plot_FES_for_paper.ipynb) Lifts same code from string MSM FES above, just that the F and extent.npy files are saved to disc so that we don't have to calculate the TICA and MSM from scratch. Then, compare to known models.
   - Figure 2D and 2E come from here



### Useful functions used occasionally
* [plot_parameters_for_paper](https://github.com/semccomas/GLUT5_string/blob/master/string/analysis/scripts/plot_parameters_for_paper.py) for consistent figure making
* [gate_functions](https://github.com/semccomas/GLUT5_string/blob/master/string/analysis/scripts/gate_functions.py) for easy measurements of gate distances. A bit on the old side




## 5. Auxillary scripts
* It can be useful to measure the salt-bridge distances across iterations, such as in Figure S6. Use [influx_v_efflux_salt_bridge_for_paper](https://github.com/semccomas/GLUT5_string/blob/master/string/analysis/scripts/influx_v_efflux_salt_bridge_for_paper.ipynb) for this










