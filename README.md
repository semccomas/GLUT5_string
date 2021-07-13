# GLUT5 string simulations project - July 2020 - August 2021



# Pipeline:
 - GLUT5 atomistic simulations 


# Code annotation
## GLUT5_atomistic:


## steered



## string:
### string analysis and free energy surface generation:

### `confout.gro` analysis:
These scripts are to analyze the actual coordinates of the system during the simulations. Here we can measure how each individual bead is drifing between iterations, or how the overall transition pathway is looking between iterations.
- Drift of individual beads between iterations:

- Overall transition pathway comparison:
   - [per_iteration_confout_SLURM](https://github.com/semccomas/GLUT5_string/blob/master/string/analysis/scripts/per_iteration_confout_SLURM.sh): Edit `trajext`, `trajdir`, `iteration`, and `max_bead` for each simulation. This will read the `confout.gro` file for each bead in the given iteration and then concatenate them into one pdb file. 
   - Submitting this script to TCB partition is not so necessary, it is a pretty low-level computation
