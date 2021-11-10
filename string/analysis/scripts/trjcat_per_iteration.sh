#!/bin/bash

# Submit to the tcb partition
#SBATCH -p tcb

# The name of the job in the queue
#SBATCH -J conf
# wall-clock time given to this job
#SBATCH -t 20:00:00
# Number of nodes and number of MPI processes per node
#SBATCH -N 1 --ntasks-per-node=4


# Output file names for stdout and stderr
#SBATCH -e error.err -o output.out

module load gromacs/2020.2

indir="/mnt/cephfs/projects/2021091701_GLUT5_string_influx_BFRU_TMD/GLUT5_string/string/string_sims/TMD_initial_path/influx_BFRU_gate_CV/md"



min_iterations=1
max_iterations=270
n_beads=14
n_swarms=31  ## is 32, but with original code it starts at s0

for iteration in $(seq $min_iterations $max_iterations); do
	eval "gmx trjcat -f $indir/$iteration/{1..$n_beads}/s{0..$n_swarms}/traj_comp.xtc -cat -o $indir/$iteration/$iteration.all_beads_swarms.xtc"

done
