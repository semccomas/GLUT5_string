#!/bin/bash

# Submit to the tcb partition
#SBATCH -p tcb

# The name of the job in the queue
#SBATCH -J conf
# wall-clock time given to this job
#SBATCH -t 00:10:00
#SBATCH -w cpu01
# Number of nodes and number of MPI processes per node
#SBATCH -N 1 --ntasks-per-node=4


# Output file names for stdout and stderr
#SBATCH -e error.err -o output.out

module load gromacs/2020.2

### the point of this is to combine all iterations for one bead to see how an individual bead is behaving over iterations. There should be some drift but no big changes!!! 
## Therefore, the purpose of this script is to compare overall stability of one bead to measure convergence and relate to FES

trajext='../../../../../2021071200_GLUT5_string_influx_TMD/GLUT5_string/string/string_sims/TMD_initial_path'
trajdir='influx_apo_gate_CV'

#trajext='../../../../../2021071201_GLUT5_string_efflux_TMD/GLUT5_string/string/string_sims/TMD_initial_path'
#trajdir='efflux_apo_gate_CV'



bead=12  #choose 1-14 (0 and 15 are fixed...)
min_iteration=1
max_iteration=300


for i in $(seq $min_iteration $max_iteration); do 
gmx editconf -f $trajext/$trajdir/md/$i/$bead/restrained/confout.gro -o ../confout_files/pdb_clips/$trajdir.bead_$bead.iteration_$i.pdb -n $trajext/$trajdir/topology/index.ndx -ndef << EOF
0
EOF

done

## allow user to use n_beads variable and trajname variable while still going in pdb order
eval "cat ../confout_files/pdb_clips/$trajdir.bead_$bead.iteration_{$min_iteration..$max_iteration}.pdb > ../confout_files/measure_per_bead/$trajdir.bead_$bead.string.pdb"



