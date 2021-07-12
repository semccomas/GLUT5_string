#!/bin/bash

# Submit to the tcb partition
#SBATCH -p tcb

# The name of the job in the queue
#SBATCH -J conf
# wall-clock time given to this job
#SBATCH -t 00:10:00

# Number of nodes and number of MPI processes per node
#SBATCH -N 1 --ntasks-per-node=4


# Output file names for stdout and stderr
#SBATCH -e error.err -o output.out

module load gromacs/2020.2

### the point of this is to combine all iterations for one bead to see how an individual bead is behaving over iterations. There should be some drift but no big changes!!! 
## Therefore, the purpose of this script is to compare overall stability of one bead to measure convergence and relate to FES

#trajdir='influx_BFRU_demyst_CV'
trajdir="efflux_gate_CV"
trajext="../../string_sims"
#trajext="~/Desktop/mnt_TCB"

bead=32  #choose 1-32 (0 and 33 are fixed...)
min_iteration=1
max_iteration=112
if [ $trajdir = "influx_BFRU_demyst_CV" ]; then
outoption=18
index=../confout_files/index_influx_BFRU.ndx

else
outoption=1
index=../confout_files/index_efflux.ndx

fi


for i in $(seq $min_iteration $max_iteration); do 
gmx editconf -f $trajext/$trajdir/md/$i/$bead/restrained/confout.gro -o ../confout_files/pdb_clips/$trajdir.bead_$bead.iteration_$i.pdb -n $index -ndef << EOF
$outoption
EOF

done

## allow user to use n_beads variable and trajname variable while still going in pdb order
eval "cat ../confout_files/pdb_clips/$trajdir.bead_$bead.iteration_{$min_iteration..$max_iteration}.pdb > ../confout_files/measure_per_bead/$trajdir.bead_$bead.string.pdb"



