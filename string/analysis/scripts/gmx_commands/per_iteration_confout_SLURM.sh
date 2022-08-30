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



### the point of this is to combine all beads across one iteration to see how they progress. All beads in one iteration would be one full transition. ie. Inopen --> OutOpen. 
## Therefore, the purpose of this script is to compare overall transitions between iterations


trajext='../../../../../2021071200_GLUT5_string_influx_TMD/GLUT5_string/string/string_sims/TMD_initial_path'
trajdir='influx_apo_gate_CV'

#trajext='../../../../../2021051900_GLUT5_string_efflux_strings/GLUT5_string/string/string_sims'
#trajdir='efflux_demyst_CV'


iteration=0
max_beads=15
other_beads_max=$((max_beads-1))

## can make start and stop since they're fixed
gmx editconf -f $trajext/$trajdir/md/0/0/restrained/confout.gro -o ../confout_files/pdb_clips/$trajdir.$iteration.0.pdb -n $trajext/$trajdir/topology/index.ndx -ndef <<EOF 
0
EOF

gmx editconf -f $trajext/$trajdir/md/0/$max_beads/restrained/confout.gro -o ../confout_files/pdb_clips/$trajdir.$iteration.$max_beads.pdb -n $trajext/$trajdir/topology/index.ndx -ndef << EOF
0
EOF

for i in $(seq 1 $other_beads_max); do 
gmx editconf -f $trajext/$trajdir/md/$iteration/$i/restrained/confout.gro -o ../confout_files/pdb_clips/$trajdir.$iteration.$i.pdb -n $trajext/$trajdir/topology/index.ndx -ndef << EOF
0
EOF

done

## allow user to use n_beads variable and trajname variable while still going in pdb order
eval "cat ../confout_files/pdb_clips/$trajdir.$iteration.{0..$max_beads}.pdb > ../confout_files/measure_per_iteration/$trajdir.$iteration.string.pdb"



