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


trajext='../../../../../2021051901_GLUT5_string_influx_strings/GLUT5_string/string/string_sims'
trajdir='influx_BFRU_demyst_CV'

#trajext='../../../../../2021051900_GLUT5_string_efflux_strings/GLUT5_string/string/string_sims'
#trajdir='efflux_demyst_CV'


iteration=25
if [ $trajdir = "influx_BFRU_demyst_CV" ]; then
outoption=18
index=../confout_files/index_influx_BFRU.ndx

else
outoption=1
index=../confout_files/index_efflux.ndx

fi

## can make start and stop since they're fixed
gmx editconf -f $trajext/$trajdir/md/0/0/restrained/confout.gro -o ../confout_files/pdb_clips/$trajdir.$iteration.0.pdb -n $index -ndef <<EOF 
$outoption
EOF

gmx editconf -f $trajext/$trajdir/md/0/33/restrained/confout.gro -o ../confout_files/pdb_clips/$trajdir.$iteration.33.pdb -n $index -ndef << EOF
$outoption
EOF

for i in $(seq 1 32); do 
gmx editconf -f $trajext/$trajdir/md/$iteration/$i/restrained/confout.gro -o ../confout_files/pdb_clips/$trajdir.$iteration.$i.pdb -n $index -ndef << EOF
$outoption
EOF

done

## allow user to use n_beads variable and trajname variable while still going in pdb order
eval "cat ../confout_files/pdb_clips/$trajdir.$iteration.{0..33}.pdb > ../confout_files/measure_per_iteration/$trajdir.$iteration.string.pdb"



