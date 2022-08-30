#!/bin/bash -l

#SBATCH -A xxxx-x-xx

#SBATCH -J free_OutOcc

#SBATCH --nodes 8 --ntasks-per-node=32

#time in hours
#SBATCH -t 00:30:00

#SBATCH --mail-user=email@email.com --mail-type=BEGIN,END
#SBATCH -e error.log
#SBATCH -o output.log

export OMP_NUM_THREADS=1
APRUN_OPTIONS="-n 64 -d 1 -cc none"
module swap PrgEnv-cray PrgEnv-gnu
module unload gromacs
module add gromacs/2019.5-plumed


## note that we dont use aprun anymore since beskow upgrade. If you want to use non parallel calcs like grompp you can just type them without srun in front and itll run like normal .sh script


##  gmx grompp -f ../topology/TMD.mdp -n ../topology/index.ndx -p ../topology/topol.top -c OutOpen-OutOcc.8500ps.gro -r OutOpen-OutOcc.8500ps.gro -o OutOcc-Occ.tpr


srun gmx_mpi mdrun -s OutOcc-Occ.tpr -deffnm OutOcc-Occ -cpi OutOcc-Occ.cpt -plumed plumed1.dat



