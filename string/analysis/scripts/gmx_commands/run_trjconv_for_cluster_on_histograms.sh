#!/bin/bash
#SBATCH -p tcb
#SBATCH -J conf
#SBATCH -t 23:00:00
#SBATCH -N 1 --ntasks-per-node=4
#SBATCH -e error.err -o output.out

## hopefully this title is self-explanatory

module unload gromacs
module load gromacs/2020.2

tprdir=../confout_files/tpr_files

## done up to 1600
##for sim in ../confout_files/FES_grids_confouts/influx_BFRU_gate_CV/histogram*;
for i in $(seq 298 299); do
	sim=../confout_files/FES_grids_confouts/influx_BFRU_gate_CV/histogram_$i
	gmx trjconv -f $sim/FES_grid_all.xtc -o $sim/FES_grid_all.protonly.cluster.xtc -n $tprdir/index.ndx -s $tprdir/influx_BFRU_gate_CV.wholesys.tpr -pbc cluster << EOF
18
18
EOF

	gmx trjconv -f $sim/FES_grid_all.protonly.cluster.xtc -n $tprdir/index.ndx -s $sim/FES_grid_all.start.pdb -o $sim/FES_grid_all.protonly.cluster.fit.xtc -fit rot+trans << EOF
18
18
EOF
	gmx trjconv -f $sim/FES_grid_all.start.pdb -o $sim/FES_grid_all.start.protonly.pdb -s $sim/FES_grid_all.start.pdb -n $tprdir/index.ndx << EOF
18
EOF
done
