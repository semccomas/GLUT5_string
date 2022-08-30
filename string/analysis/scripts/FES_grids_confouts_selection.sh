#!/bin/bash
#SBATCH -p tcb
#SBATCH -J conf
#SBATCH -t 03:00:00
#SBATCH -N 1 --ntasks-per-node=4
#SBATCH -e error.err -o output.out


mkdir ../confout_files/FES_grids_confouts/influx_BFRU_gate_CV/32
echo "removing all old pdbs in ../confout_files/FES_grids_confouts/influx_BFRU_gate_CV/32"
rm ../confout_files/FES_grids_confouts/influx_BFRU_gate_CV/32/*pdb 

module unload gromacs
module load gromacs/2020.4

gmx editconf -f /mnt/cephfs/projects/2021091701_GLUT5_string_influx_BFRU_TMD/GLUT5_string/string/string_sims/TMD_initial_path/influx_BFRU_gate_CV/md/312/13/s13/confout.gro -o ../confout_files/FES_grids_confouts/influx_BFRU_gate_CV/32/312-13-s13.clu0.pdb
cp ../confout_files/FES_grids_confouts/influx_BFRU_gate_CV/32/312-13-s13.clu0.pdb ../confout_files/FES_grids_confouts/influx_BFRU_gate_CV/32/FES_grid_all.start.pdb
gmx editconf -f /mnt/cephfs/projects/2021091701_GLUT5_string_influx_BFRU_TMD/GLUT5_string/string/string_sims/TMD_initial_path/influx_BFRU_gate_CV/md/104/13/s21/confout.gro -o ../confout_files/FES_grids_confouts/influx_BFRU_gate_CV/32/104-13-s21.clu1.pdb
gmx editconf -f /mnt/cephfs/projects/2021091701_GLUT5_string_influx_BFRU_TMD/GLUT5_string/string/string_sims/TMD_initial_path/influx_BFRU_gate_CV/md/107/13/s26/confout.gro -o ../confout_files/FES_grids_confouts/influx_BFRU_gate_CV/32/107-13-s26.clu2.pdb
eval "cat ./../confout_files/FES_grids_confouts/influx_BFRU_gate_CV/32/*.clu{0..2}.pdb > ./../confout_files/FES_grids_confouts/influx_BFRU_gate_CV/32/FES_grid_all.pdb"

gmx trjconv -f ../confout_files/FES_grids_confouts/influx_BFRU_gate_CV/32/FES_grid_all.pdb -o ../confout_files/FES_grids_confouts/influx_BFRU_gate_CV/32/FES_grid_all.xtc

rm ../confout_files/FES_grids_confouts/influx_BFRU_gate_CV/32/FES_grid_all.pdb

gmx pairdist -f ../confout_files/FES_grids_confouts/influx_BFRU_gate_CV/32/FES_grid_all.xtc -s ../confout_files/tpr_files/influx_BFRU_gate_CV.wholesys.tpr -ref "com of resid 30 to 37" -sel "com of resid 289 to 295" -o ../confout_files/FES_grids_confouts/influx_BFRU_gate_CV/32/FES_grid_all.EC_gate.xvg -xvg none

gmx pairdist -f ../confout_files/FES_grids_confouts/influx_BFRU_gate_CV/32/FES_grid_all.xtc -s ../confout_files/tpr_files/influx_BFRU_gate_CV.wholesys.tpr -ref "com of resid 136 to 145" -sel "com of resid 386 to 394" -o ../confout_files/FES_grids_confouts/influx_BFRU_gate_CV/32/FES_grid_all.IC_gate.xvg -xvg none