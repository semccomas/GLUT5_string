#!/bin/bash -l

#SBATCH -A 2020-3-28

#SBATCH -J outocc-occ

#SBATCH --nodes 8 --ntasks-per-node=32

#time in hours
#SBATCH -t 00:30:00

#SBATCH --mail-user=sarah.mccomas@scilifelab.se --mail-type=BEGIN,END
#SBATCH -e error.log
#SBATCH -o output.log

export OMP_NUM_THREADS=1
APRUN_OPTIONS="-n 64 -d 1 -cc none"
module swap PrgEnv-cray PrgEnv-gnu
module unload gromacs
module add gromacs/2019.5-plumed


cntmax=20  #max number of plumed files to run
simname=OutOcc-Occ
## grompp we used to make this: (gmx 2019.5)
##gmx grompp -f TMD.mdp -o OutOpen-OutOcc.tpr -n index.ndx -p topol.top -c OUT.equilib.200ns.gro -r OUT.equilib.200ns.gro
###gmx grompp -f TMD.mdp -o OutOcc-Occ.tpr -n index.ndx -p topol.top -c OutOpen-OutOcc.17_5ns.gro -r OutOpen-OutOcc.17_5ns.gro

plumed_latest=$(find -name 'plumed?.dat' | sort -V | tail -1)
pcnt=${plumed_latest:8:1}  #take number from most recent plumed run
cnt=$((pcnt+1))   #add one, this is our new count


## first check if RMSD is reasonable yet, if less than 0.01nm then you make the cnt bigger than cntmax so you skip the rest of this script
rmsd=$(tail -n 1 COLVAR_$pcnt | tr -s ' ' | cut -d ' ' -f 3)  #this is the RMSD field
##rmsd=0.144748
cutoff=0.01
echo "Final RMSD from $plumed_latest run was $rmsd"
if (( $(echo "$rmsd < $cutoff" |bc -l) )); then
cnt=$((cntmax+1))
echo 'Done!!'
fi


## don't need a else statement because this will automatically break the loop below
while [ "$cnt" -le "$cntmax" ]
do

cp plumed_template.dat plumed$cnt.dat


### GET STEPS
### Process .log from trajectory file
###take "Step" info and line below it | make everything split by 1 white space | take the second field of this | take the last value

###step_num=$(grep 'Step' -A1 $simname.log | tr -s ' ' | cut -d ' ' -f 2 | tail -n 1)    # this was the old way, when I would take the final step but sometimes the checkpoint file hasn't been written in a while when these sims stops
step_num=$(grep 'Writing checkpoint' $simname.log | tr -s ' ' | cut -d ' ' -f 4 | tail -n 1)
echo $step_num
echo $cnt
step_num_start=$(expr $step_num + 5)
step_num_mid=$(expr $step_num + 5000) 
step_num_end=$(expr $step_num + 150000) #about 0.3ns in - for a 30m simulation this is about half of the wall time

sed -i "s/STEP_0/$step_num_start/" plumed$cnt.dat 
sed -i "s/STEP_1/$step_num_mid/" plumed$cnt.dat 
sed -i "s/STEP_2/$step_num_end/" plumed$cnt.dat 




### GET KAPPA
if [ "$cnt" -eq 1 ] ## need special params for 1
then
	kappa_mid=2500
	kappa_end=5000

# i know this is a bit counterintuative to define end first, but it will always be bigger than 2500 by 2 to the power of the count (eg. plumed2.dat should have 5000, 10000, plumed3 should have 10000 and 20000. 

else
        kappa_mid=$(grep -oP '(?<=KAPPA2=)\w+' plumed$pcnt.dat)   #start kappa where we left off
        kappa_end=$(expr 2 \* $kappa_mid) 
fi

sed -i "s/KAPPA_1/$kappa_mid/" plumed$cnt.dat  
sed -i "s/KAPPA_2/$kappa_end/" plumed$cnt.dat  


sed -i "s/COLVAR_X/COLVAR_$cnt/" plumed$cnt.dat   

#echo
srun gmx_mpi mdrun -s $simname.tpr -deffnm $simname -cpi $simname.cpt -plumed plumed$cnt.dat
done
