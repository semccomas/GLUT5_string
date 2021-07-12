## Here I am making xvg files for gate COM so I can generate a numpy array for string_analysis

#files are formatted as: md/$iteration/$bead/$swarm/traj_comp.xtc

min_iteration=1
max_iteration=200

mddir=../../string_sims/influx_BFRU_demyst_CV/md

## can use TPR for anything, the topol.tpr should be the SAME, therefore, you dont have to do this on TCB

TM1='com of resid 30 to 37'
TM7='com of resid 289 to 295'
TM4='com of resid 136 to 145'
TM10='com of resid 386 to 394' 

## from GLUT5_atomistic
#for traj in ../input_f/*xtc;
#do
#t=`basename $traj`
#t=${t%%.*}
#gmx pairdist -f $traj -s ../input_f/$t.tpr -ref "$TM1" -sel "$TM7" -xvg none -o extracellular/$t.EC.xvg

#gmx pairdist -f $traj -s ../input_f/$t.tpr -ref "$TM4" -sel "$TM10" -xvg none -o intracellular/$t.IC.xvg


for iteration in $(seq $min_iteration $max_iteration); do

	for bead in $(seq 1 32); do   #34 beads, but 0 and 33 are fixed, so no calcs

		for swarm in $(seq 0 31); do   #32 swarms per bead

			#ls $mddir/$iteration/$bead/s$swarm/*xtc
			echo $mddir/$iteration/$bead/s$swarm
			gmx pairdist -f $mddir/$iteration/$bead/s$swarm/traj_comp.xtc -s $mddir/244/1/s1/topol.tpr -ref "$TM1" -sel "$TM7" -xvg none -o $mddir/$iteration/$bead/s$swarm/EC_gate_COM.xvg 
			##gmx pairdist -f $mddir/$iteration/$bead/s$swarm/traj_comp.xtc -s $mddir/295/1/s1/topol.tpr -ref "$TM4" -sel "$TM10" -xvg none -o $mddir/$iteration/$bead/s$swarm/IC_gate_COM.xvg 

done
done
done

