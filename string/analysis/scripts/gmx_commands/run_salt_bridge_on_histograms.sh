#!/bin/bash
#SBATCH -p tcb
#SBATCH -J conf
#SBATCH -t 23:00:00
#SBATCH -N 1 --ntasks-per-node=4
#SBATCH -e error.err -o output.out


module unload gromacs
module load gromacs/2020.2

condition=influx_BFRU_gate_CV

for i in $(seq 298 299); do

	indir=../confout_files/FES_grids_confouts/$condition/histogram_$i
	outdir=../textfiles_out/salt_bridge_GMX/$condition
	tprdir=../confout_files/tpr_files

	if [ ! -d "$outdir" ]
		then
		mkdir $outdir
	else
		echo 'directory exists: ' $outdir
	fi

	if [ -f $indir/FES_grid_all.xtc ]; then
		mkdir $outdir/histogram_$i

		gmx pairdist -f $indir/FES_grid_all.xtc -s $tprdir/$condition.wholesys.tpr -ref "resid 145" -sel "resid 91" -xvg none -o $outdir/histogram_$i/E145-R91.xvg
		gmx pairdist -f $indir/FES_grid_all.xtc -s $tprdir/$condition.wholesys.tpr -ref "resid 145" -sel "resid 401" -xvg none -o $outdir/histogram_$i/E145-R401.xvg
		gmx pairdist -f $indir/FES_grid_all.xtc -s $tprdir/$condition.wholesys.tpr -ref "resid 394" -sel "resid 334" -xvg none -o $outdir/histogram_$i/E394-R334.xvg
		gmx pairdist -f $indir/FES_grid_all.xtc -s $tprdir/$condition.wholesys.tpr -ref "resid 394" -sel "resid 152" -xvg none -o $outdir/histogram_$i/E394-R152.xvg
		gmx pairdist -f $indir/FES_grid_all.xtc -s $tprdir/$condition.wholesys.tpr -ref "resid 330" -sel "resid 334" -xvg none -o $outdir/histogram_$i/E330-R334.xvg

		rm $outdir/#*

	else
		echo "no confouts found in histo $i"
	fi

done

