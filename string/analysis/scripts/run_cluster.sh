#!/bin/bash
#SBATCH -p tcb
#SBATCH -J conf
#SBATCH -t 23:00:00
#SBATCH -N 1 --ntasks-per-node=4
#SBATCH -e error.err -o output.out


module unload gromacs
module load gromacs/2020.2

for i in $(seq 298 299); do

	indir=../confout_files/FES_grids_confouts/influx_BFRU_gate_CV/histogram_$i
	outdir=../textfiles_out/gmx_cluster_files/histogram_$i
	ndxfile=../confout_files/tpr_files/index_BFRU.ndx

	if [ ! -d "$outdir" ]
		then
		mkdir $outdir
	else
		echo 'directory exists: ' $outdir
	fi

	if [ -f $indir/FES_grid_all.protonly.cluster.fit.xtc ]; then
		
		gmx cluster -f $indir/FES_grid_all.protonly.cluster.fit.xtc -s $indir/FES_grid_all.start.protonly.pdb -nofit -cl $outdir/clusters.pdb -n $ndxfile -clid $outdir/clust-id.xvg -g $outdir/cluster.log  -xvg none -method jarvis-patrick -cutoff 0.08 -dist $outdir/rmsd-dist.xvg -o $outdir/rmsd-clust.xpm << EOF
13
18
EOF
	
		rm $outdir/#*

	else
		echo "no confouts found in histo $i"
	fi

done
