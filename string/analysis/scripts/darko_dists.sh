#!/bin/sh

###$1 = condition
###$2 = iteration # 
## $3 = which CV group
## $4 = contact 1 resid
## $5 = contact 1 atom name
## $6 = contact 2 resid
## $7 contact 2 atom name

indir="../confout_files/measure_per_iteration"
tprdir="../confout_files/tpr_files"

## Outward Open
outdir='../textfiles_out/darko_dists_GMX'

gmx pairdist -f $indir/$1.$2.string.pdb -o $outdir/$3/$1.$2.$4-$6.xvg -ref "resid $4 and atomname $5" -sel "resid $6 and atomname $7" -xvg none -s $tprdir/$1.protonly.tpr


gmx pairdist -f ../confout_files/models_conf/OUT-IN.pdb -o $outdir/$3/OUT-IN.$4-$6.xvg -ref "resid $4 and atomname $5" -sel "resid $6 and atomname $7" -xvg none -s $tprdir/influx_apo_gate_CV.protonly.tpr
