#!/bin/bash

# These are the SLURM commands
#SBATCH -o mcmc-%j-cleanup.out
#SBATCH -c 1
: '
This is the script that runs runs postprocessing on the MCMC. It waits for all instances of qsub.sh to finish. It gets submitted to SLURM by the qsub_wrapper.sh
'

# What is the config file
cfg=$1

# Excecute the config file to inherit the variables for the run
. $cfg

cd $workdir/$outname
cfg=$(basename $cfg)
Rscript $sourcedir/postprocess-mcmc.R $cfg 1> results.cleanup.out 2> results.cleanup.err
