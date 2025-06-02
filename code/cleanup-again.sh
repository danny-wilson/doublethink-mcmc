#!/bin/bash

# These are the SLURM commands
#SBATCH -o %j_cleanup_again.out
#SBATCH -c 1
: '
This is the script that runs runs postprocessing-again on the MCMC. It waits for all instances of qsub-summarize-again.sh to finish. It gets submitted to SLURM by the qsub_wrapper_summarize_again.sh
'

# What is the config file
cfg=$1

# Burnin
new_burnin=$2

# Thinning
new_thinning=$3

# Excecute the config file to inherit the variables for the run
. $cfg


# Run the postprocessing and dump the running log in results.cleanup.out
module purge
module load R/3.6.2-foss-2019b

cd $workdir/$outname
cfg=$(basename $cfg)
Rscript $sourcedir/postprocess-mcmc-again.R $cfg $new_burnin $new_thinning 1> results.cleanup-again.out 2> results.cleanup-again.err
