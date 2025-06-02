#!/bin/bash

# Based on qsub.sh

# These are the SLURM commands
#SBATCH -J sumagain
#SBATCH -o mcmc-%j.out
#SBATCH -e mcmc-%j.err
#SBATCH -c 1
: '
This is the script that runs summarize-again.R. It gets submitted as a single instance of an array job to SLURM by the qsub_wrapper_summarize_again.sh
'

# Which chain is this
chain_id=${SLURM_ARRAY_TASK_ID}

# What is the config file
cfg=$1

# Burnin
new_burnin=$2

# Thinning
new_thinning=$3

# Execute the config file to inherit the variables for the run
. $cfg

# Output some descriptors of the run for error handling
echo "------------------------------------------------"
echo "SGE Job ID: "$JOB_ID
echo "SGE Task ID: "$SLURM_ARRAY_TASK_ID
echo "Run on host: "`hostname`
echo "Username: "$USER
echo "Started at: "`date`
echo "Output file: ${wdir}/${outname}_chain_${chain_id}.out"
echo "------------------------------------------------"
echo "name: ${outname}"
echo "chain number: ${chain_id}"
echo "Number iterations per chain : "$niter
echo "Script dir: "$sourcedir
echo "Work dir: "$workdir
echo "new_thinning: "$new_thinning
echo "new_burnin: "$new_burnin
echo "------------------------------------------------"

# Add R from the modules. If you don't have the moduling system remove these lines. 
module purge
module add R/3.6.2-foss-2019b

# Define work dir
wdir="${workdir}/${outname}"
mkdir -p $wdir

# Get into the working dir
cd $wdir

# Get the local cfg file
cfg=$(basename $cfg)

# Put this start string in the config file as R requires it
echo "NOW STARTING"
# Run summarize-again.R
Rscript $sourcedir/summarize-again.R $cfg $chain_id $new_burnin $new_thinning 2> ${wdir}/${outname}_chain_${chain_id}_burnin_${new_burnin}_thinning_${new_thinning}.err 1> ${wdir}/${outname}_chain_${chain_id}_burnin_${new_burnin}_thinning_${new_thinning}.out
