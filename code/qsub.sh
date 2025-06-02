#!/bin/bash

# These are the SLURM commands
#SBATCH -J lr-mcmc
#SBATCH -o mcmc-%j.out
#SBATCH -e mcmc-%j.err
#SBATCH -c 1
: '
This is the script that runs run-mcmc.R. It gets submitted as a single instance of an array job to SLURM by the qsub_wrapper.sh
'

# Which chain is this
chain_id=${SLURM_ARRAY_TASK_ID}

# What is the config file
cfg=$1

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
echo "------------------------------------------------"

# Define work dir
wdir="${workdir}/${outname}"

# Basename of the config file
cfg_basename=$(basename $cfg)

# Get into the working dir
cd $wdir

echo "NOW STARTING"
# Run run-mcmc.R
Rscript $sourcedir/run-mcmc.R $cfg_basename $chain_id 2> ${wdir}/${outname}_chain_${chain_id}.err 1> ${wdir}/${outname}_chain_${chain_id}.out
