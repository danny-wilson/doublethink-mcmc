#!/bin/bash
: '
This is the script that that submits qsub_summarize_again.sh as an array job to SLURM (or SGE/UGE). Summarizing doublethink analysis again starts here, by starting the script like this:
./qsub_wrapper_summarize_again.sh <cfg file> <new burnin> <new thinning>
'

 timestamp=$(date +"%Y_%d_%m__%H_%M")

 # This is the config file
cfg=$1

# New burnin
new_burnin=$2

# New thinning
new_thinning=$3

# Run config file to get the variable
. $cfg
chain_number=$(($chain_numbers))
wdir="${workdir}/${outname}"

# Submit qsub-summarize-again.sh as an array job to SLURM for SGE or UGE uncomment the lines after
echo "sbatch --array 1-${chain_number} -o $workdir/$outname/ -J ${outname}_$timestamp qsub-summarize-again.sh $cfg $new_burnin $new_thinning"
output=$(sbatch --array 1-${chain_number} -J ${outname}_$timestamp qsub-summarize-again.sh $cfg $new_burnin $new_thinning)
jid="${output##* }"
# Run cleanup-again.sh that waits for all instances of qsub-summarize-again.sh to finish.
sbatch -J ${outname}_$timestamp --dependency=afterany:${jid} cleanup-again.sh ${cfg} ${new_burnin} ${new_thinning}

# This is for SGE or UGE. Uncomment these lines and comment the above SLURM commands
#qsub -N lr_mcmc_${timestamp} -t 1-${chain_numbers} qsub-summarize-again.sh ${cfg} ${new_burnin} ${new_thinning}
#qsub -N cleanup_${timestamp} -hold_jid cleanup-again.sh ${cfg} ${new_burnin} ${new_thinning}

