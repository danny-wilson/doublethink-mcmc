#!/bin/bash
: '
This is the script that that submits qsub.sh as an array job to SLURM as well as cleanup.sh which waits for all instances of qsub.sh to finish. The doublethink analysis starts here, by starting the script like this:
./qsub_wrapper.sh <cfg file>
'

timestamp=$(date +"%Y_%d_%m__%H_%M")

 # This is the config file
cfg=$1

# Run config file to get the variable
. $cfg
chain_number=$(($chain_numbers))
wdir="${workdir}/${outname}"

# Create the working directory if it does not already exist
mkdir -p $wdir

# Remove holding files if they already exist
rm -f ${wdir}/hold_file
rm -f ${wdir}/hold_file1
rm -f ${wdir}/done_file

# Copy the config file into the work dir for tracing the variables of each run
cp $1 $wdir
cfg_basename=$(basename $cfg)

# Put this start string in the copy of the config file as R requires it
cd $wdir
startstring="#file:$cfg_basename\n"
sed -i "1s/^/${startstring}/" $cfg_basename

# Submit qsub.sh as an array job to SLURM
echo "sbatch --array 1-${chain_number} -J ${outname}_$timestamp $sourcedir/qsub.sh $cfg"
output=$(sbatch --array 1-${chain_number} -J ${outname}_$timestamp $sourcedir/qsub.sh $cfg)
jid="${output##* }"
# Run cleanup.sh that waits for all instances of qsub.sh to finish.
sbatch -J ${outname}_$timestamp --dependency=afterany:${jid} $sourcedir/cleanup.sh ${cfg}
