#!/bin/bash

# Use the Bash Shell
#$ -S /bin/bash

# Put a unique job name here
# so that you can monitor it in the queue
#$ -N DeerPop

# Notifying user by email at the beginning,end,abort,suspensions of the job run
#$ -M kblack@clarkson.edu
#$ -m eas

# Tell GE to run the job from the current working directory
#$ -cwd 

# Uncomment to pass all current environment variables
#$ -V
# Uncomment to pass a single environment variable
# #$ -v VAR

# Redirecting standard output / error
#$ -o qoutput
#$ -e qerrors

# The max walltime for this job is 31 minutes
#$ -l h_rt=600:31:00

# Parallel execution request for MPICH. Set your number of processors here.
#$ -pe mpich 2

local_dir=$(pwd)
echo "This job was started in ${local_dir}"
echo "The start time was $(date)"
echo "The job id is $JOB_ID"
echo "Got $NSLOTS processors."
echo "Temp directory: $TMPDIR"
echo "Machines:"
cat $TMPDIR/machines
echo ${PE}

# Executing the program
/share/apps/bin/mpirun -np $NSLOTS -machinefile $TMPDIR/machines ./mpiFileOperation

# Submit using qsub -pe orte 16 submit.sh

