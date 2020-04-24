#!/bin/sh --login
# -- Request number of processes to use
#PBS -l procs=13
# Maximize scheduling through HH:MM:SS (max 8 hour)
#PBS -l walltime=00:30:00
# -- Specify the job queue
#PBS -q batch
# -- Specify under which account a job should run
#PBS -A coastal
# -- Set the name of the job, or moab will default to STDIN
#PBS -N a70_CHI_ATM_WAV2OCNv2.0

# -- tells the batch system to remember all of your environmental variables
# --  #PBS -V

# -- Specify where to put stdout and stderr
#PBS -o job.${PBS_JOBID}.out
#PBS -e job.${PBS_JOBID}.err

# change directory to the working directory of the job
# Use the if clause so that this script stays portable

if [ x$PBS_O_WORKDIR != x ]; then
   cd $PBS_O_WORKDIR
fi

# -- load ENV variables
source  $appfolder/modulefiles/$machine/$modulefilename

mpirun -np 13 ./NEMS.x 



