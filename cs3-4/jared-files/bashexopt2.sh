#!/bin/bash

#SBATCH --time=00:01:00   # walltime
#SBATCH --ntasks=4   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=1024M   # memory per CPU core
#SBATCH -J "test flowfarm distributed "   # job name
#SBATCH --qos=test

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
module load julia
MACHINEFILE="nodes.$SLURM_JOB_ID"
# Generate Machinefile for mpi such that hosts are in the same
#  order as if run via srun
#
srun -l /bin/hostname | sort -n | awk '{print $2}' > $MACHINEFILE

julia --machine-file MACHINEFILE example_opt_2_distributed.jl
