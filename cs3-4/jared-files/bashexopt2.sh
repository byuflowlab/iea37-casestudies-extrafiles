#!/bin/bash

#SBATCH --time=00:01:00   # walltime
#SBATCH --ntasks=4   # number of processor cores (i.e. tasks)
#SBATCH --nodes=3   # number of nodes
#SBATCH --mem-per-cpu=1024M   # memory per CPU core
#SBATCH -J "test flowfarm distributed "   # job name
#SBATCH --qos=test

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
module load julia

julia example_opt_2_distributed.jl
