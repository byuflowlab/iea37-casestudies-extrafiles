#!/bin/bash -l
#SBATCH --nodes=2
#SBATCH --ntasks=4
#SBATCH --ntasks-per-node=2
#SBATCH --time=00:10:00
#SBATCH --output=myoutput.log
#SBATCH --job-name=my-julia-job

# source $HOME/julia-v0.6/julia-environment
# cd working/folder/of/your/choice
julia distributed_practice.jl