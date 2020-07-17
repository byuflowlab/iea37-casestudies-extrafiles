#!/bin/bash -l
#SBATCH --nodes=2
#SBATCH --ntasks=4
#SBATCH --mem-per-cpu=1024M   # memory per CPU core
#SBATCH --time=00:00:20
#SBATCH -J "practice distributed "   # job name
#SBATCH --qos=test

# source $HOME/julia-v0.6/julia-environment
# cd working/folder/of/your/choice
julia distributed_practice.jl