#!/bin/bash -l
#SBATCH --nodes=10
#SBATCH --ntasks=11
#SBATCH --mem-per-cpu=1024M   # memory per CPU core
#SBATCH --time=00:10:00
#SBATCH -J "practice distributed "   # job name
#SBATCH --qos=test

# source $HOME/julia-v0.6/julia-environment
# cd working/folder/of/your/choice
julia example_opt_snopt_raytrace_4_ieacs4.jl