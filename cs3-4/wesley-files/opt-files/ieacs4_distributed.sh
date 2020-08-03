#!/bin/bash/ -l
#SBATCH --nodes=10
#SBATCH --ntasks=11
#SBATCH --mem-per-cpu=1024M  # memory per CPU core
#SBATCH --time=00:15:00
#SBATCH -J "practice distributed ieacs4"   # job name
#SBATCH --qos=test

julia example_opt_4_ieacs4_distributed.jl
