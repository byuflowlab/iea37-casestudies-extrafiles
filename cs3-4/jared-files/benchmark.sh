#!/bin/bash -l
#SBATCH --ntasks=1441
#SBATCH --mem-per-cpu=1024M   # memory per CPU core
#SBATCH --time=01:00:00
#SBATCH -J "practice distributed "   # job name
#SBATCH --qos=test

module load julia

julia call_benchmark.jl