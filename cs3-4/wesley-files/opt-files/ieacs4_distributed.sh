#!/bin/bash

#SBATCH --ntasks=12
#SBATCH --mem-per-cpu=2048M  # memory per CPU core
#SBATCH --time=08:00:00
#SBATCH -J 'ieacs4 - WEC with reduced wind rose 103'
#SBATCH --mail-user=wesleyjholt@gmail.com   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
# #SBATCH --array=1
# #SBATCH --qos=test

LAYOUT_NUMBER=103
TOL=3.8e-6

module load julia
julia example_opt_4_ieacs4_WEC_discrete.jl ${LAYOUT_NUMBER} ${TOL}
