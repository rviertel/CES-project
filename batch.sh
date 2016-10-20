#!/bin/bash
#SBATCH -J CES_PROJECT
#SBATCH -o CES_PROJECT.o%j
#SBATCH -N 16 -n 16
#SBATCH -t 10:00:00
#SBATCH --mail-user=rdviertel@gmail.com
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --account=borisyuk
#SBATCH --partition=ember-freecycle

module load openmpi gsl

mpirun uq -s
