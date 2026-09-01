#!/bin/bash
#SBATCH --account=beenome100
#SBATCH --output=PCA_%j.out
#SBATCH --error=PCA_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=gzayas97@ufl.edu
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=20gb
