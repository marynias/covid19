#!/bin/bash

#SBATCH --job-name=script
#SBATCH --partition=hmem
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=14
#SBATCH --time=08:00:00
#SBATCH --mem=20G
#SBATCH -e slurm-%j.err
#SBATCH -o slurm-%j.out


srun Rscript --vanilla $input
