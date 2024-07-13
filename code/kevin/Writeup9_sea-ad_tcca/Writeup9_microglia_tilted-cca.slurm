#!/bin/bash

#SBATCH --account=biostat
#SBATCH --partition=all-20c96g
#SBATCH --time=48:00:00

#SBATCH --job-name=tcca
#SBATCH --mem-per-cpu=80gb

#SBATCH --error=/home/users/%u/kzlinlab/projects/subject-de/out/kevin/Writeup9/slurm/slurm_%j.err
#SBATCH --output=/home/users/%u/kzlinlab/projects/subject-de/out/kevin/Writeup9/slurm/slurm_%j.out

R CMD BATCH --no-save --no-restore Writeup9_microglia_tilted-cca.R ../../../../../out/kevin/Writeup9/slurm/$SLURM_JOB_NAME-$SLURM_JOB_ID.Rout
