#!/bin/bash

#SBATCH --account=biostat
#SBATCH --partition=all-12c128g
#SBATCH --time=48:00:00

#SBATCH --job-name=nebula
#SBATCH --mem-per-cpu=50gb

#SBATCH --error=/home/users/%u/kzlinlab/projects/subject-de/out/tati/Writeup6/slurm/slurm_%j.err
#SBATCH --output=/home/users/%u/kzlinlab/projects/subject-de/out/tati/Writeup6/slurm/slurm_%j.out

R CMD BATCH --no-save --no-restore Writeup6_SEA-AD_NEBULA.R ../../../../../out/tati/Writeup6/slurm/$SLURM_JOB_NAME-$SLURM_JOB_ID.Rout
