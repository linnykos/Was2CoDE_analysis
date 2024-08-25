#!/bin/bash

#SBATCH --account=biostat
#SBATCH --partition=all-12c128g
#SBATCH --time=72:00:00

#SBATCH --job-name=seaad-was2
#SBATCH --mem-per-cpu=100gb

#SBATCH --error=/home/users/%u/kzlinlab/projects/subject-de/out/tati/Writeup6/slurm/slurm_%j.err
#SBATCH --output=/home/users/%u/kzlinlab/projects/subject-de/out/tati/Writeup6/slurm/slurm_%j.out

R CMD BATCH --no-save --no-restore Writeup6_SEA-AD_Was2.R ../../../../../out/tati/Writeup6/slurm/$SLURM_JOB_NAME-$SLURM_JOB_ID.Rout
