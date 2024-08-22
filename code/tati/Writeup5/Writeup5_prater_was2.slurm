#!/bin/bash

#SBATCH --account=biostat
#SBATCH --partition=all-12c128g
#SBATCH --time=48:00:00

#SBATCH --job-name=katie-was2
#SBATCH --mem-per-cpu=100gb

#SBATCH --error=/home/users/%u/kzlinlab/projects/subject-de/out/tati/Writeup5/slurm/slurm_%j.err
#SBATCH --output=/home/users/%u/kzlinlab/projects/subject-de/out/tati/Writeup5/slurm/slurm_%j.out

R CMD BATCH --no-save --no-restore Writeup5_prater_Was2.R ../../../../../out/tati/Writeup5/slurm/$SLURM_JOB_NAME-$SLURM_JOB_ID.Rout