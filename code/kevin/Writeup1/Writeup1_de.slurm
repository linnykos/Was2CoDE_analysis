#!/bin/bash

#SBATCH --job-name=de
#SBATCH --account=biostat
#SBATCH --partition=all-12c128g
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=25gb

#SBATCH --error=/home/users/%u/kzlinlab/projects/subject-de/out/kevin/Writeup1/slurm/slurm_%j.err
#SBATCH --output=/home/users/%u/kzlinlab/projects/subject-de/out/kevin/Writeup1/slurm/slurm_%j.out

R CMD BATCH --no-save --no-restore Writeup1_de.R ../../../../../out/kevin/Writeup1/slurm/$SLURM_JOB_NAME-$SLURM_JOB_ID.Rout