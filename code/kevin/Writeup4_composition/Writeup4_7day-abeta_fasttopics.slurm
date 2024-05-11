#!/bin/bash

#SBATCH --account=biostat
#SBATCH --partition=all-12c128g
#SBATCH --time=48:00:00

#SBATCH --job-name=aq_ft
#SBATCH --mem-per-cpu=50gb

#SBATCH --error=/home/users/%u/kzlinlab/projects/subject-de/out/kevin/Writeup4/slurm/slurm_%j.err
#SBATCH --output=/home/users/%u/kzlinlab/projects/subject-de/out/kevin/Writeup4/slurm/slurm_%j.out

R CMD BATCH --no-save --no-restore Writeup4_7day-abeta_fasttopics.R ../../../../../out/kevin/Writeup4/slurm/$SLURM_JOB_NAME-$SLURM_JOB_ID.Rout
