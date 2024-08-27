#!/bin/bash

#SBATCH --account=biostat
#SBATCH --partition=all-20c96g
#SBATCH --time=48:00:00

#SBATCH --job-name=deseq2
#SBATCH --mem-per-cpu=50gb

#SBATCH --error=/home/users/%u/kzlinlab/projects/subject-de/out/tati/Writeup4/slurm/slurm_%j.err
#SBATCH --output=/home/users/%u/kzlinlab/projects/subject-de/out/tati/Writeup4/slurm/slurm_%j.out

R CMD BATCH --no-save --no-restore Writeup4_SEA-AD_Pseudobulk-DEseq2.R ../../../../../out/tati/Writeup4/slurm/$SLURM_JOB_NAME-$SLURM_JOB_ID.Rout
