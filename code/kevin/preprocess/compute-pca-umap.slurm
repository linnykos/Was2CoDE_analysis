#!/bin/bash

#SBATCH --account=biostat
#SBATCH --partition=all-12c128g
#SBATCH --time=48:00:00

#SBATCH --job-name=preprocess
#SBATCH --mem-per-cpu=25gb

#SBATCH --error=/home/users/kzlin/lab/projects/subject-de/out/kevin/preprocess/slurm/slurm_%j.err
#SBATCH --output=/home/users/kzlin/lab/projects/subject-de/out/kevin/preprocess/slurm/slurm_%j.out

R CMD BATCH --no-save --no-restore compute-pca-umap.R ../../../../../out/kevin/preprocess/slurm/slurm_$SLURM_JOB_ID.Rout