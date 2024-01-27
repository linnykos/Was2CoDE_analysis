#!/bin/bash

#SBATCH --account=biostat
#SBATCH --partition=all-12c128g
#SBATCH --time=48:00:00

#SBATCH --job-name=esvd
#SBATCH --mem-per-cpu=50gb

#SBATCH --error=/home/users/%u/kzlinlab/projects/subject-de/out/kevin/Writeup2/slurm/slurm_%j.err
#SBATCH --output=/home/users/%u/kzlinlab/projects/subject-de/out/kevin/Writeup2/slurm/slurm_%j.out

R CMD BATCH --no-save --no-restore Writeup2_microglia-pvm_esvd.R ../../../../../out/kevin/Writeup2/slurm/slurm_$SLURM_JOB_ID.Rout
