#!/bin/bash

#SBATCH --account=biostat
#SBATCH --partition=all-12c128g
#SBATCH --time=48:00:00

#SBATCH --job-name=scvi
#SBATCH --mem-per-cpu=50gb

#SBATCH --error=/home/users/%u/kzlinlab/projects/subject-de/out/kevin/Writeup3/slurm/slurm_%j.err
#SBATCH --output=/home/users/%u/kzlinlab/projects/subject-de/out/kevin/Writeup3/slurm/slurm_%j.out

conda init
conda activate scVI
python Writeup3_microglia-pvm_resilience_scvi.py > ../../../../../out/kevin/Writeup3/slurm/$SLURM_JOB_NAME-$SLURM_JOB_ID.txt