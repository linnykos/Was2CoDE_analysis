#!/bin/bash

#SBATCH --account=biostat
#SBATCH --partition=all-12c128g
#SBATCH --time=48:00:00

#SBATCH --job-name=seaad-atac
#SBATCH --mem-per-cpu=100gb

#SBATCH --error=/home/users/%u/kzlinlab/projects/subject-de/out/kevin/Writeup6/slurm/slurm_%j.err
#SBATCH --output=/home/users/%u/kzlinlab/projects/subject-de/out/kevin/Writeup6/slurm/slurm_%j.out

source ~/miniconda3/etc/profile.d/conda.sh
conda activate scVI

conda env list
python --version
conda list

python Writeup6_seaad_atac-exploration.py > ../../../../../out/kevin/Writeup6/slurm/$SLURM_JOB_NAME-$SLURM_JOB_ID.txt