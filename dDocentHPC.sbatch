#!/bin/bash -l

#SBATCH --job-name=mkBAM
#SBATCH -o mkBAM-%j.out
#SBATCH --time=00:00:00
#SBATCH --cpus-per-task=40
#SBATCH --exclusive=user

enable_lmod
module load container_env ddocent/2.9.4

export SINGULARITY_BIND=/home/e1garcia
export PARALLEL_SHELL=/bin/bash

# crun.ddocent bwa-meme
# crun.ddocent bwa

FUNCTION=$1
CONFIG=$2

crun.ddocent bash /home/e1garcia/shotgun_PIRE/dDocentHPC/dDocentHPC.bash $FUNCTION $CONFIG
