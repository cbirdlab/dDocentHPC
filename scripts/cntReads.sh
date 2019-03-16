#!/bin/bash

#SBATCH --job-name=cntREADS
#SBATCH --time=96:00:00
#SBATCH -p normal
#SBATCH --nodes=1

module load parallel

#count reads of any file ending in fq.gz in the present directory

ls *fq.gz | parallel "echo -n {}' ' && zgrep -c '^\+$' {}" | sort -g >> NumFqGzReads.txt
