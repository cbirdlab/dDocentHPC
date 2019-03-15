#!/bin/bash

#SBATCH --job-name=fqc
#SBATCH --time=96:00:00
#SBATCH -p normal
#SBATCH --nodes=1
#SBATCH --mail-user=cbird@tamucc.edu
#SBATCH --mail-type=begin  # email me when the job starts
#SBATCH --mail-type=end    # email me when the job finish

module load fastqc
module load parallel
module load multiqc

#run fastqc in parallel
ls *fq.gz | parallel "fastqc {}"

#make a directory called FASTQC and move the FASTQC files into the directory
mkdir FASTQC
ls *fastqc.html | parallel "mv {} FASTQC"
ls *fastqc.zip | parallel "mv {} FASTQC"

multiqc FASTQC/
