#!/bin/bash

#SBATCH --job-name=trim_ref
#SBATCH --time=48:00:00
#SBATCH -p normal
#SBATCH --nodes=1

module load ddocent

#this is an example sbatch script to run dDocentHPC on a slurm supercomputer

#this will trim the fq.gz files using the settings in config.4.all
#it is assumed that the directory you run this script from has the 
#fq.gz files
bash dDocentHPC.bash trimFQ config.4.all

#this will assemble the fq.gz files in the mkREF directory that was
#created by the trim function.  We have to cd to the correct directory
#then execute dDocentHPC
cd mkREF
bash ../dDocentHPC.bash mkREF ../config.4.all
