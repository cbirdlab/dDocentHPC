#!/bin/bash

#SBATCH --job-name=trimFQ-mkVCF
#SBATCH -o trimFQ-mkREF-mkBAM-fltrBAM-mkVCF-%j.out
#SBATCH --time=48:00:00
#SBATCH -p normal
#SBATCH --nodes=1

module load ddocent

#this is an example sbatch script to run dDocentHPC on a slurm supercomputer
#it is assumed that you have copied dDocentHPC.bash and config.4.all to your project directory
#and your *F.fq.gz and *R.fq.gz files are in the project directory

#while this will run all of dDocentHPC's functions, it is recommended to run each function independently
#and evaluate the output before moving to the next step

#this will trim the fq.gz files using the settings in config.4.all
#it is assumed that the directory you run this script from has the 
#fq.gz files
bash dDocentHPC.bash trimFQ config.4.all

#or run both trims at same time
#srun --nodes=1 bash dDocentHPC.bash trimFQref config.4.all &
#srun --nodes=1 bash dDocentHPC.bash trimFQmap config.4.all &
#wait

#this will assemble the fq.gz files in the mkREF and mkBAM directories that were
#created by the trim function.  We have to cd to the mkREF directory
#then execute dDocentHPC mkREF
cd mkREF
bash ../dDocentHPC.bash mkREF ../config.4.all

#this will use dDocent mkBAM to map the reads in the mkBAM directory to the reference*fasta created by 
#dDocent mkREF.  We have to change directories, copy the reference*fasta to mkBAM, then run dDocent mkBAM 
cd ../mkBAM
cp ../mkREF/reference.*.fasta .
bash ../dDocentHPC.bash mkBAM ../config.4.all

#this will use dDocent fltrBAM to filter the BAM files
bash ../dDocentHPC.bash fltrBAM ../config.4.all

#this will use freebayes to genotype the bam files and make a VCF
bash ../dDocentHPC.bash mkVCF ../config.4.all

