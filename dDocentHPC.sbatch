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

