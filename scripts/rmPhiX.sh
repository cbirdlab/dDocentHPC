#!/bin/bash

#run this script on your data after demultiplexing.
#phiX in your reads occurs due to index swapping
#this will give you an estimate of how much index swapping is occuring
#refer to the appendix entry on index swappingin the RAD Data Processing Manual

#This script assumes that the following files are in their typical directory: dDocent224tref_hpc6.bash, config.phix, phix genome

#SBATCH --job-name=rmPhiX
#SBATCH --time=96:00:00
#SBATCH -p cbirdq
#SBATCH --nodes=1
#SBATCH --mail-user=yourname@tamucc.edu
#SBATCH --mail-type=begin  # email me when the job starts
#SBATCH --mail-type=end    # email me when the job finish

module load ddocent/2.24
module load parallel

#Because of the way I've hacked dDocent, the phiX genome
#file should be labeled as reference.174.174.fasta
cutoffs="174.174"

echo " ";echo `date` "THIS SCRIPT REMOVES PHIX READS"

#copy the phiX genome to current directory and rename
cp /work/GenomicSamples/genomes/phiX174.fasta ./reference.$cutoffs.fasta

#change naming convention on fq.gz files for mapping to occur
rename .F. .R1. *.fq.gz
rename .R. .R2. *.fq.gz

echo " ";echo `date` "DdOCENT BEGAN MAPPING READS TO PHIX GENOME"

#map reads to phiX reference
bash /work/GenomicSamples/GCL/scripts/dDocent224tref_hpc6.bash /work/GenomicSamples/GCL/scripts/config.phix

echo " ";echo `date` "CONVERTING BAM TO SAM"

#convert bam to sam so it can be grepped
ls *.$cutoffs-RG.bam | sed 's/.bam//g' | parallel "samtools view {}.bam > {}.sam "

echo " ";echo `date` "COUNTING THE NUMBER OF PHIX READS"

#count the number of lines in the sam file, which is the number of reads that are phiX contamination
ls *.sam | sed 's/.sam//g' | parallel "echo -n {}' ' && wc -l {}.sam " | sort -g > NumPhixReads.txt

#grab list of phiX read names, these can be used to grep the fq.gz files
ls *.sam | sed 's/.sam//g' | parallel "grep -Eo '^\w+:\w+:\w+:\w+:\w+:\w+:\w+' {}.sam > {}.phiX "  

echo " ";echo `date` "PURIFYING fq.gz FILES OF PHIX"

#use list of phiX reads to purify the fq.gz files
ls *R1.fq.gz | sed 's/R1.fq.gz//g' | parallel "zgrep -x -v -A3 -f {}$cutoffs-RG.phiX {}R1.fq.gz > noPhiX-{}R1.fq.gz "
ls *R2.fq.gz | sed 's/R2.fq.gz//g' | parallel "zgrep -x -v -A3 -f {}$cutoffs-RG.phiX {}R2.fq.gz > noPhiX-{}R2.fq.gz "

#change naming convention on phiX free fq.gz files
rename .R1. .F. noPhiX*fq.gz
rename .R2. .R. noPhiX*fq.gz

#remove the sam files
rm *.sam

echo " ";echo `date` "SCRIPT COMPLETED SUCCESSFULLY!"

