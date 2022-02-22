#!/bin/bash

#SBATCH --job-name=cntPhiX
#SBATCH -p cbirdq
#SBATCH --nodes=1

module load ddocent/2.24
module load parallel

#this script will generate files that quantify phix contamination then leave no trace

#This script depends upon dDocent224tref_hpc6.bash, config.phix, phix genome


##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#set this for your phiX genome.  Because of the way I've hacked dDocent, the phiX genome
#file should be renamed as reference.174.174.fasta.
cutoffs="174.174"
##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

echo " ";echo `date` "THIS SCRIPT QUANTIFIES PHIX CONTAMINATION"

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
ls *.sam | sed 's/.sam//g' | parallel "echo -n {}' ' && wc -l {}.sam && echo -n ' ' && head -n 1  " | sort -g > NumPhixReads.txt

echo " ";echo `date` "COUNTING THE NUMBER OF READS"

#count the total number of reads per file
ls *R1.fq.gz | sed 's/.fq.gz//g' | parallel "echo -n {}' ' && zgrep -c '^\+$' {}.fq.gz " | sort -g > NumTotalReads.txt

echo " ";echo `date` "MAKING FILES WITH NAMES OF PHIX CONTAMINATED READS"

#grab list of phiX read names, these can be used to grep the fq.gz files
ls *.sam | sed 's/.sam//g' | parallel "grep -Eo '^\w+-?\w?+?:\w+:\w+-?\w?+?:\w+:\w+:\w+:\w+' {}.sam > {}.phiX "

echo " ";echo `date` "CLEANUP ON AISLE 7"

rm *.sam
rm *.bam
rm *.bam.bai
rm reference.$cutoffs.*
rm mapped.$cutoffs.bed
rm bamlist.$cutoffs.list
rm namelist*


#change naming convention on phiX free fq.gz files
rename .R1. .F. *fq.gz
rename .R2. .R. *fq.gz


echo " ";echo `date` "SCRIPT COMPLETED SUCCESSFULLY!"

